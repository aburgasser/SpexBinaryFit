; -----------------------------------------------------------------------------------------------------
 PRO SPEXBINARYFIT, $
;
; Purpose: Fit binary templates of late M, L and T dwarfs to a near-infrared spectrum 
;  (preferably spex prism data); output best fits and mean component properties
; 
;  See Burgasser et al. (2010) for details on operation of program
;
;  See spexbinaryfit_help.txt file for more information on usage of this program 
;
;  EXTERNAL FUNCTIONS REQUIRED (to be integrated)
; 	append
;	plotps
;	pop
;	readarr
;	rndoff
;	specindex
;  NOTE: REQUIRES GSFC LIBRARY
;
;  HISTORY
;  CREATED: 19 Nov 2009 by Adam J. Burgasser based on Burgasser et al. (2010)
;  UPDATES:
;  2010 Dec 9 (AJB) 
;	improved output with spexbinaryfit_printout, made template construction easier
;  2011 Mar 25 (AJB)
;	updated with redefined template structure
;  2012 June 12 (AJB)
;	fixed plotting error
;  2013 Jan 27 (AJB)
;	updated help file
;	changed single template update to /regenerate and updated templates
;	fixed scaling problem in fitting
;	added additional reporting options	
;	streamlined operation by running singles first
;	added /noyoung and /nosubdwarfs for rejecting some singles 
;  2014 Feb 11 (AJB)
;	finished off output structure output
;	added bestfitsingle/bestfitbinary outputs to store [[wave],[flux],[noise]] of 
;		best single and binary templates 
;	replaced CONSTANT with WEIGHT (constant weights are now default)
;	added name variable for template rejection
; 2014 Oct 27 (AJB)
;	added showdiff keyword
;	eliminated binary template structures and forced direct modeling at run
;	introduced new spexbinaryfit_mask_templates function to reject before comparison
;	introduced new spexbinaryfit_make_binaries function
;	changed baseline plot colors to color blind friendly
;	fixed rejection procedure to actually reject
;	added Dupuy & Liu 2012 MKO JHK relations
;
; TO DO
; - integrate external functions (see above)
; - add parameter input file to set presets
; - output template information file (in spexbinaryfit_library_singles)
; -----------------------------------------------------------------------------------------------------

prismfile, results, folder=folder, datafolder=datafolder, nbest=nbest, regenerate=regenerate, reset=reset, bright=bright, faint=faint, liu06=liu06, bur07=bur07, loo08=loo08, fah12=fah12, dup12=dup12, jmag=jmag, hmag=hmag, kmag=kmag, single=single, smooth=smooth, inset=inset, prefix=prefix, reject=reject, deltamag=deltamag, e_deltamag=e_deltamag, deltafilter=deltafilter, sigma=sigma, nirspt=nirspt, optspt=optspt, spexspt=spexspt, note=note, pcolors=pcolors, lcolors=lcolors, fitregions=fitregions, yrange=yrange, help=help, xrange=xrange, psptrange=psptrange, ssptrange=ssptrange, difference=difference, short=short, extras=extras, right=right, nonoise=nonoise, nfact=nfact, full=full, mko=mko, nicmos=nicmos, nirc2=nirc2, wfc3=wfc3, noyoung=noyoung, nosubdwarfs=nosubdwarfs, noblue=noblue, nored=nored, onlyyoung=onlyyoung, onlysubdwarfs=onlysubdwarfs, onlyblue=onlyblue, onlyred=onlyred, nocruz=nocruz, noplot=noplot, silent=silent, outputstr=outputstr, bestfitsingle=bestfitsingle, bestfitbinary=bestfitbinary, weight=weight, name=name, latex=latex, showdiff=showdiff

; ************************************************
; SETUP: YOU SHOULD ONLY NEED TO CHANGE THE FOLLOWING
; ************************************************
pfold = '/Users/adam/idl/spexbinaryfit/'		; folder containing these program files
dfold = '/Users/adam/spectra/spex/prism/'		; folder containing your spex data
; ************************************************


; HELP FILE 
if keyword_set(help) then begin
 f=''
 openr, unit, pfold+'spexbinaryfit_help.txt', /get_lun
 while (eof(unit) eq 0) do begin
  readf, unit, f
  print, f
 endwhile
 goto, FINISH
endif

; PROGRAM DEFAULTS
tb='	'
suffix = ''
calline = ''
sptstring = append('M'+strtrim(string(indgen(10)),2),$
 append('L'+strtrim(string(indgen(10)),2),'T'+strtrim(string(indgen(10)),2)))
templatesingfile = pfold+'/templates/templates_singles.dat'
fitregions0 = [[0.95,1.35],[1.45,1.8],[2.,2.35]]

time = systime(1)


; PARAMETER CHECK
if (n_elements(folder) eq 0) then folder = pfold+'fits/'
folder = folder+'/'
if (file_search(folder) eq '') then spawn, 'mkdir '+folder
if (n_elements(datafolder) eq 0) then datafolder = dfold 
if (n_elements(nbest) eq 0) then nbest=10
if (n_elements(prefix) eq 0) then prefix=''
if (n_elements(name) eq 0) then name=''
if (n_elements(note) eq 0) then note=name
if (n_elements(fitregions) eq 0) then fitregions=fitregions0
if (n_elements(normregions) eq 0) then normregions=fitregions0
if (n_elements(xrange) eq 0) then xrange=0
if (keyword_set(weight)) then statstr='G!Dk!N' else statstr = '!9c!3!U2!N'
if (keyword_set(difference) eq 1) then statstr = '!9s!3!U2!N'
if (n_elements(psptrange) lt 2) then psptrange=[15.,40.]
if (n_elements(ssptrange) lt 2) then ssptrange=[15.,40.]
if (n_elements(pcolors) eq 0) then pcolors=[0,2,4,6,0,0]
if (n_elements(lcolors) eq 0) then lcolors=[0,2,4,0,0,0]
if (n_elements(nfact) eq 0) then nfact=0.05		; default noise scaling of 5%
plotfitrange = fitregions
;normregions = fitregions

; REGENERATE SINGLE TEMPLATE LIBRARY IF NEEDED
; REQUIRES ACCESS TO ENTIRE LIBRARY
; fix this to also generate binary templates
print, templatesingfile
if (file_search(templatesingfile) eq '' or keyword_set(regenerate)) then begin
 spexbinaryfit_generate_library, dbfile=dbfile, datafolder=datafolder, pfold=pfold
endif

 ; PROCESS INPUT DATA
DATA: if (n_elements(prismfile) eq 0) then goto, FINISH
sz = size(prismfile)
case 1 of

; NORMAL FILE NAME
 sz(n_elements(sz)-2) eq 7: begin
  f = file_search(prismfile)
  if (f eq '') then begin
   f = file_search(datafolder+prismfile)
   if (f(0) eq '') then begin
    print, 'SPEXBINARYFIT: Data file '+prismfile+' not found locally or in '+datafolder
    goto, FINISH
   endif else prismfile = datafolder+prismfile
  endif
  spec = spexbinaryfit_readspectrum(prismfile)
  lam_raw = spec(*,0)
  flx_raw = spec(*,1)  
  ns_raw = spec(*,2)
 end

; INPUT ARRAY
 (sz(n_elements(sz)-2) eq 4 or sz(n_elements(sz)-2) eq 5) and sz(0) ge 2 : begin
  lam_raw = prismfile(*,0)
  flx_raw = prismfile(*,1)  
  if (sz(2) gt 2) then ns_raw = prismfile(*,2) else ns_raw = flx_raw*0.
 end
 
; STRUCTURE
 sz(n_elements(sz)-2) eq 8: begin
  names = tag_names(prismfile)
  w = where(strpos(names,'LAM') gt -1,cnt)
  if cnt gt 0 then lam_raw = prismfile.(w(0))
  w = where(strpos(names,'WAVE') gt -1,cnt)
  if cnt gt 0 then lam_raw = prismfile.(w(0))
  w = where(strpos(names,'FLX') gt -1,cnt)
  if cnt gt 0 then flx_raw = prismfile.(w(0))
  w = where(strpos(names,'FLUX') gt -1,cnt)
  if cnt gt 0 then flx_raw = prismfile.(w(0))
  w = where(strpos(names,'FLAM') gt -1,cnt)
  if cnt gt 0 then flx_raw = prismfile.(w(0))
  w = where(strpos(names,'NOISE') gt -1,cnt)
  if cnt gt 0 then ns_raw = prismfile.(w(0))
  w = where(strpos(names,'NS') gt -1,cnt)
  if cnt gt 0 then ns_raw = prismfile.(w(0))
  w = where(strpos(names,'ERR') gt -1,cnt)
  if cnt gt 0 then ns_raw = prismfile.(w(0))
  w = where(strpos(names,'UNC') gt -1,cnt)
  if cnt gt 0 then ns_raw = prismfile.(w(0))
  if (n_elements(lam_raw) eq 0) then begin
    print, 'SPEXBINARYFIT: could not recognize wavelength array in input structure file; please use LAM or WAVE'
    goto, FINISH
  endif
  if (n_elements(flx_raw) eq 0) then begin
    print, 'SPEXBINARYFIT: could not recognize flux array in input structure file; please use FLUX or FLAM'
    goto, FINISH
  endif
  if (n_elements(ns_raw) eq 0) then ns_raw = flx_raw*0.
 end 

 else: begin
   print, 'SPEXBINARYFIT: could not recognize input; please use filename, data array or structure'
   goto, FINISH
 end
endcase

; NOISE - if not present, or S/N <= 1 then assume 5% flux 
if (total(ns_raw) eq 0. or median(ns_raw/flx_raw) ge 1.) then nonoise=1
if (keyword_set(nonoise)) then begin
 if (n_elements(nfact) eq 0) then nfact = 0.05
 ns_raw = flx_raw*nfact
endif

; if just computing difference, set noise = 1
if (keyword_set(difference)) then ns_raw=ns_raw*0.+1.

; SMOOTH DATA IF DESIRED (not prism mode)
if (keyword_set(smooth)) then begin
 deres2, lam_raw, flx_raw, 120, flx_sm
 deres2, lam_raw, ns_raw, 120, ns_sm
endif else begin
 flx_sm = flx_raw
 ns_sm = ns_raw
endelse

; NORMALIZE - MAY NOT WANT TO DO THIS IN THE FUTURE
nflag = lam_raw*0
for i=0,n_elements(normregions(0,*))-1 do begin
 w = where(lam_raw ge normregions(0,i) and lam_raw le normregions(1,i),cnt)
 if (cnt gt 0) then nflag(w) = 1
endfor
wnorm = where(nflag eq 1,cnt)
scl = max(flx_sm(wnorm))
flx_sm = flx_sm/scl
ns_sm = ns_sm/scl

; -----------------------------------------------
; PREP TEMPLATES
; -----------------------------------------------

; READ IN SINGLE TEMPLATES
restore, file= templatesingfile
stdlam = db_singles.wavelength

; INTERPOLATE DATA ONTO STANDARD WAVELENGTH SCALE
flx = interpol(flx_sm,lam_raw,stdlam)
ns = interpol(ns_sm,lam_raw,stdlam)

; SELECT SPECTRAL TYPE FORM TO USE
case 1 of
 keyword_set(optspt): begin
  refspt = db_singles.opt_type
  refsptn = db_singles.osptn
 end
 keyword_set(nirspt): begin
  refspt = db_singles.nir_type
  refsptn = db_singles.nsptn
 end
 keyword_set(spexspt): begin
  refspt = db_singles.spex_type
  refsptn = db_singles.spex_sptn
 end
 else: begin
  refspt = db_singles.spts
  refsptn = db_singles.sptn
 end
endcase
db_original = create_struct(db_singles,'refspt',refspt,'refsptn',refsptn)

; -----------------------------------------------
; TEMPLATE REJECTION 
; -----------------------------------------------

singlerejflag = intarr(db_original.ntemplates)*0		; same as primaryrejflag
secondaryrejflag = intarr(db_original.ntemplates)*0

; REJECT TEMPLATES IF NAME AGREES WITH INPUT NAME
wrej = where(name ne '' and (strpos(db_original.shortname,name) gt 0 or $
  strpos(strupcase(db_original.name),strupcase(name)) gt 0),cntrej)
if (cntrej gt 0) then begin
 singlerejflag(wrej) = 1
 secondaryrejflag(wrej) = 1
endif

; REJECT TEMPLATES BASED ON REJECT INPUT SHORTNAME
if (n_elements(reject) ne 0) then begin
 if (n_elements(reject) eq 1) then reject = [reject]
 for j=0,n_elements(reject)-1 do begin
  wrej = where(db_original.shortname eq reject(j) or db_original.name eq reject(j),cntrej)
  if (cntrej gt 0) then begin
   singlerejflag(wrej) = 1
   secondaryrejflag(wrej) = 1
 endif
endfor
endif

; REJECT TEMPLATES BASED ON SPECTRAL TYPES
wrej = where(db_original.refsptn lt psptrange(0) or db_original.refsptn gt psptrange(1),cntrej)
if (cntrej gt 0) then singlerejflag(wrej) = 1
wrej = where(db_original.refsptn lt ssptrange(0) or db_original.refsptn gt ssptrange(1),cntrej)
if (cntrej gt 0) then secondaryrejflag(wrej) = 1

; REJECT YOUNG SOURCES IF DESIRED
if (keyword_set(noyoung)) then begin
 wrej = where(strpos(db_original.library,'young') gt -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrej) = 1
  secondaryrejflag(wrej) = 1
 endif
endif

; REJECT ALL BUT YOUNG SOURCES IF DESIRED
if (keyword_set(onlyyoung)) then begin
 wrej = where(strpos(db_original.library,'young') eq -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrej) = 1
  secondaryrejflag(wrej) = 1
 endif
endif

; REJECT SUBDWARFS IF DESIRED
if (keyword_set(nosubdwarfs)) then begin
 wrej = where(strpos(db_original.library,'subdwarf') gt -1 or strpos(db_original.opt_type,'sd') gt -1 or strpos(db_original.nir_type,'sd') gt -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrej) = 1
  secondaryrejflag(wrej) = 1
 endif
endif

; REJECT ALL BUT SUBDWARFS IF DESIRED
if (keyword_set(onlysubdwarfs)) then begin
 wrej = where(strpos(db_original.library,'subdwarf') eq -1 and strpos(db_original.opt_type,'sd') eq -1 and strpos(db_original.nir_type,'sd') eq -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrej) = 1
  secondaryrejflag(wrej) = 1
 endif
endif

; REJECT BLUE DWARFS IF DESIRED
if (keyword_set(noblue)) then begin
 wrejp = where(strpos(db_original.library,'blue') gt -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrejp) = 1
  secondaryrejflag(wrejp) = 1
 endif
endif

; REJECT ALL BUT BLUE DWARFS IF DESIRED
if (keyword_set(onlyblue)) then begin
 wrejp = where(strpos(db_original.library,'blue') eq -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrejp) = 1
  secondaryrejflag(wrejp) = 1
 endif
endif

; REJECT RED DWARFS IF DESIRED
if (keyword_set(nored)) then begin
 wrejp = where(strpos(db_original.library,'red') gt -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrejp) = 1
  secondaryrejflag(wrejp) = 1
 endif
endif

; REJECT ALL BUT RED DWARFS IF DESIRED
if (keyword_set(onlyred)) then begin
 wrejp = where(strpos(db_original.library,'red') eq -1,cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrejp) = 1
  secondaryrejflag(wrejp) = 1
 endif
endif

; REJECT UNPUBLISHED CRUZ SPECTRA IF DESIRED
if (keyword_set(nocruz)) then begin
 wrejp = where(db_original.data_reference eq 'CRU-NP',cntrej)
 if (cntrej gt 0) then begin
  singlerejflag(wrejp) = 1
  secondaryrejflag(wrejp) = 1
 endif
endif

primaryrejflag = singlerejflag

; -----------------------------------------------
; FIT SPECTRA TO SINGLE TEMPLATES
; -----------------------------------------------

; COMPARE DATA TO TEMPLATES
if (keyword_set(silent) eq 0) then begin
 print, ''
 print, 'Single fits:'
endif

; calculate minimum deviation between spectra with scaling
db_singles = spexbinaryfit_mask_templates(db_original,singlerejflag)
comp = db_singles.flux
if (keyword_set(silent) eq 0) then print, 'Using '+strtrim(string(db_singles.ntemplates),2)+' single templates'
nbests = nbest < db_singles.ntemplates

devsing = spexbinaryfit_fitting(stdlam, flx, comp, noise=ns, dof=dof, weight=weight, fitregions=fitregions)
if (keyword_set(difference)) then devsing(*,0) = 100.*devsing(*,0)
ssing = sort(devsing(*,0))

; print & plot best cases
  if (keyword_set(silent) eq 0) then $
  spexbinaryfit_printout, db_singles, devsing, /short, nbest=nbests, dof=dof, /single, full=full, mko=mko, nicmos=nicmos, nirc2=nirc2, wfc3=wfc3
  if (keyword_set(noplot) eq 0) then $
  spexbinaryfit_printout, db_singles, devsing, outfile=folder+prefix+'singlefits'+suffix+'.txt', nbest=nbests, dof=dof, /single, full=full, mko=mko, nicmos=nicmos, nirc2=nirc2, wfc3=wfc3, latex=latex

for i=0,nbests-1 do begin

 labels=[$
  note,$
  db_singles.name(ssing(i))+' '+db_singles.refspt(ssing(i)),$
  statstr+' = '+ strmid(strtrim(string(rndoff(devsing(ssing(i),0)*100.)/100.),2),0,4)]

 if keyword_set(difference) then nsp = ns*0. else nsp=ns
 if keyword_set(short) then labels(1) = [db_singles.shortname(ssing(i))+' '+db_singles.refspt(ssing(i))]

 if (keyword_set(noplot) eq 0) then $
  spexbinaryfit_plot, stdlam, flx, db_singles.flux(*,ssing(i))*devsing(ssing(i),1), noise=nsp, labels=labels, pcolors = pcolors, lcolors=[0,2,0,0], yrange=yrange, outfile=folder+prefix+'singlefit_'+strtrim(string(i+1),2)+suffix+'.eps', inset=inset, xrange=xrange, plotfitrange=plotfitrange, right=right, showdiff=showdiff
endfor

; create output structure and best fit variable
outputstr = create_struct($
 'wave',stdlam,$
 'flux',flx,$
 'noise',ns,$
 'note',note,$ 
 'dof',dof,$ 
 'filters',db_singles.filters,$ 
 'filter_names',db_singles.filter_names,$ 
 'n_singles',n_elements(db_singles.shortname),$
 'singles_name',db_singles.name(ssing[0:nbests-1]),$
 'singles_shname',db_singles.shortname(ssing[0:nbests-1]),$
 'singles_spt',refspt(ssing[0:nbests-1]),$
 'singles_flux',db_singles.flux(*,ssing[0:nbests-1]),$
 'singles_chi2',reform(devsing(ssing[0:nbests-1],0)),$
 'singles_scale',reform(devsing(ssing[0:nbests-1],1)))
bestfitsingle = [[stdlam],[db_singles.flux(*,ssing(0))*devsing(ssing(0),1)]] 

if (keyword_set(single)) then goto, DIAGNOSTIC
 

; -----------------------------------------------
; FIT SPECTRA TO BINARY TEMPLATES
; -----------------------------------------------

; GENERATE STRUCTURES OF COMPONENTS AND MAKE BINARIES

db_primaries = db_singles

db_secondaries = spexbinaryfit_mask_templates(db_original,secondaryrejflag)
if (keyword_set(silent) eq 0) then print, 'Using '+strtrim(string(db_secondaries.ntemplates),2)+' secondary templates'

spexbinaryfit_make_binaries, db_primaries, db_secondaries, db_binaries, jmag=jmag, hmag=hmag, kmag=kmag, liu06=liu06, bur07=bur07, loo08=loo08, fah12=fah12, dup12=dup12, faint=faint, bright=bright, deltamag=deltamag, e_deltamag=e_deltamag, deltafilter=deltafilter, nirspt=nirspt, optspt=optspt, spexspt=spexspt, sigma=sigma, calline=calline


; CALCULATE MINIMUM DEVIATION BETWEEN SPECTRA WITH SCALING
print, ''
print, 'Binary fits:'

comp = db_binaries.flux
if (keyword_set(silent) eq 0) then print, 'Using '+strtrim(string(db_binaries.ntemplates),2)+' binary templates'
nbestb = nbest < db_binaries.ntemplates

devbin = spexbinaryfit_fitting(stdlam, flx, comp, noise=ns, dof=dof, weight=weight, fitregions=fitregions)
if (keyword_set(difference)) then devbin(*,0) = 100.*devbin(*,0)
del = fltarr(n_elements(devbin(*,0)))*1.
if (n_elements(deltamag) gt 0) then del = 10.^(-0.4*(db_binaries.relmags(0,*)))
sbin = sort(devbin(*,0))
;stop


; F-TEST PROBABILITY THAT BEST BINARY IS A BETTER FIT THAN BEST SINGLE
prob = 100.*f_pdf(min(devsing(*,0))/min(devbin(*,0)),dof,dof)

; PRINT & PLOT BEST CASES
  if (keyword_set(silent) eq 0) then $
   spexbinaryfit_printout, db_binaries, devbin, /short, nbest=nbestb, dof=dof, prob=prob, full=full, mko=mko, nicmos=nicmos, nirc2=nirc2, wfc3=wfc3
  if (keyword_set(noplot) eq 0) then $
   spexbinaryfit_printout, db_binaries, devbin, outfile=folder+prefix+'binaryfits'+suffix+'.txt', nbest=nbestb, dof=dof, prob=prob, full=full, mko=mko, nicmos=nicmos, nirc2=nirc2, wfc3=wfc3, latex=latex


for i=0,nbestb-1 do begin

 labels=[$
  note,$
  db_binaries.primaries.name(db_binaries.pairs(sbin(i),0))+' '+$
   db_binaries.primaries.refspt(db_binaries.pairs(sbin(i),0)),$
  db_binaries.secondaries.name(db_binaries.pairs(sbin(i),1))+' '+$
   db_binaries.secondaries.refspt(db_binaries.pairs(sbin(i),1)),$
  statstr+' = '+ strmid(strtrim(string(rndoff(devbin(sbin(i),0)*100.)/100.),2),0,4)]

 if keyword_set(short) then begin
  labels(1:2) = [$
    db_binaries.primaries.shortname(db_binaries.pairs(sbin(i),0))+' '+$
     db_binaries.primaries.refspt(db_binaries.pairs(sbin(i),0)),$
    db_binaries.secondaries.shortname(db_binaries.pairs(sbin(i),1))+' '+$
     db_binaries.secondaries.refspt(db_binaries.pairs(sbin(i),1))]
 endif

 if keyword_set(difference) then nsp = ns*0. else nsp=ns

; NOTE: REMOVED SCALING
  if (keyword_set(noplot) eq 0) then $
   spexbinaryfit_plot, stdlam, flx, db_binaries.primaries.flux(*,db_binaries.pairs(sbin(i),0))*devbin(sbin(i),1),$
    db_binaries.secondaries.flux(*,db_binaries.pairs(sbin(i),1))*devbin(sbin(i),1)*10.^(-0.4*del),$
    db_binaries.flux(*,sbin(i))*devbin(sbin(i),1), noise=nsp, labels=labels, pcolors = pcolors, lcolors=lcolors, $
    yrange=yrange, outfile=folder+prefix+'binaryfit_'+strtrim(string(i+1),2)+suffix+'.eps', inset=inset, $
    xrange=xrange, plotfitrange=plotfitrange, right=right, showdiff=showdiff

; -----------------------------------------------
; PLOTTING EXTRAS
; -----------------------------------------------

if (keyword_set(extras)) then begin

; plot spectrum on its own
 spexbinaryfit_plot, stdlam, flx, noise=nsp, labels=[''], pcolors = pcolors, lcolors=lcolors, yrange=yrange, outfile=folder+prefix+'source.eps', inset=inset, xrange=xrange, plotfitrange=plotfitrange

; plot just primary source fit
 spexbinaryfit_plot, stdlam, flx, db_binaries.primaries.flux(*,db_binaries.pairs(sbin(i),0))*devbin(sbin(i),1), $
  noise=nsp, labels=['',labels(1)], pcolors = pcolors, lcolors=lcolors, yrange=yrange, $
  outfile=folder+prefix+'primary_'+strtrim(string(i+1),2)+suffix+'.eps', inset=inset, xrange=xrange, $
  plotfitrange=plotfitrange, showdiff=showdiff

; plot just primary and secondary source fit
 spexbinaryfit_plot, stdlam, flx, db_binaries.primaries.flux(*,db_binaries.pairs(sbin(i),0))*devbin(sbin(i),1), $
  db_binaries.secondaries.flux(*,db_binaries.pairs(sbin(i),1))*devbin(sbin(i),1)*del(sbin(i)),   $ *10.^(-0.4*del), $
  noise=nsp,labels=['',labels(1:2)], pcolors = pcolors, lcolors=lcolors, yrange=yrange, $
  outfile=folder+prefix+'secondary_'+strtrim(string(i+1),2)+suffix+'.eps', inset=inset, xrange=xrange,$
  plotfitrange=plotfitrange, showdiff=showdiff

endif
 
endfor

; -----------------------------------------------
; CREATE OUTPUT STRUCTURE FOR FURTHER ANALYSIS
; -----------------------------------------------

; update output structure
outputstr = create_struct(outputstr, $
 'n_binaries',db_binaries.ntemplates,$
 'calibration',calline,$
 'binaries_primary_name',reform(db_binaries.primaries.name(db_binaries.pairs(sbin[0:nbestb-1],0))),$
 'binaries_primary_shname',reform(db_binaries.primaries.shortname(db_binaries.pairs(sbin[0:nbestb-1],0))),$
 'binaries_secondary_name',reform(db_binaries.secondaries.name(db_binaries.pairs(sbin[0:nbestb-1],1))),$
 'binaries_secondary_shortname',reform(db_binaries.secondaries.shortname(db_binaries.pairs(sbin[0:nbestb-1],1))),$
 'binaries_primary_spt',reform(db_binaries.primaries.refspt(db_binaries.pairs(sbin[0:nbestb-1],0))),$
 'binaries_secondary_spt',reform(db_binaries.secondaries.refspt(db_binaries.pairs(sbin[0:nbestb-1],1))),$
 'binaries_primary_flux',db_binaries.primaries.flux(*,db_binaries.pairs(sbin[0:nbestb-1],0)),$
 'binaries_secondary_flux',db_binaries.secondaries.flux(*,db_binaries.pairs(sbin[0:nbestb-1],1)),$
 'binaries_flux',db_binaries.flux(*,sbin[0:nbest-1])*devbin(sbin[0:nbestb-1],1),$
 'binaries_chi2',reform(devbin(sbin[0:nbestb-1],0)),$
 'binaries_scale',reform(devbin(sbin[0:nbestb-1],1)),$
 'binaries_relmags',db_binaries.relmags(*,sbin[0:nbestb-1]))
bestfitbinary = [[stdlam],[db_binaries.flux(*,sbin(0))*devbin(sbin(0),1)]] 
 
if (keyword_set(noplot) eq 0) then save, outputstr, bestfitsingle, bestfitbinary, file=folder+prefix+'fits.dat'

; -----------------------------------------------
; DIAGNOSTIC INFORMATION
; -----------------------------------------------

DIAGNOSTIC: if (keyword_set(silent) eq 0) then begin
s = sort(db_original.shortname)
u = uniq(db_original.shortname(s))
print, ''
print, 'Number of orginal spectra used: '+strtrim(string(db_original.ntemplates),2)+' spectra of '+strtrim(string(n_elements(u)),2)+' sources'
s = sort(db_singles.shortname)
u = uniq(db_singles.shortname(s))
print, ''
print, 'Number of single templates used: '+strtrim(string(db_singles.ntemplates),2)+' spectra of '+strtrim(string(n_elements(u)),2)+' sources'
if (keyword_set(single) eq 0) then begin
 s = sort(db_binaries.secondaries.shortname)
 u = uniq(db_binaries.secondaries.shortname(s))
 print, 'Number of secondary templates used: '+strtrim(string(db_binaries.secondaries.ntemplates),2)+' spectra of '+strtrim(string(n_elements(u)),2)+' sources'
 print, 'Number of binary templates used: '+strtrim(string(db_binaries.ntemplates),2)
 print, calline
endif 
print, 'Degrees of Freedom: '+strtrim(string(dof),2)
print, 'Data saved to '+folder
print, 'Processing time: '+strtrim(string(systime(1)-time),2)+' seconds'
print, ''
endif

FINISH: return
end

