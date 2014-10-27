;-----------------------------------------------------------------------------------------------------
  PRO SPEXBINARYFIT_MAKE_BINARIES, $
;
; Purpose: produce a flux calibrated spectral library of binary templates
; 
; CREATED: 2014 Oct 27 (based on spexbinaryfit_library_binaries)
; -----------------------------------------------------------------------------------------------------

db_primaries, db_secondaries, db_systems, calline=calline, savefile=savefile, nirspt=nirspt, optspt=optspt, spexspt=spexspt, bright=bright, faint=faint, liu06=liu06, bur07=bur07, loo08=loo08, fah12=fah12, dup12=dup12, jmag=jmag, hmag=hmag, kmag=kmag, deltamag=deltamag, e_deltamag = e_deltamag, deltafilter=deltafilter, sigma=sigma 

on_error, 0

tb='	'
f=''
resolved = 0

if (keyword_set(jmag) eq 0 and keyword_set(hmag) eq 0 and keyword_set(kmag) eq 0) then jmag = 1

pairs = intarr(db_primaries.ntemplates*db_secondaries.ntemplates,2)
hybflx = fltarr(n_elements(db_primaries.wavelength),db_primaries.ntemplates*db_secondaries.ntemplates)
relmags = fltarr(n_elements(db_primaries.filters(*,0)),db_primaries.ntemplates*db_secondaries.ntemplates)*0.-99.

; choose which spectral type to base combination on 
case 1 of
 keyword_set(optspt): begin
  refsptp = db_primaries.opt_type
  refsptpn = db_primaries.osptn
  refspts = db_secondaries.opt_type
  refsptsn = db_secondaries.osptn
 end
 keyword_set(nirspt): begin
  refsptp = db_primaries.nir_type
  refsptpn = db_primaries.nsptn
  refspts = db_secondaries.nir_type
  refsptsn = db_secondaries.nsptn
 end
 keyword_set(spexspt): begin
  refsptp = db_primaries.spex_type
  refsptpn = db_primaries.spex_sptn
  refspts = db_secondaries.spex_type
  refsptsn = db_secondaries.spex_sptn
 end
 else: begin
  refsptp = db_primaries.spts
  refsptpn = db_primaries.sptn
  refspts = db_secondaries.spts
  refsptsn = db_secondaries.sptn
 end
endcase


; COMBINE BINARIES BY RELATIVE MAGNITUDES

if (n_elements(deltamag) ne 0 and n_elements(deltafilter) ne 0) then begin
 resolved = 1
 filters_used = db_primaries.filters*0

; assume uncertainty of 0.1 mag if not provided
 if (n_elements(e_deltamag) eq 0) then e_deltamag = deltamag*0.+0.1

; assume sigma = 3 if not provided
 if (n_elements(sigma) eq 0) then sigma=3

; determine correct filters - only do for primary structure (assume the same)
 for i=0,n_elements(deltafilter)-1 do begin
  wf = where(db_primaries.filters eq deltafilter(i),cnt)
  if (cnt gt 0) then filters_used(wf(0)) = 1 else begin
   print, 'WARNING! Missing filter in measurements not used - re-run templates'
   pop, deltamag, i, /index
   pop, e_deltamag, i, /index
   pop, deltafilter, i, /index
 endelse
 endfor
 wfok = where(filters_used eq 1)

; create relative scalings 
; secondary SpT >= primary SpT - 2, identical source pairs OK
 kkk=0
 photflag = intarr(db_primaries.ntemplates*db_secondaries.ntemplates)*0
 for i=0,db_primaries.ntemplates-1 do begin
  w = where(refsptsn ge refsptpn(i)-2,cnt)
  if (cnt gt 0 and refsptp(i) ne '') then begin
   for j=0,cnt-1 do begin
    relmags(*,kkk+j) = db_secondaries.mags(*,w(j))-db_primaries.mags(*,i)
    offset = total((relmags(wfok,kkk+j)-deltamag)/(e_deltamag^2))/total(1./(e_deltamag)^2)
    relmags(*,kkk+j) = relmags(*,kkk+j) - offset
    hybflx(*,kkk+j) = db_primaries.flux(*,i)+db_secondaries.flux(*,w(j))*10.^(0.4*offset)	; primary + scaled secondary
    pairs(kkk+j,*) = [i,w(j)]
    if (max(abs(relmags(wfok,kkk+j)-deltamag)/e_deltamag) le sigma) then photflag(kkk+j) = 1
   endfor
   kkk=kkk+cnt
  endif
 endfor 
 filters_used_p = intarr(db_primaries.ntemplates)*0+filters_used(wfok(0))
 filters_used_s = intarr(db_secondaries.ntemplates)*0+filters_used(wfok(0))
goto, FINISH

endif else begin


; COMBINE BINARIES BY ABSOLUTE MAGNITUDE RELATION

absmagsp = fltarr(db_primaries.ntemplates)
absmagss = fltarr(db_secondaries.ntemplates)
filters_used_p = fltarr(db_primaries.ntemplates)*0.-99.
filters_used_s = fltarr(db_secondaries.ntemplates)*0.-99.

; J (2MASS) for M5 - L5  from Cruz et al. 2003 - a default
coeff = [-4.410,5.043,-0.6193,0.03453,-0.0006892]
sptrng = [15,25.]
fid = 7

w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
wf = where(db_primaries.filters eq fid)
filters_used_p(w) = wf(0)
absmagsp(w) = poly(refsptpn(w)-10.,coeff)
w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
wf = where(db_secondaries.filters eq fid)
filters_used_s(w) = wf(0)
absmagss(w) = poly(refsptsn(w)-10.,coeff)

case 1 of

; J (MKO) Dupuy & Liu 2012 - for M6-T9
 (keyword_set(dup12) eq 1 and keyword_set(jmag) eq 1): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-M6, '
  calline = calline+'Using Dupuy & Liu (2012) M_J/SpT  relation for M6-T9'
  coeff = [-2.83129e1,1.63986e1,-2.74405e0,2.32771e-1,-1.03332e-2,2.27641e-4,-1.94920e-6]
  sptrng = [16.,39.]
  fid = 41
  sptoff = 10.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
 
; H (MKO) Dupuy & Liu 2012 - for M6-T9
 (keyword_set(dup12) eq 1 and keyword_set(hmag) eq 1): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-M6, '
  calline = calline+'Using Dupuy & Liu (2012) M_H/SpT  relation for M6-T9'
  coeff = [-2.97306e-1,1.69138e1,-2.85705e0,2.45209e-1,-1.10960e-2,2.51601e-4,-2.24083e-6]
  sptrng = [16.,39.]
  fid = 42
  sptoff = 10.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
 
; K (MKO) Dupuy & Liu 2012 - for M6-T9
 (keyword_set(dup12) eq 1 and keyword_set(kmag) eq 1): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-M6, '
  calline = calline+'Using Dupuy & Liu (2012) M_K/SpT  relation for M6-T9'
  coeff = [-1.52200e1,1.01248e1,-1.63930e0,1.351771-1,-5.84342e-3,1.25731e-4,-1.04935e-6]
  sptrng = [16.,39.]
  fid = 41
  sptoff = 10.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
 
; J (MKO) Faherty et al. 2012 - for >L2
 (keyword_set(fah12) eq 1 and keyword_set(jmag) eq 1): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L2, '
  calline = calline+'Using Faherty et al. (2012) M_J/SpT faint relation for L2-T8'
  coeff = [14.7427,-1.95991,0.272793,-0.0128331,0.000202557]
  sptrng = [22.,38.]
  fid = 41
  sptoff = 10.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
 
; H (MKO) Faherty et al. 2012 - for >L2
 (keyword_set(fah12) eq 1 and keyword_set(jmag) eq 1): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L2, '
  calline = calline+'Using Faherty et al. (2012) M_H/SpT faint relation for L2-T8'
  coeff = [13.0639,-1.55464,0.223912,-0.0107219,0.000174458]
  sptrng = [22.,38.]
  fid = 42
  sptoff = 10.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
 
; K (MKO) Faherty et al. 2012 - for >L2
 (keyword_set(fah12) eq 1 and keyword_set(jmag) eq 0 and keyword_set(hmag) eq 0): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L2, '
  calline = calline+'Using Faherty et al. (2012) M_K/SpT faint relation for L2-T8'
  coeff = [9.57924,-0.436982,0.0903377,-0.00457965,0.0000805685]
  sptrng = [22.,38.]
  fid = 43
  sptoff = 10.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
 
 ; J (2MASS) Looper et al. 2008 - for >L2
; from Looper et al. 2008
 (keyword_set(loo08) eq 1 and keyword_set(jmag) eq 1): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L2, '
  calline = calline+'Using Looper et al. (2008) M_J/SpT faint relation for L2-T8'
  coeff = [11.817,1.255e-1,3.690e-2,1.663e-2,-3.915e-3,2.595e-4,-5.462e-6]
  sptrng = [22.,38.]
  fid = 7
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
 
; Ks (2MASS) Looper et al. 2008 - for >L2
; from Looper et al. 2008
 (keyword_set(loo08) eq 1 and keyword_set(jmag) eq 0): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L2, '
  calline = calline+'Using Looper et al. (2008) M_Ks/SpT faint relation for L2-T8'
  coeff = [10.521,7.369e-2,2.565e-2,1.299e-2,-2.864e-3,1.911e-4,-4.104e-6]
  sptrng = [22.,38.]
  fid = 9
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
  
; J (MKO) Liu et al. 2006 "faint"- for >L5
; from Liu et al. 2006, excluding known + possible binaries
 (keyword_set(bright) eq 0 and keyword_set(jmag) eq 1): begin
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L5, '
  calline = calline+'Using Liu et al. (2006) M_J/SpT faint relation for L5-T8'
  coeff = [11.359,3.174e-1,1.102e-1,-1.877e-2,9.169e-4,-1.233e-5]
  sptrng = [25.,38.]
  fid = 41
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end

; J (MKO) Liu et al. 2006 "bright"-  for >L5
; from Liu et al. 2006, excluding known binaries
 (keyword_set(bright) eq 1and keyword_set(jmag) eq 1): begin
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L5, '
  calline = calline+'Using Liu et al. (2006) M_J/SpT bright relation for L5-T8'
  coeff = [11.746,-2.259e-1,3.229e-1,-5.155e-2,2.966e-3,-5.648e-5]
  sptrng = [25.,38.]
  fid = 41
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end

; K (MKO) Liu et al. 2006 "faint"-  for >L5
; from Liu et al. 2006, excluding known and possible binaries
 (keyword_set(bright) eq 0 and keyword_set(liu06) eq 1 and keyword_set(jmag) eq 0): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L5, '
  calline = calline+'Using Liu et al. (2006) M_K/SpT faint relation for L5-T8'
  coeff = [10.182,3.688e-1,2.028e-2,-4.488e-3,2.301e-4,-2.466e-6]
  sptrng = [25.,38.]
  fid = 43
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end

; K (MKO) Liu et al. 2006 "bright"-  for >L5
; from Liu et al. 2006, excluding known binaries
 (keyword_set(bright) eq 1 and keyword_set(liu06) eq 1 and keyword_set(jmag) eq 0): begin 
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L5, '
  calline = calline+'Using Liu et al. (2006) M_K/SpT bright relation for L5-T8'
  coeff = [10.731,-3.964e-1,3.150e-1,-4.912e-2,2.994e-3,-6.179e-5]
  sptrng = [25.,38.]
  fid = 43
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end

; K (MKO) Burgasser 2007 - "bright" -  for >L5
; based on polynomial fit to MKO K data (e.g., Knapp et al. 2004 and references
; therein/thereout) and parallax data from Dahn et al. 2002, Tinney et al. 2003
; and Vrba et al. 2004; uncertainty < 0.2 mag (36 sources)
; includes component magnitudes for Eps Indi Bab (McCaughrean et al. 2004)
; SD0423-0414 & SD1021-0304 (Burgasser et al. 2006), 
; and Kelu 1 AB (Liu & Leggett 2005)
 (keyword_set(bright) eq 1 and keyword_set(bur07) eq 1 and keyword_set(jmag) eq 0): begin  
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L5, '
  calline = calline+'Using Burgasser (2007) bright M_K/SpT relation for L5-T8'
  coeff = [10.2389, 0.912754, -0.553959, 0.185108, -0.0283012, 0.00214514, -7.88129e-05, 1.12414e-06]
  sptrng = [25.,38.]
  fid = 43
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end
  
; K (MKO) Burgasser 2007 "faint" -  for >L5 - DEFAULT
; based on polynomial fit to MKO K data (e.g., Knapp et al. 2004 and references
; therein/thereout) and parallax data from Dahn et al. 2002, Tinney et al. 2003
; and Vrba et al. 2004; uncertainty < 0.2 mag (34 sources)
; includes component magnitudes for Eps Indi Bab (McCaughrean et al. 2004)
; SD0423-0414 & SD1021-0304 (Burgasser et al. 2006), 
; and Kelu 1 AB (Liu & Leggett 2005)
; excludes SD1254-0122, 2M0559-1404, and 2M0937+29
 else: begin
  calline = 'Using Cruz et al. (2003) M_J/SpT relation for M5-L5, '
  calline = calline+'Using Burgasser (2007) faint M_K/SpT relation for L5-T8'
  coeff = [10.4458, 0.232154, 0.0512942, -0.0402365, 0.0141398, -0.00227108,  0.000180674, -6.98501e-06, 1.05119e-07]
  sptrng = [25.,38.]
  fid = 43
  sptoff = 20.

  w = where(refsptpn ge sptrng(0) and refsptpn le sptrng(1))
  wf = where(db_primaries.filters eq fid)
  filters_used_p(w) = wf(0)
  absmagsp(w) = poly(refsptpn(w)-sptoff,coeff)
  w = where(refsptsn ge sptrng(0) and refsptsn le sptrng(1))
  wf = where(db_secondaries.filters eq fid)
  filters_used_s(w) = wf(0)
  absmagss(w) = poly(refsptsn(w)-sptoff,coeff)
 end

endcase
filters_used = [fid]

; SCALE THE SPECTRA
w = where(filters_used_p ne -99.,cnt)
for i=0,cnt-1 do begin
 db_primaries.flux(*,w(i)) = db_primaries.flux(*,w(i))*10.^(0.4*(db_primaries.mags(filters_used_p(w(i)),w(i))-absmagsp(w(i))))
 db_primaries.mags(*,w(i)) = db_primaries.mags(*,w(i))-db_primaries.mags(filters_used_p(w(i)),w(i))+absmagsp(w(i))
endfor

w = where(filters_used_s ne -99.,cnt)
for i=0,cnt-1 do begin
 db_secondaries.flux(*,w(i)) = db_secondaries.flux(*,w(i))*10.^(0.4*(db_secondaries.mags(filters_used_s(w(i)),w(i))-absmagss(w(i))))
 db_secondaries.mags(*,w(i)) = db_secondaries.mags(*,w(i))-db_secondaries.mags(filters_used_s(w(i)),w(i))+absmagss(w(i))
endfor

; CREATE BINARY TEMPLATES 
; secondary SpT >= primary SpT - 2, identical source pairs rejected
kkk=0
for i=0,db_primaries.ntemplates-1 do begin
 w = where(refsptsn ge refsptpn(i)-2 and db_secondaries.shortname ne db_primaries.shortname(i) and filters_used_s ne -99.,cnt)
 if (cnt gt 0 and refsptp(i) ne '' and filters_used_p(i) ne -99.) then begin
  for j=0,cnt-1 do begin
   relmags(*,kkk+j) = db_secondaries.mags(*,w(j))-db_primaries.mags(*,i)
   hybflx(*,kkk+j) = db_primaries.flux(*,i)+db_secondaries.flux(*,w(j))
   pairs(kkk+j,*) = [i,w(j)]
  endfor
  kkk=kkk+cnt
 endif
endfor 

endelse


FINISH: hybflx = hybflx(*,0:kkk-1)
pairs = pairs(0:kkk-1,*)
relmags = relmags(*,0:kkk-1)

; report statistics
print, ' Total number of binary combinations created: '+strtrim(string(kkk),2)
ntemplates = kkk

if (resolved) then begin
 w = where(photflag eq 1,cnt)
 if (cnt gt 0) then begin
  hybflx = hybflx(*,w)
  pairs = pairs(w,*)
  relmags = relmags(*,w)
  print, ' Total number of binary combinations permitted: '+strtrim(string(cnt),2)
  ntemplates = cnt
 endif else message, 'No binaries match photometric constraints'
endif

; CREATE STUCTURE 
db_systems = create_struct('ntemplates',ntemplates,'primaries', db_primaries, 'secondaries', db_secondaries, 'wavelength', db_primaries.wavelength, 'filter_names', db_primaries.filter_names, 'flux', hybflx, 'pairs', pairs, 'filters_used', filters_used, 'relmags', relmags)

; SAVE IF DESIRED
if (n_elements(savefile) ne 0) then $
 save, db_systems, file=savefile

return
end

