;-----------------------------------------------------------------------------------------------------
 PRO SPEXBINARYFIT_GENERATE_LIBRARY, $
;
; Purpose: produce a library of single templates
; 
; HERSTORY
; created: 2014 Oct 27 by Adam J. Burgasser (based on spexbinaryfit_singles_library
;
; TO DO
; - verify that filters are correct
; - add output of template information to latex table
; -----------------------------------------------------------------------------------------------------

savefile, pfold=pfold, dbfile=dbile, datafolder=datafolder

on_error, 0

print, ''
print, 'Building single spectral library'
tb='	'
f=''

; PLOTTING PARAMETERS
!p.font=0
!p.thick=4
!x.thick=3
!y.thick=3

; CHECK INPUT PARAMETERS
if (n_elements(pfold) eq 0) then pfold = '/Users/adam/idl/spexbinaryfit/'
if (n_elements(dbfile) eq 0) then dbfile = '/Users/adam/web/browndwarfs/spexprism/db/db_spexprism.txt'
if (n_elements(datafolder) eq 0) then datafolder = '/Users/adam/spectra/spex/prism/'
savefile = pfold+'templates/template_singles.dat'
plotfile = pfold+'templates/templates_singles.ps'
latexfle = pfold+'templates/templates_unpublished.tex'

; FILTERS USED FOR SOURCE MAGNITUDES
; 2MASS JHKs, MKO JHKKsKp, WFCAM Y, 
; NIRC2 JHHcKKcKsKp,
; WIRC CH4s, CH4l,
; WFPC2 F1042xQE,
; NICMOS NIC1 F090M, F108N,F110M,F110W,F113N,F140W,F145M,F160W,F165M,F170M,
; WFC3 F127M, F139M, F164N, F167N
; VLT HKs

filters = [7,8,9,41,42,43,44,45,77,82,83,84,85,86,87,88,98,99,36,55,58,59,60,61,62,63,64,66,68,78,79,80,81,100,101]

filter_names = ['2MASS J','2MASS H','2MASS Ks','MKO J','MKO H','MKO K','MKO Ks','MKO Kp','WFCAM Y','NIRC2 J','NIRC2 H','NIRC2 Hc','NIRC2 K','NIRC2 Kc','NIRC2 Ks','NIRC2 Kp','WIRC CH4s','WIRC CH4l','WFPC2 F1042','NIC1 F090M','NIC1 F108N','NIC1 F110M','NIC1 F110W','NIC1 F113N','NIC1 F140W','NIC1 F145M','NIC1 F160W','NIC1 F165M','NIC1 F170M','WFC3 F127M','WFC3 F139M','WFC3 F164N','WFC3 F167N','VLT H','VLT Ks']

; resort filters
s = sort(filters)
filters = filters(s)
filter_names = filter_names(s)

; read in spex prism database file
READDB: db = spexbinaryfit_readstruct(dbfile)

; set numeric spectral types and shortnames - nir sptn for no optical classification or >=L9
db = create_struct(db,$
 'shortname',spexbinaryfit_designation2shortname(db.designation),$
 'osptn',spexbinaryfit_numericspt(db.opt_type), $
 'nsptn',spexbinaryfit_numericspt(db.nir_type), $
 'spts', db.opt_type,  $
 'sptn',spexbinaryfit_numericspt(db.opt_type))
w = where(db.osptn le 0 or db.nsptn ge 29.)
db.spts(w) = db.nir_type(w)
db.sptn(w) = db.nsptn(w)

; INITIALIZE NEW SET OF TEMPLATES

; set up template flag if needed
 if (max(strpos(tag_names(db),'TEMPLATE')) eq -1) then $
  db = create_struct(db,'template',strarr(n_elements(db.data_file))+'Y')

; set templates: SpT >= M4, no poor data, no binaries, no giants, no subdwarfs
 wtemplate = where($
;  strupcase(db.template) eq 'Y' and $		; already selected as possible template - NOW OFF
  db.quality_flag eq 'OK' and $				; good spectra
  db.data_reference ne 'KIR-NP' and $		; no unpublished Kirkpatrick spectra
  db.data_reference ne 'CRU-NP' and $		; no unpublished Cruz spectra
  strpos(db.discovery_reference,'-NP') eq -1 and $	; no "unreported" sources 
  strpos(db.opt_type,'galaxy') eq -1 and $		; not a galaxy
  strpos(db.nir_type,'galaxy') eq -1 and $		
  strpos(db.opt_type,'WD') eq -1 and $		; not a white dwarf
  strpos(db.nir_type,'WD') eq -1 and $		
  strpos(db.opt_type,'B') eq -1 and $			; not a B star
  strpos(db.opt_type,'C') eq -1 and $			; not a carbon star
  strpos(db.opt_type,'K') eq -1 and $			; not a K star
  strpos(db.library,'giant') eq -1 and $		; not a giant
  strpos(db.library,'binary') eq -1 and $		; not a binary or spectral binary (candidate)
  strlen(db.designation) gt 10 and $			; full designation
  db.sptn ge 15. and db.sptn le 39.,ntemplates)				; M5 <= SpT <= T9

; READ IN DATA AND NORMALIZE
; set up a structure for just the template sources
db_singles = create_struct('ntemplates',ntemplates)
names = tag_names(db)
for i=0, n_elements(names)-1 do db_singles = create_struct(db_singles,names[i],db.(i)[wtemplate])

print, ' Number of single templates: '+strtrim(string(ntemplates),2)
stdmags = fltarr(n_elements(filters),ntemplates)
spex_sptn = fltarr(ntemplates)
spex_sptn_e = fltarr(ntemplates)
spex_type = strarr(ntemplates)


; initiate plot file
if (n_elements(plotfile) ne 0) then begin
 plotps, plotfile, /open, /color
 !p.multi=[0,2,2]
endif

; print out a latexfile of unpublished templates
if (n_elements(latexfile) ne 0) then begin
 s = sort(db.designation(wtemplate)+db.name(wtemplate)+db.observation_date(wtemplate))
 spexbinaryfit_latextable, db, latexfile, index=wtemplate(s), title='SpeX Spectral Templates'
endif

; READ IN SPECTRA, FLUX CALIBRATE, ASSIGN MAGNITUDES & PLOT
for i=0,ntemplates-1 do begin

; read in
 spec = spexbinaryfit_readspectrum(datafolder+db.data_file(wtemplate(i)))

; initialize
 if (i eq 0) then begin
  wave = spec(*,0)
  flux = fltarr(n_elements(wave),ntemplates)
  noise = fltarr(n_elements(wave),ntemplates)
  wnorm = where((wave gt 0.9 and wave lt 1.35) or (wave gt 1.5 and wave lt 1.8) or (wave gt 2.0 and wave lt 2.4))
 endif
 
; interpolate flux and noise onto wave scale, determine flux scaling
 flux(*,i) = interpol(spec(*,1),spec(*,0),wave)
 scl = max(smooth(flux(wnorm,i),11))
 flux(*,i) = flux(*,i)/scl
 if (n_elements(spec(0,*) gt 2)) then noise(*,i) = interpol(spec(*,2),spec(*,0),wave)/scl $
   else noise(*,i) = flux(*,i)*0.+0.05		; guess for no noise measure

; DETERMINE MAGNITUDES   
 stdmags(*,i) = spexbinaryfit_filtflux(wave,flux(*,i),filters,1,/force,/silent)

; DETERMINE SPEX SPECTRAL TYPE   
 spexbinaryfit_nirclassify, wave, flux(*,i), sp, spe, sps, sptinit = db.sptn(wtemplate(i))
 spex_type(i) = sps
 spex_sptn(i) = sp
 spex_sptn_e(i) = spe

; PLOT SPECTRUM  
  plot, wave, flux(*,i), xrange=[0.8,2.5], yrange=[-0.02,1.2], /xsty, /ysty, xtitle='!3Wavelength (!7m!3m)', ytitle='Normalized !3F!D!9l!3!N', charsize=1., title=db.data_file(wtemplate(i))
  oplot, wave, noise(*,i), color=100
  xyouts, 2.45, 1.1, db.shortname(wtemplate(i)), align=1
  xyouts, 2.45, 1.0, 'SpT = '+db.opt_type(wtemplate(i))+'/'+db.nir_type(wtemplate(i)), align=1
  xyouts, 0.9, 1.1, db.quality_flag(wtemplate(i)), charsize=1.5, align=0 
 
endfor

; CREATE STRUCTURE AND SAVE
db_singles = create_struct(db_singles,$
 'wavelength',wave,$
  'flux',flux,$ 
  'noise',noise,$
  'mags',stdmags,$
  'filters',filters,$
  'filter_names',filter_names,$
  'spex_sptn',spex_sptn,$
  'spex_sptn_e',spex_sptn_e)
db_singles.spex_type = spex_type

save, db_singles, file=savefile

if (n_elements(plotfile) ne 0) then begin
 plotps, plotfile, /close
 !p.multi=0
endif

FINISH: return
end

