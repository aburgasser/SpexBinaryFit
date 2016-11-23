; -----------------------------------------------------------------------------------------------------
; PRO SPEXBINARYFIT_FITTING
;
; Purpose: do actual fit of spectral data to templates
; 
; history: 
;  2009 Nov 19 (AJB): initial file
;  2014 Feb 11 (AJB):
;	reordered inputs
;	changes WEIGHTS input to WEIGHT flag, made constant weights the default
; -----------------------------------------------------------------------------------------------------
function spexbinaryfit_fitting, lam, flx, cflx, noise=noise, dof=dof, weight=weight, fitregions=fitregions

on_error, 0

if (n_elements(noise) eq 0) then noise=flx*1.
if (n_elements(nbest) eq 0) then nbest=10

; regions to fit spectra data
if (n_elements(fitregions) eq 0) then fitregions=[[0.95,1.35],[1.45,1.8],[2.0,2.35]]
fitflag = lam*0
for i=0,n_elements(fitregions(0,*))-1 do begin
 w = where(lam ge fitregions(0,i) and lam le fitregions(1,i),cnt)
 if (cnt gt 0) then fitflag(w) = 1
endfor
wfit = where(fitflag eq 1,cnt)
if (cnt eq 0) then message, 'No regions available for fitting'

; weighting schemes - there is room for expansion
wts = lam*1.
if (keyword_set(weight)) then begin
   wts = abs(shift(lam,1)-lam)		; weight by size of wavelength bin
   wts(0) = median(wts)
endif
  
; DEGREES OF FREEDOM
if (n_elements(dof) eq 0) then begin
 dof = fix(total(wts(wfit))/max(wts(wfit)))-1.	; degrees of freedom  = # points -1 
 dof = dof/5. 								; to roughly compensate for slit width
endif

; run comparisons
ncomp = n_elements(cflx(0,*))
dev = fltarr(ncomp,2)
for i=long(0),ncomp-1 do begin
 tscl = total(wts(wfit)*flx(wfit)*cflx(wfit,i)/(noise(wfit)^2))/total(wts(wfit)*(cflx(wfit,i)^2)/(noise(wfit)^2))
 dev(i,*) = [total(wts(wfit)*((flx(wfit)-cflx(wfit,i)*tscl)^2)/(noise(wfit)^2))/total(wts(wfit)),tscl]
endfor

return, dev
end
