; -----------------------------------------------------------------------------------------------------
; PRO SPEXBINARYFIT_READSPECTRUM
;
; Purpose: read in spectral data file
; 
; last update: 18 Nov 2009
; -----------------------------------------------------------------------------------------------------

function spexbinaryfit_readspectrum, file   ; , normrng=normrng

on_error, 0

;if (n_elements(normrng) eq 0) then normrng = [[0.9,1.4],[1.5,1.7],[2.0,2.3]]
normrng = [[0.9,1.4],[1.5,1.7],[2.0,2.3]]

; default S/N
snscl = 0.05

if (file_search(file) eq '') then message, file+' not found'

; DETERMINE FILE FORMAT
tmp = str_sep(file,'.')
file_format = tmp(n_elements(tmp)-1)

; READ FITS OR TEXT FILES
if (strpos(file_format,'fit') ge 0) then $
 fits_read, file, data, header else $
 data = readarr(file,/comment,/flt,/space)
; transpose if necessary
if n_elements(data(0,*)) gt 5 then data = transpose(data)
lam = data(*,0)
flx = data(*,1)
if (n_elements(data(0,*)) gt 2) then noise = data(*,2) else noise = abs(flx*snscl)

; PROCESS ARRAYS
wnan = where(strlowcase(strtrim(string(lam),2)) ne 'nan' and strlowcase(strtrim(string(flx),2)) ne 'nan',cntnan)
lam = lam (wnan)
flx = flx(wnan)
noise = noise(wnan)

; NORMALIZE SPECTRUM
w = where($
 lam ge normrng(0,0) and lam le normrng(1,0) or $
 lam ge normrng(0,1) and lam le normrng(1,1) or $
 lam ge normrng(0,2) and lam le normrng(1,2))
scl = max(smooth(flx(w),11))
flx = flx/scl
noise = noise/scl

; some spectra are in SNR instead of noise
if (median(flx/noise) lt 1.) then noise=flx/noise

; get rid of infinities and zero noise spots
w = where(noise le 0. or finite(noise) eq 0,cnt)
if (cnt gt 0) then noise(w) = median(noise)

FINISH: return, [[lam],[flx],[noise]]
end


