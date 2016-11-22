;------------------------------------------------------------------------
function SPEXBINARYFIT_FILTFLUX, $
;
; purpose: computes magnite flux from a spectrum using library of filter
;   passbands and spectrum of vega
;
; usage:
;   IDL> mag = filtflux(lambda, flux, filter, units, /info, /energy, 
;	  /silent, /photons, /ab, /force, /mean, 
;	  PROFILE=[filter transmission profile])
;
; input:
;   lambda = wavelength, prefered microns
;   flux = flux in units listed below, erg/cm2/s/um prefered
;   filter = filter number or name, from list below, can be array
;   units = units number or name, from list below
;
; keywords:
;   INFO - provide information on use of program
;	SILENT - keep silent
;	ENERGY - report back luminosity in given filter band in units of erg/cm2/s
;	PHOTONS - report number of photons/s in given filter band
;	AB - for AB magnitudes
;	FORCE - forces filter profile to fit within spectrum
;	MEAN - report mean value in given units
;	PROFILE - input filter transmission profile, must be array of form 
;	 [[wavelength],[transmission]
;
; output:
;   mag = logarithmic mag scaled to vega
;
; external functions:
;   readcol - from the GSFC astronomy library, an associated programs
;
; authored: 
;   1/15/01  Adam Burgasser @ Caltech
;   5/7/01   implemented int_tabulated and fixed errors
;   9/17/01  added AB mag calculation 
;   10/21/02 added FORCE keyword
;   2/8/06   added MEAN keyword
;   8/25/06  fixed flux calibration to include lambda term as per
;    Cushing et al. (2005) - not necessary for 2MASS RSR
;   1/2/08   added PHOTONS and PROFILE keywords
;   6/17/08  repackaged for use by UROP students, used Kurucz Vega spectrum,
;	 changed filter and units to accept strings
;   1/20/09  add WFCAM Z and Y filters
;   10/11/14 added WISE filters
;   10/27/14 ported to spexbinaryfit_filtflux
;------------------------------------------------------------------------

 lambda, flux, filter, units, info=info, $
 silent=silent, photons=photons, energy=energy, ab=ab, force=force, $
 mean=mean, profile=profile, filtfolder=filtfolder

; **********************************

; SET THE FOLLOWING LINE TO THE FOLDER IN WHICH THE FILTERS FILES ARE CONTAINED
if (keyword_set(filtfolder) eq 0) then filtfolder = '/Users/adam/idl/spexbinaryfit/filters/'

; **********************************

; check for filter folder
if (file_search(filtfolder) eq '') then begin
 print, 'FILTFLUX: Cannot find filter folder '+filtfolder+'.  Please sure this is correctly specified in filtflux.pro; EXITING'
 return, -99.
endif

; basic setup
resp = ''
unitsnames = ['','erg/cm2/s/um','erg/cm2/s/A','W/m2/um','W/m2/A','erg/cm2/s/Hz','W/m2/Hz','Jy']
rejflag = 0
vfile = 'vega_kurucz.txt'

; ---------------------------
; FILTER DEFINITIONS
; ---------------------------
filtnames = ['',$
 '2MASS J filt',$
 '2MASS H filt',$
 '2MASS Ks filt',$
 '2MASS J optical resp',$
 '2MASS H optical resp',$
 '2MASS Ks optical resp',$
 '2MASS J optical+atm resp (RSR)',$
 '2MASS H optical+atm resp (RSR)',$
 '2MASS Ks optical+atm resp (RSR)',$
 'LCO J',$
 'LCO H',$
 'LCO K',$
 'LCO Ks',$
 'Bessel B',$
 'Bessel V',$
 'Bessel R',$
 'Bessel I',$
 'Bessel J',$
 'Bessel H',$
 'Bessel K',$
 'Bessel L',$
 'Bessel Lp',$
 'Bessel M',$
 'N',$
; 'Sloan i',$
; 'Sloan z',$
 'Gunn r',$
 'Gunn i',$
 'Gunn z',$
 'Gunn r X CCD13 QE',$
 'Gunn i X CCD13 QE',$
 'Gunn z X CCD13 QE',$
 'Cousins i',$
 'Cousins i X CCD13 QE',$
 'WFPC F814',$
 'WFPC F1042',$
 'WFPC F814 X QE',$
 'WFPC F1042 X QE',$
 'Y P60',$
 'Y NIRC',$
 'UFTI I + ATM',$
 'UFTI Z + ATM',$
 'MKO J + ATM',$
 'MKO H + ATM',$
 'MKO K + ATM',$
 'MKO Ks',$
 'MKO Kp',$
 'MKO Lp + ATM',$
 'MKO Mp + ATM',$
 'MKO J',$
 'MKO H',$
 'MKO K',$
 'IRAC 1',$
 'IRAC 2',$
 'IRAC 3',$
 'IRAC 4',$
 'NICMOS F090m (NIC1)',$
 'NICMOS F095n (NIC1)',$
 'NICMOS F097n (NIC1)',$
 'NICMOS F108n (NIC1)',$
 'NICMOS F110m (NIC1)',$
 'NICMOS F110w (NIC1)',$
 'NICMOS F113n (NIC1)',$
 'NICMOS F140w (NIC1)',$
 'NICMOS F145m (NIC1)',$
 'NICMOS F160w (NIC1)',$
 'NICMOS F164n (NIC1)',$
 'NICMOS F165m (NIC1)',$
 'NICMOS F166n (NIC1)',$
 'NICMOS F170m (NIC1)',$
 'NICMOS F187n (NIC1)',$
 'NICMOS F190n (NIC1)',$
 'SDSS u',$
 'SDSS g',$
 'SDSS r',$
 'SDSS i',$
 'SDSS z',$
 'WFCAM Z',$
 'WFCAM Y',$
 'WFC3 F127M',$
 'WFC3 F139M',$
 'WFC3 F164N',$
 'WFC3 F167N',$
 'NIRC2 J',$
 'NIRC2 H',$
 'NIRC2 Hcont',$
 'NIRC2 K',$
 'NIRC2 Kcont',$
 'NIRC2 Ks',$
 'NIRC2 Kp',$
 'NIRC2 Lp',$
 'NIRC2 Ms',$
 'WIRC Jcont',$
 'WIRC Hcont',$
 'WIRC Kcont',$
 'WIRC Pa beta',$
 'WIRC Br gamma',$
 'WIRC Fe II',$
 'WIRC CO',$
 'WIRC CH4s',$
 'WIRC CH4l',$
 'VLT CONICA H',$
 'VLT CONICA Ks',$
 'WISE W1 (RSR)',$
 'WISE W2 (RSR)',$
 'WISE W3 (RSR)',$
 'WISE W4 (RSR)']

filtfiles = ['',$
 'j_2mass_filter.txt',$
 'h_2mass_filter.txt',$
 'ks_2mass_filter.txt',$
 'j_2mass_optresp.txt',$
 'h_2mass_optresp.txt',$
 'ks_2mass_optresp.txt',$
 'j_2mass_rsr.txt',$
 'h_2mass_rsr.txt',$
 'ks_2mass_rsr.txt',$
 'j_lco.txt',$
 'h_lco.txt',$
 'k_lco.txt',$
 'ks_lco.txt',$
 'B_Bessel.txt',$
 'V_Bessel.txt',$
 'r_bessel.txt',$
 'i_bessel.txt',$
 'j_bessel.txt',$
 'h_bessel.txt',$
 'k_bessel.txt',$
 'l_bessel.txt',$
 'lp_bessel.txt',$
 'm_bessel.txt',$
 'newN_band.txt',$
; 'i_sloan.txt',$
; 'z_sloan.txt',$
 'gunnr.txt',$
 'gunni.txt',$
 'gunnz.txt',$
 'ccd13gunnr.txt',$
 'ccd13gunni.txt',$
 'ccd13gunnz.txt',$
 'i_cousins.txt',$
 'i_cousins_ccd13.txt',$
 'f814w.txt',$
 'f1042m.txt',$
 'wfpc2_f814w.txt',$
 'wfpc2_f1042m.txt',$
 'y_p60.txt',$
 'y_nirc.txt',$
 'ufti_i_dich_atm_det.txt',$
 'ufti_z_dich_atm.txt',$
 'mko_j_atm.txt',$
 'mko_h_atm.txt',$
 'mko_k_atm.txt',$
 'mko_ks.txt',$
 'mko_kp.txt',$
 'mko_lp_atm.txt',$
 'mko_mp_atm.txt',$
 'mko_j.txt',$
 'mko_h.txt',$
 'mko_k.txt',$
 'irac1.dat',$
 'irac2.dat',$
 'irac3.dat',$
 'irac4.dat',$
 'nic1_f090m.txt',$
 'nic1_f095n.txt',$
 'nic1_f097n.txt',$
 'nic1_f108n.txt',$
 'nic1_f110m.txt',$
 'nic1_f110w.txt',$
 'nic1_f113n.txt',$
 'nic1_f140w.txt',$
 'nic1_f145m.txt',$
 'nic1_f160w.txt',$
 'nic1_f164n.txt',$
 'nic1_f165m.txt',$
 'nic1_f166n.txt',$
 'nic1_f170m.txt',$
 'nic1_f187n.txt',$
 'nic1_f190n.txt',$
 'u_sdss.txt',$
 'g_sdss.txt',$
 'r_sdss.txt',$
 'i_sdss.txt',$
 'z_sdss.txt',$
 'z_wfcam.txt',$
 'y_wfcam.txt',$
 'wfc3_F127M.txt',$
 'wfc3_F139M.txt',$
 'wfc3_F164N.txt',$
 'wfc3_F167N.txt',$
 'nirc2-j.txt',$
 'nirc2-h.txt',$
 'nirc2-hcont.txt',$
 'nirc2-k.txt',$
 'nirc2-kcont.txt',$
 'nirc2-ks.txt',$
 'nirc2-kp.txt',$
 'nirc2-lp.txt',$
 'nirc2-ms.txt',$
 'wirc_jcont.txt',$
 'wirc_hcont.txt',$
 'wirc_kcont.txt',$
 'wirc_pabeta.txt',$
 'wirc_brgamma.txt',$
 'wirc_feii.txt',$
 'wirc_co.txt',$
 'wirc_ch4s.txt',$
 'wirc_ch4l.txt',$
 'CONICA_h.txt',$
 'CONICA_ks.txt',$
 'RSR-W1.txt',$
 'RSR-W2.txt',$
 'RSR-W3.txt',$
 'RSR-W4.txt']

; ---------------------------
; INFO
; ---------------------------
if (keyword_set(info) or n_params() lt 1) then begin
  print, 'FILTFLUX function'
  print, ''
  print, 'Usage:'
  print, ' IDL> mag = filtflux(wavelength, flux, filter [string name or number in set], units [string name or number in set], PROFILE=[filter transmission], /silent, /energy, /photons, /ab, /force, /mean)'
  print, ''
  print, 'Filters:'
  for i=1,n_elements(filtnames)-1 do print, strtrim(string(i))+' '+filtnames(i)
  print, ''
  print, 'Units:'
  for i=1,n_elements(unitsnames)-1 do print, strtrim(string(i))+' '+unitsnames(i)
  print, ''
  print, 'Examples:'
  print, '(1) If you have read in an array of wavelengths (in units of microns, preferred) and flux values (in units of erg/cm2/s/micron) into the variables wavelength and flux, than the following will determine the magnitude of the spectrum in the 2MASS J band (including optical response from the telescope and atmospheric absorption):'
  print, ''
  print, 'IDL> mag = filtflux(wavelength, flux, 6, 1)'
  print, ''
  print, 'Alternately, the same magnitude will be returned with: '
  print, ''
  print, "IDL> mag = filtflux(wavelength, flux, '2MASS J optical+atm resp (RSR)', 'erg/cm2/s/um')"  
  print, ''
  print, '(2) To compute the total energy emitted from a source in the Sloan r-band, where flux is measured in Jy, and also suppress comments, use:'
  print, ''
  print, "IDL> energy = filtflux(wavelength, flux, 'SDSS r', 'Jy', /silent, /energy)"
  print, ''
  print, '(3) If you read in your own filter profile, it can be fed to filtflux through the PROFILE keyword, as long as it is placed into an appropriate array format:'
  print, 'IDL> profile = [[filter_wavelength],[filter_transmission]'
  print, 'IDL> mag = filtflux(wavelength, flux, PROFILE=profile)'
  print, ''
  return, -1
endif


; ---------------------------
; SET UP FOR CALCULATION
; ---------------------------


; lam units - default is microns!!! 
; original in um
if (max(lambda) lt 1000.) then lam = lambda
; original in Ang
if (max(lambda) gt 1000. and min(lambda) lt 1.e12) then lam=lambda/1.e4 
; original in Hz
if (min(lambda) gt 1.e12) then lam = 3.e14/lambda

; filter choice processing
if (n_elements(profile) eq 0) then begin
 if (n_elements(filter) eq 0) then begin

  print, 'Filters:'
  for i=1,n_elements(filtnames)-1 do print, strtrim(string(i))+' '+filtnames(i)
  print, ''
  read, resp, prompt='Enter a filter number (1) or name: '
  if (resp ne '') then begin
   sz = size(resp,type=type)
   if (type eq 7) then begin
    filter=resp
    goto, FILTERSTR
   endif else filter=fix(resp) 
  endif else filter=1
 endif

 type = size(filter,/type)
 if (type eq 7) then begin
FILTERSTR: w = where(strlowcase(filter) eq strlowcase(filtnames),cnt)
  if (cnt lt 1) then goto, BADFILTER
  filter = w(0)
 endif
 sz = size(filter)
 if (sz(0) eq 0) then filter = [filter]  
 if (max(filter) ge n_elements(filtnames) or min(filter) le 0) then begin
BADFILTER: print, 'FILTFLUX: Attempting to use an undefined filter; EXITING'
  return, -99.
 endif
 if (keyword_set(silent) eq 0) then print, 'Using filter '+filtnames(filter)

; process custom filter
endif else begin
 sz = size(profile)
 if (sz(0) ne 2 or sz(2) lt 2) then begin
  print, 'FILTFLUX: Input filter profile array must be in the form [[wavelength],[transmission]]; EXITING'
  return, -99.
 endif
 filter = [0]
 if (keyword_set(silent) eq 0) then print, 'Using input profile'
endelse

; flux units - default is erg/cm2/s/um
if (n_elements(units) eq 0) then units=1
type = size(units,/type)
if (type eq 7) then begin
 w = where(strlowcase(units) eq strlowcase(unitsnames),cnt)
 if (cnt lt 1) then goto, BADUNITS
 units = w(0)
endif
if (units ge n_elements(unitsnames) or units le 0) then begin
BADUNITS: print, 'FILTFLUX: Attempting to use an undefined flux unit; EXITING'
 return, -99.
endif
if (keyword_set(silent) eq 0) then print, 'Using units of '+unitsnames(units)

; convert flux to erg/s/cm2/micron
case (units) of
  2: flx = flux*1.e4
  3: flx = flux*1.e3
  4: flx = flux*1.e7
  5: flx = flux*3.e14/(lam*lam)
  6: flx = flux*3.e17/(lam*lam)
  7: flx = flux*(3.e-9)/(lam*lam)
  else: flx=flux
endcase

; AB MAGNITUDE PREPARATION

if (keyword_set(ab) or keyword_set(mean)) then begin
 lamnu = 3.e14/lam
 flxnu = flx*lam*lam/3.e14
endif else $

; VEGA MAGNITUDE PREPARATION

; vega data - units converted to ergs/cm2/s/um 
 if (file_search(filtfolder+vfile) eq '') then begin
  print, 'FILTFLUX: Cannot find '+vfile+' in folder '+filtfolder+', please check; EXITING'
  return, -99.
 endif
 readcol, filtfolder+vfile, vlam, vflx, comment='#', format='f,f', /silent



; --------------------------------
; MAIN CALCULATION
; --------------------------------

; loop calculation for different filters
mag = filter*double(0)
for j=0,n_elements(filter)-1 do begin

; read in filter or use input profile
 if (n_elements(profile) eq 0) then begin
  if (file_search(filtfolder+filtfiles(filter(j))) eq '') then begin
   print, 'FILTFLUX: Cannot find filter file '+filtfiles(filter(j))+' in folder '+filtfolder+', please check; skipping measurement'
   rejflag = 0
  endif
  readcol, filtfolder+filtfiles(filter(j)), a, b, format='f,f', comment='#', /silent 
  filt = [[a],[b]]
 endif else filt = profile
 szf = size(filt)

 limit = [min(filt(*,0)),max(filt(*,0))]
 if ((limit(0) lt min(lam)) or (limit(1) gt max(lam))) then begin
   print, 'FILTFLUX: warning! filter profile extends beyond data values'
   w = where((filt(*,0) ge min(lam)) and (filt(*,0) le max(lam)), cnt)
   if (keyword_set(force) eq 0 or cnt lt 1) then begin
     rejflag = 1
     print, 'FILTFLUX: spectrum does not span filter; setting mag value to -99.'
   endif else begin
     rejflag = 0
     filt = filt(w,*)
     limit = [min(filt(*,0)),max(filt(*,0))]
   endelse
 endif

; check for uniq filter values
  u = uniq(filt(*,0))
  flam = filt(u,0)
  ftrans = filt(u,1)/max(filt(u,1))
  
  w = where((lam ge limit(0)) and (lam le limit(1)),cnt)

; weighting term for photon counts (except RSR filters)
  if (strpos(filtnames(j),'(RSR)') gt -1) then $
    flamwt = flam*0.+1. else flamwt = flam
  

; reject mag if necessary
  if (rejflag eq 1) then mag(j) = -99. else begin

;   w = append([w(0)-1],append(w,[w(cnt-1)+1]))
   w = [[w(0)-1],w,[w(cnt-1)+1]]
   flxcorr = interpol(flx(w),lam(w),flam)

; AB mags
   if (keyword_set(ab)) then begin
    if (keyword_set(silent) eq 0) then print, 'FILTFLUX: Computing AB magnitudes'
    flamnu = 3.e14/flam
    flxcorr = double(interpol(flxnu(w),lamnu(w),flamnu))
    a = int_tabulated(alog10(flamnu),flxcorr*ftrans,/double,/sort)
    b = int_tabulated(alog10(flamnu),ftrans,/double,/sort)
    mag(j) = -2.5*alog10(a/b)-48.6
    goto, ENDLOOP
   endif 

; mean flux density
   if (keyword_set(mean)) then begin
    if (keyword_set(silent) eq 0) then print, 'FILTFLUX: Computing mean flux density in erg/cm2/s/micron'
;   mag(j) = max(ftrans*flxcorr) 
    mag(j) = int_tabulated(flam,flamwt*ftrans*flxcorr,/double)/$
     int_tabulated(flam,flamwt*ftrans,/double) 
    goto, ENDLOOP
   endif

; energy
   if (keyword_set(energy)) then begin
    if (keyword_set(silent) eq 0) then print, 'FILTFLUX: Computing total energy in erg/cm2/s'
    mag(j) = int_tabulated(flam,ftrans*flxcorr,/double) 
    goto, ENDLOOP
   endif

; photons
   if (keyword_set(photons)) then begin
    if (keyword_set(silent) eq 0) then print, 'FILTFLUX: Computing total number of photons/s'
    mag(j) = 5.0277e11*int_tabulated(flam,ftrans*flxcorr*flam,/double)/int_tabulated(flam,ftrans*flam,/double)
    goto, ENDLOOP
   endif

; vega mags
   if (keyword_set(silent) eq 0) then print, 'FILTFLUX: Computing Vega magnitudes'
   w = where((vlam ge limit(0)) and (vlam le limit(1)),cnt)
   w = [[w(0)-1],w,[w(cnt-1)+1]]
   flxcorr2 = interpol(vflx(w),vlam(w),flam)
   mag(j) = -2.5*alog10(int_tabulated(flam,flamwt*ftrans*flxcorr,/double)/int_tabulated(flam,flamwt*ftrans*flxcorr2,/double))

ENDLOOP: dummy=0
 endelse 
endfor

; output results
if (n_elements(filter) eq 1) then return, mag(0) else return, mag


end
