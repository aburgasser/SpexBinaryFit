; -----------------------------------------------------------------------------------------------------
; PRO SPEXBINARYFIT_PRINTOUT
;
; Purpose: printout data from fits
; 
; HISTORY
; 2010 Nov 29: updated for clearer outputs
; 2012 Jan 28: added "full" keyword to prevent excess text on screen
; 2014 Feb 11: added latex outputs
; -----------------------------------------------------------------------------------------------------

pro spexbinaryfit_printout, db, dev, outfile=outfile, nbest=nbest, delim=delim, short=short, single=single, dof=dof, prob=prob, full=full, mko=mko, nicmos=nicmos, wfc3=wfc3, nirc2=nirc2, latex=latex

on_error, 0

tb='	'
pmstr = '+/-'

if (n_elements(nbest) eq 0) then nbest=10
if (keyword_set(latex)) then begin
 delim=' & '
 pmstr='$\pm$'
endif
if (n_elements(delim) eq 0) then delim=tb
if (n_elements(dof) eq 0) then dof=100.
if (n_elements(prob) eq 0) then prob='N/A'


; check if it has mass as a unit - for model fitting
tagnames = tag_names(db)
w = where(tagnames eq 'MASS',cnt)
if (cnt gt 0) then begin
 physstr = ['Teff','logg','Mass','Age','Lbol','Radius']
 sigdig = [0,2,2,1,2,2]
 physparam = [[db.teff],[db.logg],[db.mass],[alog10(db.age)+9.],[db.lbol],[db.radius]]
 pflag = 1 
endif else pflag = 0

; SET UP FILTER REPORTING
case 1 of
 keyword_set(nicmos): wfilt = where($
  db.filter_names eq 'NIC1 F090M' or $
  db.filter_names eq 'NIC1 F110W' or $
  db.filter_names eq 'NIC1 F145M' or $
  db.filter_names eq 'NIC1 F160W' or $
  db.filter_names eq 'NIC1 F165M' or $
  db.filter_names eq 'NIC1 F170M',fcnt)
 keyword_set(wfc3): wfilt = where($
  db.filter_names eq 'NIC1 F110W' or $
  db.filter_names eq 'WFC3 F127M' or $
  db.filter_names eq 'WFC3 F139M' or $
  db.filter_names eq 'NIC1 F165M' or $
  db.filter_names eq 'WFC3 F164N',fcnt)
 keyword_set(mko): wfilt = where($
  db.filter_names eq 'MKO J' or $
  db.filter_names eq 'MKO H' or $
  db.filter_names eq 'MKO K',fcnt)
 keyword_set(nirc2): wfilt = where($
  db.filter_names eq 'NIRC2 J' or $
  db.filter_names eq 'NIRC2 H' or $
;  db.filter_names eq 'NIRC2 K' or $
  db.filter_names eq 'NIRC2 Ks' ,fcnt)
 keyword_set(full): wfilt = where($
  db.filter_names ne '', fcnt)
 else:  wfilt=where($
  db.filter_names eq '2MASS J' or $
  db.filter_names eq '2MASS H' or $
  db.filter_names eq '2MASS Ks' or $
  db.filter_names eq 'MKO J' or $
  db.filter_names eq 'MKO H' or $
  db.filter_names eq 'MKO K',fcnt)
endcase

; INDIVIDUAL FITS

s = sort(dev(*,0))

; HEADER
if (n_elements(outfile) ne 0) then begin
  openw, unit, outfile, /get_lun
  if (keyword_set(latex)) then begin
    if (keyword_set(single)) then begin
      printf, unit, '\begin{deluxetable}{llc}'
      printf, unit, '\tablecaption{Single Template Fits}'
      printf, unit, '\tablewidth{0pt}'
      printf, unit, '\tablehead{'
      printf, unit, '\colhead{Primary} &'
      printf, unit, '\colhead{SpT} &'
      printf, unit, '\colhead{$\chi^2$} \\'
      printf, unit, '}'
      printf, unit, '\startdata'
    endif else begin
      ccode = 'llll'
      for i=0,fcnt-1 do ccode=ccode+'c'
      printf, unit, '\begin{deluxetable}{'+ccode+'c}'
      printf, unit, '\tablecaption{Binary Template Fits}'
      printf, unit, '\tablewidth{0pt}'
      printf, unit, '\tablehead{'
      printf, unit, '\colhead{Primary} &'
      printf, unit, '\colhead{SpT} &'
      printf, unit, '\colhead{Secondary} &'
      printf, unit, '\colhead{SpT} &'
      for i=0,fcnt-1 do printf, unit, '\colhead{$\Delta$'+db.filter_names(wfilt(i))+'} &'
      printf, unit, '\colhead{$\chi^2$} \\'
      printf, unit, '}'
      printf, unit, '\startdata'
   endelse
 endif 
endif

line = '# Primary'+delim+'SpT'
if (keyword_set(single) eq 0) then begin
;   if (pflag) then for i=0,n_elements(physstr)-1 do line=line+delim+physstr(i)+'_A'
  line = line+delim+'Secondary'+delim+'SpT'
;  if (pflag) then for i=0,n_elements(physstr)-1 do line=line+delim+physstr(i)+'_B'
  for j=0,fcnt-1 do line=line+delim+'d'+db.filter_names(wfilt(j))
endif ; else begin
; if (pflag) then for i=0,n_elements(physstr)-1 do line=line+delim+physstr(i)
;endelse
if (n_elements(outfile) eq 0) then print, line else begin
 if (keyword_set(latex) eq 0) then printf, unit, line+delim+'chi2 \\'
endelse


for i=0,nbest-1 do begin 
; for binaries
 if (keyword_set(single) eq 0) then begin
; primary
  if (keyword_set(short) eq 1) then line= db.primaries.shortname(db.pairs(s(i),0))+delim else line= db.primaries.name(db.pairs(s(i),0))+delim
  line=line+db.primaries.refspt(db.pairs(s(i),0))+delim
;  if (pflag) then begin
;   for j=0,n_elements(physstr)-1 do begin
;    pstr = strtrim(string(rndoff(physparam(pairs(s(i),0),j)*10^sigdig(j))/10.^sigdig(j)),2)
;    line=line+delim+strmid(pstr,0,strpos(pstr,'.')+sigdig(j)+1)
;   endfor
;  endif
; secondary
  if (keyword_set(short) eq 1) then line=line+db.secondaries.shortname(db.pairs(s(i),1))+delim else line=line+db.secondaries.name(db.pairs(s(i),1))+delim
  line=line+db.secondaries.refspt(db.pairs(s(i),1))
;  if (pflag) then begin
;   for j=0,n_elements(physstr)-1 do begin
;    pstr = strtrim(string(rndoff(physparam(pairs(s(i),1),j)*10^sigdig(j))/10.^sigdig(j)),2)
;    line=line+delim+strmid(pstr,0,strpos(pstr,'.')+sigdig(j)+1)
;   endfor
;  endif
  for j=0,fcnt-1 do line=line+delim+strmid(strtrim(string(rndoff(db.relmags(wfilt(j),s(i))*1000.)/1000.),2),0,5)
  line=line+delim+strmid(strtrim(string(rndoff(dev(s(i),0)*100.)/100.),2),0,4)

; for singles
 endif else begin
  if (keyword_set(short) eq 1) then line= db.shortname(s(i))+delim else line= db.name(s(i))+delim
  line=line+db.refspt(s(i))+delim
;  if (pflag) then begin
;   for j=0,n_elements(physstr)-1 do begin
;    pstr = strtrim(string(rndoff(physparam(s(i),j)*10^sigdig(j))/10.^sigdig(j)),2)
;    line=line+delim+strmid(pstr,0,strpos(pstr,'.')+sigdig(j)+1)
;   endfor
;  endif
  line=line+strmid(strtrim(string(rndoff(dev(s(i),0)*100.)/100.),2),0,4)
 endelse
 
 if (n_elements(outfile) eq 0) then print, '>> '+line else begin
   if (keyword_set(latex)) then line=line+' \\'
   printf, unit, line
 endelse
endfor

; SUMMARY OF BEST FITS

eta = dev(*,0)/min(dev(*,0),loc)
wt = 2.*(1.-f_pdf(eta,dof,dof))
w = where(wt gt 0.1,cnt)


; for binaries
if (keyword_set(single) eq 0) then begin
 pspt = total(wt(w)*db.primaries.refsptn(db.pairs(w,0)))/total(wt(w))
 pspt_e = sqrt(total(wt(w)*(db.primaries.refsptn(db.pairs(w,0))-pspt)^2)/total(wt(w))) 
 sspt = total(wt(w)* db.secondaries.refsptn(db.pairs(w,1)))/total(wt(w))
 sspt_e = sqrt(total(wt(w)*(db.secondaries.refsptn(db.pairs(w,1))-sspt)^2)/total(wt(w)))
 dmagmn = fltarr(fcnt,2)
 for j=0,fcnt-1 do begin
  dm = total(wt(w)*(db.relmags(wfilt(j),w)))/total(wt(w))
  dm_e = sqrt(total(wt(w)*(db.relmags(wfilt(j),w)-dm)^2)/total(wt(w)))
  dmagmn(j,*) = [dm,dm_e] 
 endfor
; if (pflag) then begin
;  pphysmn = fltarr(n_elements(physstr),2)
;  sphysmn = fltarr(n_elements(physstr),2)
;  for j=0,n_elements(physstr)-1 do begin
;   dm = total(wt(w)*(physparam(db.pairs(w,0),j)))/total(wt(w))
;   dm_e = sqrt(total(wt(w)*(physparam(db.pairs(w,0),j)-dm)^2)/total(wt(w)))
;   pphysmn(j,*) = [dm,dm_e] 
;   dm = total(wt(w)*(physparam(db.pairs(w,1),j)))/total(wt(w))
;   dm_e = sqrt(total(wt(w)*(physparam(db.pairs(w,1),j)-dm)^2)/total(wt(w)))
;   sphysmn(j,*) = [dm,dm_e] 
;  endfor
; endif

; for singles
endif else begin
 pspt = total(wt(w)*db.refsptn(w))/total(wt(w))
 pspt_e = sqrt(total(wt(w)*(db.refsptn(w)-pspt)^2)/total(wt(w)))  
; if (pflag) then begin
;  pphysmn = fltarr(n_elements(physstr),2)
;  for j=0,n_elements(physstr)-1 do begin
;   dm = total(wt(w)*(physparam(w,j)))/total(wt(w))
;   dm_e = sqrt(total(wt(w)*(physparam(w,j)-dm)^2)/total(wt(w)))
;   pphysmn(j,*) = [dm,dm_e] 
;  endfor
; endif
endelse

; primary
line = 'Primary'+delim+$
 spexbinaryfit_numericspt(pspt, /reverse)+pmstr+strmid(strtrim(string(rndoff(pspt_e*10.)/10.),2),0,3)
; if (pflag) then begin
;  for j=0,n_elements(physstr)-1 do begin
;   pstr = strtrim(string(rndoff(pphysmn(j,0)*10^sigdig(j))/10.^sigdig(j)),2)
;   pstr_e = strtrim(string(rndoff(pphysmn(j,1)*10^sigdig(j))/10.^sigdig(j)),2)
;   line=line+delim+strmid(pstr,0,strpos(pstr,'.')+sigdig(j)+1)+pmstr+strmid(pstr_e,0,strpos(pstr_e,'.')+sigdig(j)+1)
;  endfor
; endif
 line=line+delim
 
; secondary
if (keyword_set(single) eq 0) then begin
 line = line+'Secondary'+delim+$
 spexbinaryfit_numericspt(sspt, /reverse)+pmstr+strmid(strtrim(string(rndoff(sspt_e*10.)/10.),2),0,3)
; if (pflag) then begin
;  for j=0,n_elements(physstr)-1 do begin
;   pstr = strtrim(string(rndoff(sphysmn(j,0)*10^sigdig(j))/10.^sigdig(j)),2)
;   pstr_e = strtrim(string(rndoff(sphysmn(j,1)*10^sigdig(j))/10.^sigdig(j)),2)
;   line=line+delim+strmid(pstr,0,strpos(pstr,'.')+sigdig(j)+1)+pmstr+strmid(pstr_e,0,strpos(pstr_e,'.')+sigdig(j)+1)
;  endfor
; endif
 for j=0,n_elements(dmagmn(*,0))-1 do line=line+delim+strmid(strtrim(string(rndoff(dmagmn(j,0)*1000.)/1000.),2),0,5)+pmstr+strmid(strtrim(string(rndoff(dmagmn(j,1)*1000.)/1000.),2),0,5)
 line=line+delim+strtrim(string(prob),2)
endif 

if (n_elements(outfile) eq 0) then print, '>> '+line else begin
   if (keyword_set(latex)) then line=line+' \\'
   printf, unit, line
endelse

if (n_elements(outfile) ne 0 and keyword_set(latex) ne 0) then begin
  printf, unit, '\enddata'
  printf, unit, '\end{deluxetable}'
endif

if (n_elements(outfile) ne 0) then begin
 close, unit
 free_lun, unit
endif
   
return
end
