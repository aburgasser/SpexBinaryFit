function spexbinaryfit_numericspt, spt, reverse=reverse, roundoff=roundoff, flags=flags

on_error, 0

if (n_elements(spt) eq 1) then spt = [spt]
if (n_elements(flags) eq 0) then flags = ''
if (n_elements(flags) eq 1) then flags = [flags]

letters = ['K','M','L','T','Y']
offsets = [0.,10.,20.,30.,40.]
pclasses = ['d/sd','sd','esd','usd']
sclasses = ['subdwarf','giant','young','alpha','beta','gamma']

; NUMBER -> LETTER
if keyword_set(reverse) then begin
 if keyword_set(roundoff) then spt2 = rndoff(spt*2.)/2. else spt2=spt
 sptn = strarr(n_elements(spt))+'UNK'
 for i=0,n_elements(letters)-1 do begin
  w = where(spt2 ge offsets(i) and spt2 lt offsets(i)+10.,cnt)
  if (cnt gt 0) then sptn(w) = letters(i)+strmid(strtrim(string(spt2(w)-offsets(i)),2),0,3)
 endfor
; luminosity or subdwarf class
 for i=0,n_elements(pclasses)-1 do begin
  w = where(flags eq pclasses(i),cnt)
  if (cnt gt 0) then sptn(w) = pclasses(i)+sptn(w)
 endfor
 for i=0,n_elements(sclasses)-1 do begin
  w = where(flags eq sclasses(i),cnt)
  if (cnt gt 0) then sptn(w) = sptn(w)+' '+sclasses(i)
 endfor

; LETTER -> NUMBER
endif else begin

 spt_work = strcompress(spt,/remove_all)

; initiallize
 sptn = fltarr(n_elements(spt))*0.-99.

; remove all extraneum
 for i=0,n_elements(spt_work)-1 do begin 
  spt_work(i) = strreplace(spt_work(i),'V','')
  spt_work(i) = strreplace(spt_work(i),'I','')
  for j=0,n_elements(pclasses)-1 do spt_work(i) = strreplace(spt_work(i),pclasses(j),'')
  for j=0,n_elements(sclasses)-1 do spt_work(i) = strreplace(spt_work(i),sclasses(j),'')
 endfor
 
; fill in based on first letter
 for i=0,n_elements(letters)-1 do begin
  pos = strpos(strupcase(spt_work),letters(i))
  w = where(pos gt -1,cnt)
  if (cnt gt 0) then begin 
   for j=0,cnt-1 do sptn(w(j)) = float(strmid(spt_work(w(j)),pos(w(j))+1,3))+offsets(i)
  endif
 endfor

; luminosity or subdwarf class
 for i=0,n_elements(pclasses)-1 do begin
  pos = strpos(strupcase(spt),pclasses(i))
  w = where(pos gt -1,cnt)
  if (cnt gt 0) then flags(w) = pclasses(i)
 endfor
 for i=0,n_elements(sclasses)-1 do begin
  pos = strpos(strupcase(spt),sclasses(i))
  w = where(pos gt -1,cnt)
  if (cnt gt 0) then flags(w) = sclasses(i)
 endfor

; w = where(strmid(spt_work,0,1) eq 'L',cnt)
; if (cnt gt 0) then sptn(w) = float(strmid(spt_work(w),1,3))+20.
; w = where(strmid(spt_work,0,1) eq 'T',cnt)
; if (cnt gt 0) then sptn(w) = float(strmid(spt_work(w),1,3))+30.

endelse

FINISH: return, sptn
end
