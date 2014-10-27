; -------------------------------------------------------------------------------
; PRO DESIGNATION2SHORTNAME
;
; Purpose: Convert source designations into "short names"; e.g., 1234+3456
; -------------------------------------------------------------------------------

function spexbinaryfit_designation2shortname, designation, shortname

on_error, 0

; first strip off up to J
pos_start = strpos(designation,'J',/reverse_search)+1

; identify where +/- lie
pos_mid = strpos(designation,'+',/reverse_search)
w = where(pos_mid eq -1,cnt)
if (cnt gt 0) then pos_mid(w) = strpos(designation(w),'-',/reverse_search)
w = where(pos_mid eq -1,cnt)
if (cnt gt 0) then pos_mid(w) = 0

; create shortnames
shortname = strarr(n_elements(designation))
for i=long(0),n_elements(designation)-1 do shortname(i) = strmid(designation(i),pos_start(i),4)+strmid(designation(i),pos_mid(i),5)

FINISH: return, shortname
end

