;-----------------------------------------------------------------------------------------------------
FUNCTION spexbinaryfit_dateconvert, $
;
; Purpose: convert numerical date code into a YYYY MMM DD string
; 
; HERSTORY
; created: 2013 Sep 11 by Adam J. Burgasser
; -----------------------------------------------------------------------------------------------------
instring, long=long

monstr = ['January','February','March','April','May','June','July','August','September','October','November','December']

if (keyword_set(long)) then mstr = monstr else mstr = strmid(monstr,0,3)

case 1 of
 strmid(instring,0,0) eq '0': return, '20'+strmid(instring,0,2)+' '+mstr(fix(strmid(instring,2,2))-1)+' '+strmid(instring,4,2)
 strmid(instring,0,0) eq '9': return, '19'+strmid(instring,0,2)+' '+mstr(fix(strmid(instring,2,2))-1)+' '+strmid(instring,4,2)
 else: return, strmid(instring,0,4)+' '+mstr(fix(strmid(instring,4,2))-1)+' '+strmid(instring,6,2)
endcase

end


