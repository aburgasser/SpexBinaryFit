; -----------------------------------------------------------------------------
; FUNCTION READSTRUCT
;
; Purpose: Read in tab-delimited database file to string array
; output is a structure containing header keywords and parameters
; -----------------------------------------------------------------------------

function spexbinaryfit_readstruct, file_in, delim=delim

on_error, 1

if (n_elements(delim) eq 0) then delim='	'

if (file_search(file_in) eq '') then begin
 message, 'Database file '+file_in+' not present'
 goto, FINISH
endif

arr = readarr(file_in,delim=delim,max=50000)
header = reform(arr(0,*))
; replace bad symbols
header = strrep(header,'(','')
header = strrep(header,')','')
header = strrep(header,'#','')
header = strrep(header,'&','')
header = strrep(header,'?','')
header = strrep(header,'-','')
header = strrep(header,';','')
header = strrep(header,'"','')
header = strtrim(header,2)

nrows = n_elements(arr(*,0))

output_structure = create_struct(header(0), strrep(arr(1:nrows-1,0),'"',''))

for i=long(1),n_elements(header)-1 do output_structure = create_struct(output_structure, header(i), strrep(arr(1:nrows-1,i),'"',''))

FINISH: return, output_structure
end

