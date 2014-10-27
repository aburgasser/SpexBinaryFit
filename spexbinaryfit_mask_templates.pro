; -----------------------------------------------------------------------------------------------------
; PRO SPEXBINARYFIT_MASK_TEMPLATES
;
; Purpose: create a new structure with templates properly masked
; 
; history: 
;  2014 Oct 27 (AJB): initial file
; -----------------------------------------------------------------------------------------------------
function spexbinaryfit_mask_templates, db, mask

on_error, 0
special_tags = ['NTEMPLATES','WAVELENGTH','FLUX','NOISE','MAGS','FILTERS','FILTER_NAMES']

if (n_elements(mask) eq 0) then return, db
tags = strupcase(tag_names(db))

w = where(mask eq 0,cnt)
if (cnt eq 0) then message, 'Warning, all templates rejected (all mask = 1)' else begin

; create with specific tags
 db_out = create_struct($
  'ntemplates',cnt,$
  'wavelength',db.wavelength,$
  'flux',db.flux(*,w),$
  'noise',db.flux(*,w),$
  'mags',db.mags(*,w),$
  'filters',db.filters,$
  'filter_names',db.filter_names)
  
; general tags
 for i=0,n_elements(tags)-1 do begin
  wt = where(special_tags eq tags(i),cnt)
  if (cnt eq 0) then db_out = create_struct(db_out,tags(i),db.(i)[w])
 endfor
; specifc tags
 
 return, db_out
endelse

end

