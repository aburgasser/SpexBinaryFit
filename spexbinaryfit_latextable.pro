;-----------------------------------------------------------------------------------------------------
pro SPEXBINARYFIT_LATEXTABLE, $ 
;
; Purpose: create a latex table summarizing information from spex sources
; 
; HERSTORY
; created: 2013 Sep 11 by Adam J. Burgasser
; -----------------------------------------------------------------------------------------------------
indb, outfile, index=index, title=title, newobs=newobs, standalone=standalone, allrefs=allrefs

if (n_elements(indb) eq 0) then message, 'no db'
if (n_elements(outfile) eq 0) then outfile='/Users/adam/idl/spexbinaryfit/templates/latex_table.tex'
if (n_elements(title) eq 0) then title='SpeX Spectra'
snrflag = 0

; make sure required tags are present, add dummies as necessary
db = indb
tnames = strlowcase(tag_names(db))
if (max(strpos(tnames,'name')) lt 0) then message, 'no name variable in db'
if (max(strpos(tnames,'designation')) lt 0) then $
 db = create_struct(db,'designation',strarr(n_elements(db.name)))
if (max(strpos(tnames,'opt_type')) lt 0) then $
 db = create_struct(db,'opt_type',strarr(n_elements(db.name)))
if (max(strpos(tnames,'spex_type')) lt 0) then $
 db = create_struct(db,'spex_type',strarr(n_elements(db.name)))
if (max(strpos(tnames,'jmag')) lt 0) then $
 db = create_struct(db,'jmag',strarr(n_elements(db.name)))
if (max(strpos(tnames,'kmag')) lt 0) then $
 db = create_struct(db,'kmag',strarr(n_elements(db.name)))
if (max(strpos(tnames,'observation_date')) lt 0) then $
 db = create_struct(db,'observation_date',strarr(n_elements(db.name)))
if (max(strpos(tnames,'discovery_reference')) lt 0) then $
 db = create_struct(db,'discovery_reference',strarr(n_elements(db.name)))
if (max(strpos(tnames,'optical_reference')) lt 0) then $
 db = create_struct(db,'optical_reference',strarr(n_elements(db.name)))
if (max(strpos(tnames,'nir_reference')) lt 0) then $
 db = create_struct(db,'nir_reference',strarr(n_elements(db.name)))
if (max(strpos(tnames,'data_reference')) lt 0) then $
 db = create_struct(db,'data_reference',strarr(n_elements(db.name)))
if (max(strpos(tnames,'snr')) lt 0) then $
 db = create_struct(db,'snr',strtrim(string(intarr(n_elements(db.name))),2))
w = where(db.data_reference eq '',cnt)
if (cnt gt 0) then db.data_reference(w) = 'UNPUB'

; index for table lines
if (n_elements(index) eq 0) then index = indgen(n_elements(db.name))

; create mapping of references
if keyword_set(allrefs) then $
 refs = [db.discovery_reference(index),db.optical_reference(index),db.nir_reference(index),db.data_reference(index)] else $
 refs = db.data_reference(index)
if (n_elements(refs) gt 0) then begin
 refs = strcompress(refs,/remove_all)
 w = where(refs eq '')
 pop, refs, w, /index
 s = sort(refs) 
 u = uniq(refs(s))
 refs = refs(s(u))
endif

openw, unit, outfile, /get_lun

if (keyword_set(standalone)) then begin
 printf, unit, '\documentclass[12pt,preprint]{aastex}'
 printf, unit, '\begin{document}
 printf, unit, ''
endif

printf, unit, '% '+strtrim(string(n_elements(index)),2)+' sources'

if (keyword_set(newobs)) then $
 printf, unit, '\begin{deluxetable}{lllllcccccl}' else $
 printf, unit, '\begin{deluxetable}{lllllccl}'
printf, unit, '\tablecaption{'+title+'\label{tab:spex}}
printf, unit, '\tabletypesize{\scriptsize}
printf, unit, '\tablewidth{0pt}
printf, unit, '\rotate
printf, unit, '\tablehead{
printf, unit, ' & & \multicolumn{3}{c}{Spectral Type} & \multicolumn{2}{c}{2MASS} \\'
printf, unit, '\cline{3-5} \cline{6-7}'
printf, unit, '\colhead{Source} &
printf, unit, '\colhead{Designation} &
printf, unit, '\colhead{Opt} &
printf, unit, '\colhead{NIR} &
printf, unit, '\colhead{SpeX\tablenotemark{a}} &
printf, unit, '\colhead{$J$} &
printf, unit, '\colhead{$J-K_s$} &
if (keyword_set(newobs)) then begin
 printf, unit, '\colhead{Date} &'
 printf, unit, '\colhead{$\lambda/\Delta\lambda$} &'
 printf, unit, '\colhead{SNR} &'
endif
printf, unit, '\colhead{Ref\tablenotemark{b}} \\
printf, unit, '}
printf, unit, '\startdata

for i=0,n_elements(index)-1 do begin
 if (i gt 0 and db.name(index(i)) eq db.name(index((i-1)>0)) and db.designation(index(i)) eq db.designation(index((i-1)>0))) then $
  line=' & & & & & & & ' else begin
  line = strreplace(strreplace(db.name(index(i)),'[',''),']','')+' & '+db.designation(index(i))+' & '
  if (db.opt_type(index(i)) eq '') then line=line+'\nodata'
  line = line+db.opt_type(index(i))+' & '
  if (db.nir_type(index(i)) eq '') then line=line+'\nodata'
  line = line+db.nir_type(index(i))+' & '
  if (db.spex_type(index(i)) eq '') then line=line+'\nodata'
  line=line+db.spex_type(index(i))+' & '
  if (db.jmag(index(i)) ne '') then $
   line=line+strmid(strtrim(string(rndoff(float(db.jmag(index(i)))*100.)/100.),2),0,5)+' & ' else $
   line=line+'\nodata & '
  if (db.jmag(index(i)) ne '' and db.kmag(index(i)) ne '') then begin
   jk = float(db.jmag(index(i)))-float(db.kmag(index(i)))
   if (jk lt 0) then len = 5 else len=4
   line=line+strmid(strtrim(string(rndoff(jk*100.)/100.),2),0,len)+' & '
  endif else line=line+'\nodata & '
 endelse
 if (keyword_set(newobs)) then begin
  line=line+spexbinaryfit_dateconvert(db.observation_date(index(i)))+' & '+db.resolution(index(i))+' & '+strtrim(string(fix(db.snr(index(i)))),2)
  if fix(db.snr(index(i))) eq 0 then begin
   line=line+'\tablenotemark{c}'
   snrflag = 1
  endif
  line=line+' & '
 endif

; references
 if n_elements(refs) gt 0 then begin
  if keyword_set(allrefs) then begin
   if (db.discovery_reference(index(i)) ne '') then line=line+strtrim(string(where(refs eq db.discovery_reference(index(i)))+1),2) else line=line+'\nodata'
   line=line+'; '
   if (db.optical_reference(index(i)) eq '' and db.nir_reference(index(i)) eq '') then line=line+'\nodata' else begin
    if (db.optical_reference(index(i)) ne '') then line=line+strtrim(string(where(refs eq db.optical_reference(index(i)))+1),2)
    if (db.optical_reference(index(i)) ne '' and db.nir_reference(index(i)) ne '' and db.nir_reference(index(i)) ne db.optical_reference(index(i))) then line=line+', '
    if (db.nir_reference(index(i)) ne '' and db.nir_reference(index(i)) ne db.optical_reference(index(i))) then line=line+strtrim(string(where(refs eq db.nir_reference(index(i)))+1),2)
   endelse
   line=line+'; '
  endif
  line=line+strtrim(string(where(refs eq db.data_reference(index(i)))+1),2)
 endif
 printf, unit, line+' \\'
; print, line
endfor

printf, unit, '\enddata
printf, unit, '\tablenotetext{a}{Near-infrared classification from SpeX data based on index method described in \citet{2007ApJ...659..655B}.}
if (keyword_set(allrefs)) then $
 printf, unit, '\tablenotetext{b}{First citation is the discovery reference; next citation(s) are classification references (optical and near-infrared); final citation is the data reference.} else $
 printf, unit, '\tablenotetext{b}{Citation for data.}
if snrflag eq 1 then printf, unit, '\tablenotetext{c}{Original data did not contain uncertainty array so no signal-to-noise ratio could be calculated.}
line='\tablerefs{'
for i=0,n_elements(refs)-1 do line=line+'('+strtrim(string(i+1),2)+') \citet{'+refs(i)+'}; '
printf, unit, strmid(line,0,strlen(line)-2)+'}'
printf, unit, '\end{deluxetable}

if (keyword_set(standalone)) then $
 printf, unit, '\end{document}

close, unit
free_lun, unit

FINISH: return
end
