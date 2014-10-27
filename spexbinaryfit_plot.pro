; -----------------------------------------------------------------------------------------------------
; PRO SPEXBINARYFIT_PLOT
;
; Purpose: plot fits in a uniform manner
; 
; updates: 
; 2009 Nov 18 - written by AJB
; 2012 Jun 12 - AJB fixed plotting bug
; 2014 Feb 11 - AJB fixed formatting of characters
; 2014 Oct 27 - AJB
; 	added showdiff keyword to show difference spectrum
; 	made sure noise spectrum was actually plotted
; -----------------------------------------------------------------------------------------------------

pro spexbinaryfit_plot, lambda, flux1, flux2, flux3, flux4, noise=noise, bw=bw, labels=labels, pcolors=pcolors, lcolors=lcolors, outfile=outfile, xtitle=xtitle, ytitle=ytitle, yrange=yrange, xrange=xrange, right=right, inset=inset, note=note, plotfitrange=plotfitrange, showdiff=showdiff

!p.font=0
!p.thick=4
!x.thick=3
!y.thick=3
!p.multi=0
csty = '!3'

if (n_elements(xtitle) eq 0) then xtitle='Wavelength (!9m'+csty+'m)'
if (n_elements(ytitle) eq 0) then ytitle='Normalized F!D!9l!7!N'
if (n_elements(xrange) lt 2) then xrange =[0.9,2.4]
if (n_elements(yrange) eq 0) then yrange =[-0.02,1.45]
if (keyword_set(showdiff)) then yrange(0) = -0.1
if (n_elements(pcolors) eq 0) then pcolors=[0,6,4,8,0,0]
if (n_elements(lcolors) eq 0) then lcolors=pcolors
lcolors = [lcolors[0:n_elements(labels)-2],0]
noisecolor=150
diffcolor = 100

if (n_elements(outfile) gt 0) then begin
  set_plot, 'ps'
  device, /encapsulated, xsize=24, ysize=18, filename=outfile,/portrait, bits_per_pixel=8, /color
endif

tek_color

plot, [0],[0], xrange=xrange,yrange= yrange,/xsty,/ysty, xtitle=csty+xtitle, ytitle=csty+ytitle, charsize=1.8

; indicate fitting range
if (n_elements(plotfitrange) ge 2) then begin
 for i=0,n_elements(plotfitrange(0,*))-1 do polyfill, [plotfitrange(0,i), plotfitrange(0,i), plotfitrange(1,i), plotfitrange(1,i), plotfitrange(0,i)],[0.98,1.0,1.0,0.98,0.98]*(yrange(1)-yrange(0))+yrange(0), color=200
endif

oplot, lambda, flux1, color=pcolors(0), thick=4
 if (n_elements(flux2) gt 0) then begin
  oplot, lambda, flux2, color=pcolors(1), thick=4
  diff = flux1-flux2
 endif
 if (n_elements(flux3) gt 0) then begin
  oplot, lambda, flux3, color=pcolors(2), thick=4
  diff = flux1-flux3
 endif
 if (n_elements(flux4) gt 0) then begin
  oplot, lambda, flux4, color=pcolors(3), thick=4
  diff = flux1-flux4
 endif
if (keyword_set(showdiff)) then oplot, lambda, diff, color=diffcolor
if (n_elements(noise) eq n_elements(lambda)) then oplot, lambda, noise, color=noisecolor
oplot, [0,10],[0,0], linestyle=1

if (n_elements(labels) gt 0) then begin
 if (keyword_set(right)) then $
  for i=0,n_elements(labels)-1 do xyouts, xrange(1)-0.05*(xrange(1)-xrange(0)), yrange(1)-(i+2)*0.05*(yrange(1)-yrange(0)), csty+labels(i), align=1, charsize=1.4, color=lcolors(i) else $
  for i=0,n_elements(labels)-1 do xyouts, xrange(0)+0.05*(xrange(1)-xrange(0)), yrange(1)-(i+2)*0.05*(yrange(1)-yrange(0)), csty+labels(i), align=0, charsize=1.4, color=lcolors(i)
endif
 
if (keyword_set(inset)) then begin
 w = where(lambda gt 1.55 and lambda lt 1.7)
 med = median(flux1(w))
 std = stddev(flux1(w))
 xrng = [1.5,1.74]
 yrng = [-1,1]*3.*std+med
  plot, [0],[0], xrange=xrng, yrange=yrng,/xsty,/ysty, xtitle=csty+xtitle, charsize=0.9, xmargin=[75,10], ymargin=[37,7], /noerase
  polyfill, [xrng(0),xrng(0),xrng(1),xrng(1),xrng(0)], [yrng(0),yrng(1),yrng(1),yrng(0),yrng(0)], color=1
 oplot, lambda, flux1, color=pcolors(0), thick=4
 if (n_elements(flux2) gt 0) then begin
  oplot, lambda, flux2, color=pcolors(1), thick=4
  diff = flux1-flux2
 endif
 if (n_elements(flux3) gt 0) then begin
  oplot, lambda, flux3, color=pcolors(2), thick=4
  diff = flux1-flux3
 endif
 if (n_elements(flux4) gt 0) then begin
  oplot, lambda, flux4, color=pcolors(3), thick=4
  diff = flux1-flux4
 endif
endif
  
if (n_elements(outfile) gt 0) then begin
  device, /close
  set_plot, 'x'
endif

return
end

