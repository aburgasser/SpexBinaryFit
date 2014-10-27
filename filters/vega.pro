pro vega

cd, 'C:\My Documents\tdwarfs\data\'

vega = readarr('vega.txt',/comment,/space,type='f')
lam = vega(*,0)*1.e-4
flx = vega(*,1)*1.e4

jfile = 'j_2mass_filter.txt'
hfile = 'h_2mass_filter.txt'
ksfile = 'ks_2mass_filter.txt'
ibfile= 'i_bessel.txt'
jbfile= 'j_bessel.txt'
hbfile= 'h_bessel.txt'
kbfile= 'k_bessel.txt'
lbfile= 'l_bessel.txt'
lpbfile= 'lp_bessel.txt'
mbfile= 'm_bessel.txt'

jfilt = readarr(jfile,/comment,/space,type='f')
hfilt = readarr(hfile,/comment,/space,type='f')
ksfilt = readarr(ksfile,/comment,/space,type='f')
ibfilt = readarr(ibfile,/comment,/space,type='f')
jbfilt = readarr(jbfile,/comment,/space,type='f')
hbfilt = readarr(hbfile,/comment,/space,type='f')
kbfilt = readarr(kbfile,/comment,/space,type='f')
lbfilt = readarr(lbfile,/comment,/space,type='f')
lpbfilt = readarr(lpbfile,/comment,/space,type='f')
mbfilt = readarr(mbfile,/comment,/space,type='f')

; jfilter
filt=jfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,1)-filt(i,1)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,1)),max(filt(*,1))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,1))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,2)*flxin(i)*dellam

jcnts = cnts

; h filter
filt=hfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,1)-filt(i,1)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,1)),max(filt(*,1))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,1))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,2)*flxin(i)*dellam

hcnts = cnts

; ks filter
filt=ksfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,1)-filt(i,1)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,1)),max(filt(*,1))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,1))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,2)*flxin(i)*dellam

kscnts = cnts

; I bess filter
filt=[[ibfilt(*,0)*1.e-4],[ibfilt(*,1)]]
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,0)-filt(i,0)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,0)),max(filt(*,0))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,0))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,1)*flxin(i)*dellam

ibcnts = cnts

; J bess filter
filt=[[jbfilt(*,0)*1.e-4],[jbfilt(*,1)]]
;filt=jbfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,0)-filt(i,0)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,0)),max(filt(*,0))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,0))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,1)*flxin(i)*dellam

jbcnts = cnts

; H bess filter
filt=[[hbfilt(*,0)*1.e-4],[hbfilt(*,1)]]
;filt=hbfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,0)-filt(i,0)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,0)),max(filt(*,0))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,0))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,1)*flxin(i)*dellam

hbcnts = cnts

; K bess filter
filt=[[kbfilt(*,0)*1.e-4],[kbfilt(*,1)]]
;filt=kbfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,0)-filt(i,0)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,0)),max(filt(*,0))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,0))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,1)*flxin(i)*dellam

kbcnts = cnts

; L bess filter
filt=[[lbfilt(*,0)*1.e-4],[lbfilt(*,1)]]
;filt=lbfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,0)-filt(i,0)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,0)),max(filt(*,0))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,0))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,1)*flxin(i)*dellam

lbcnts = cnts

; Lp bess filter
filt=[[lpbfilt(*,0)*1.e-4],[lpbfilt(*,1)]]
;filt=lpbfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,0)-filt(i,0)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,0)),max(filt(*,0))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,0))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,1)*flxin(i)*dellam

lpbcnts = cnts

; M bess filter
filt=[[mbfilt(*,0)*1.e-4],[mbfilt(*,1)]]
;filt=mbfilt
szf = size(filt)
dum = fltarr(szf(1)-1)
for i=0,szf(1)-2 do dum(i)=filt(i+1,0)-filt(i,0)
mom= moment(dum)
dellam = mom(0)
limit = [min(filt(*,0)),max(filt(*,0))]

w = where((lam ge limit(0)) and (lam le limit(1)))
flxin = interpol(flx(w),lam(w),filt(*,0))
cnts = 0.
for i=0,szf(1)-1 do cnts = cnts+filt(i,1)*flxin(i)*dellam

mbcnts = cnts

; plot it

w = where((lam ge min(ibfilt(*,0)*1e-4)) and (lam le max(mbfilt(*,0)*1.e-4)))
mx = max(flx(w))
plotps, 'vega.ps', /open
plot, [-1.e3],[-1.e3],xrange=1.e-4*[min(ibfilt(*,0))*0.8,max(mbfilt(*,0))*1.1],$
  yrange=[-1.e-2*mx,mx*1.05],$
  xstyle=1,ystyle=1,xtitle='!4k!3 (!4l!3m)',ytitle='F!D!4k!3!N (erg cm!U-2!N s!U-1!N !4l!3m!U-1!N)',$
  title='Vega'
oplot, lam, flx
oplot, jfilt(*,1),jfilt(*,2)*mx,linestyle=1
oplot, hfilt(*,1),hfilt(*,2)*mx,linestyle=1
oplot, ksfilt(*,1),ksfilt(*,1)*mx,linestyle=1
oplot, ibfilt(*,0)*1.e-4,ibfilt(*,1)*mx
oplot, jbfilt(*,0)*1.e-4,jbfilt(*,1)*mx
oplot, hbfilt(*,0)*1.e-4,hbfilt(*,1)*mx
oplot, kbfilt(*,0)*1.e-4,kbfilt(*,1)*mx
oplot, lbfilt(*,0)*1.e-4,lbfilt(*,1)*mx
oplot, lpbfilt(*,0)*1.e-4,lpbfilt(*,1)*mx
oplot, mbfilt(*,0)*1.e-4,mbfilt(*,1)*mx
xyouts, median(jfilt(*,1)), 0.5*mx, 'J band!C'+strtrim(string(jcnts),2), align=0.5,charsize=0.7
xyouts, median(hfilt(*,1)), 0.45*mx, 'H band!C'+strtrim(string(hcnts),2), align=0.5,charsize=0.7
xyouts, median(ksfilt(*,1)), 0.5*mx, 'Ks band!C'+strtrim(string(kscnts),2), align=0.5,charsize=0.7
xyouts, median(ibfilt(*,0)*1.e-4), 0.3*mx, 'I band!C'+strtrim(string(ibcnts),2), align=0.5,charsize=0.7
xyouts, median(jbfilt(*,0)*1.e-4), 0.2*mx, 'J band!C'+strtrim(string(jbcnts),2), align=0.5,charsize=0.7
xyouts, median(hbfilt(*,0)*1.e-4), 0.3*mx, 'H band!C'+strtrim(string(hbcnts),2), align=0.5,charsize=0.7
xyouts, median(kbfilt(*,0)*1.e-4), 0.2*mx, 'K band!C'+strtrim(string(kbcnts),2), align=0.5,charsize=0.7
xyouts, median(lbfilt(*,0)*1.e-4), 0.3*mx, 'L band!C'+strtrim(string(lbcnts),2), align=0.5,charsize=0.7
xyouts, median(lpbfilt(*,0)*1.e-4), 0.4*mx, 'Lp band!C'+strtrim(string(lpbcnts),2), align=0.5,charsize=0.7
xyouts, median(mbfilt(*,0)*1.e-4), 0.3*mx, 'M band!C'+strtrim(string(mbcnts),2), align=0.5,charsize=0.7
plotps, /close

cnts = [jcnts,hcnts,kscnts,ibcnts,jbcnts,hbcnts,kbcnts,lbcnts,lpbcnts,mbcnts]
for i=0,9 do print, cnts(i)

return
end