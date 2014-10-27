;-----------------------------------------------------------------------------------------------------
; PROCEDURE SPEXBINARYFIT_NIRCLASSIFY
;
; Purpose: classify SpeX spectra of L and T dwarfs
; 
; created: 8 Aug 2010
; -----------------------------------------------------------------------------------------------------

pro spexbinaryfit_nirclassify, lam, flx, spt, spt_e, spts, indices, noise=noise, sptinit=sptinit, sptind=sptind

; SPECTRAL TYPE/INDEX COEFFICIENTS FROM BURGASSER 2007
sptcoeff = [$
 [19.4914  ,   -39.1927   ,   131.159   ,  -215.619  ,    103.831],$  ; H2O-J
 [20.9821  ,   -19.7818   ,   25.2719   ,  -32.2132  ,   0.908676],$ ; CH4-J
 [27.0751  ,   -84.5044    ,  242.363   ,  -338.112  ,    149.137],$ ; H2O-H
 [20.1266  ,   -22.9054    ,  43.6100   ,  -50.6812  ,    20.8428],$ ; CH4-J
 [47.0508  ,   -221.852    ,  626.776   ,  -804.704  ,    353.621],$ ; H2O-K
 [18.8457  ,   -22.4606    ,  25.3379   ,  -4.73417  ,   -12.5855]] ; CH4-K

; M & L DWARF SPECTRAL TYPE/INDEX COEFFICIENTS FROM ????
sptcoeff2 = [$
 [23.9868, -25.8775],$  ; H2O-J sig = 0.94
 [23.5419, -27.5808],$  ; H2O-H sig = 0.86
 [24.1948, -27.1723]]  ; H2O-K sig = 0.86

; COMPUTE INDICES 
 indices = [$
   specindex(lam,flx,[[1.14,1.165],[1.26,1.285]],0,[0,1],/integr), $	; H2O-J
   specindex(lam,flx,[[1.315,1.335],[1.26,1.285]],0,[0,1],/integr), $	; CH4-J
   specindex(lam,flx,[[1.48,1.52],[1.56,1.60]],0,[0,1],/integr), $	; H2O-H
   specindex(lam,flx,[[1.635,1.675],[1.56,1.60]],0,[0,1],/integr), $	; CH4-H
   specindex(lam,flx,[[1.975,1.995],[2.08,2.1]],0,[0,1],/integr), $  	; H2O-K
   specindex(lam,flx,[[2.215,2.255],[2.08,2.12]],0,[0,1],/integr),$ 	; CH4-K
   specindex(lam,flx,[[2.06,2.10],[1.25,1.29]],0,[0,1],/integr),$	 	; K/J
   specindex(lam,flx,[[1.27,1.30],[1.30,1.33]],0,[0,1],/integr),$ 	; J-shape
   specindex(lam,flx,[[1.61,1.64],[1.56,1.59],[1.66,1.69]],0,[0,1,1],/integr),$ 	; H-dip
   specindex(lam,flx,[[2.06,2.10],[2.10,2.14]],0,[0,1],/integr)]	 	; K-shape

; ALTERNATE CLASSIFICATIONS

 typs = [$
  poly(indices(0),sptcoeff(*,0)),$
  poly(indices(1),sptcoeff(*,1)),$
  poly(indices(2),sptcoeff(*,2)),$
  poly(indices(3),sptcoeff(*,3)),$
  poly(indices(4),sptcoeff(*,4)),$
  poly(indices(5),sptcoeff(*,5)),$
  -119.,-119.,-119.,-119.]+20.
 typs2 = [$
  poly(indices(0),sptcoeff2(*,0)),$
  poly(indices(1),sptcoeff(*,1)),$
  poly(indices(2),sptcoeff2(*,1)),$
  poly(indices(3),sptcoeff(*,3)),$
  poly(indices(4),sptcoeff2(*,2)),$
  poly(indices(5),sptcoeff(*,5)),$
  -119.,-119.,-119.,-119.]+20.

 if (n_elements(sptinit) eq 0) then spt = mean(typs([0,1,2,3,5])) else spt = sptinit

; loop through
 for j=0,2 do begin
  case(1) of
   spt ge 29: begin		; T dwarfs: use H2O-J, CH4-J, H2O-H, CH4-H, CH4-K
    wtyp = [0,1,2,3,5]
   end
   spt le 22: begin		; M/L dwarfs: use alternate H2O-J, H2O-H, H2O-K
    wtyp = [0,2,4]
    typs = typs2
   end
   else: begin
    wtyp = [0,2,5]			; L dwarfs: use H2O-J, H2O-H, H2O-K
   end
  endcase
 endfor

spt = mean(typs(wtyp))
spt_e = stddev(typs(wtyp))
spts = spexbinaryfit_numericspt(spt,/reverse)

sptind = typs
x = intarr(n_elements(typs))
x(wtyp) = 1
w = where(x eq 0)
sptind(w) = -99.

if (spt_e ge 1) then spts = spts+':'
if (spt_e ge 2) then spts = spts+':'

return
end
