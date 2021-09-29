;% ======
;% jovian_jrm09_internal
;% ======
;% Code to calculate the JRM09 model of Jupiter's internal magnetic field model 
;% Reference: Connerney et al., 2018, https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2018GL077312
;%
;% Required inputs (System III Spherical, right handed):
;%  r_rj       - radial distance, in Rj.                    
;%  colat_rads - colatitude, in radians.                    Value(s) should be 0 <= colat_rads <=  pi.
;%  elong_rads - East longitude, right handed, in radians.  Value(s) should be 0 <= elong_rads <= 2pi.
;%
;% Outputs:
;%  B - Spherical Magnetic field vector the JRM09 internal magnetic field model, [Br, Btheta, Bphi], units of nT.
;%
;% This code can be used with con2020_model_rtp.pro (https://github.com/marissav06/con2020_idl/blob/main/con2020_model_rtp.pro),
;%     which gives the field produced by the Connerney et al. (2020) current sheet, to yield the full (internal+external) magnetic field in
;%     Jupiter's magnetosphere
;%
;% Usage:
;% For internal field only: B = jovian_jrm09_internal(r_rj, colat_rads, elong_rads)
;% For full field model: B = jovian_jrm09_internal(r_rj, colat_rads, elong_rads) + con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads)
;% However, note that this code currently only takes scalar inputs (i.e. r_rj, colat_rads, elong_rads should be individual points) while con2020_model_rtp takes vector inputs
;%
;% This code was written by Marissa Vogt (mvogt@bu.edu) and last updated Sept. 2021
;% It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009
;% Thanks to Masafumi Imai for providing code for his version of the JRM09 model, which was used to test and validate this code

FUNCTION jovian_jrm09_internal, r_rj, colat_rads, elong_rads

  ;Values from Connerney et al. 2020
  ;See supplemental online information Table S1, https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2018GL077312&file=grl57087-sup-0005-2018GL077312-ds01.txt

  g = [0.d, 0.d, 410244.7d, -71498.3d, 11670.4d, -56835.8d, 48689.5d, $
    4018.6d, -37791.1d, 15926.3d, -2710.5d, $
    -34645.4d, -8247.6d, -2406.1d, -11083.8d, -17837.2d, $
    -18023.6d, 4683.9d, 16160.0d, -16402.0d, -2600.7d, -3660.7d, $
    -20819.6d, 9992.9d, 11791.8d, -12574.7d, 2669.7d, 1113.2d, 7584.9d, $
    598.4d, 4665.9d, -6495.7d, -2516.5d, -6448.5d, 1855.3d, -2892.9d, 2968.d, $
    10059.2d, 1934.4d, -6702.9d, 153.7d, -4124.2d, -867.2d, -3740.6d, -732.4d, -2433.2d, $
    9671.8d, -3046.2d, 260.9d, 2071.3d, 3329.6d, -2523.1d, 1787.1d, -1148.2d, 1276.5d, -1976.8d, $
    -2299.5d, 2009.7d, 2127.8d, 3498.3d, 2967.6d, 16.3d, 1806.5d, -46.5d, 2897.8d, 574.5d, 1298.9d, replicate(0.d, 200)]

  h = [0.d, 0.d, 0.d, 21330.5d, 0.d, -42027.3d, 19353.2d, 0.d, -32957.3d, 42084.5d, -27544.2d, $
    0.d, 31994.5d, 27811.2d, -926.1d, 367.1d, 0.d, 45347.9d, -749.0d, 6268.5d, 10859.6d, 9608.4d,$
    0.d, 14533.1d, -10592.9d, 568.6d, 12871.7d, -4147.8d, 3604.4d, $
    0.d, -7626.3d, -10948.4d, 2633.3d, 5394.2d, -6050.8d, -1526.d, -5684.2d, $
    0.d, -2409.7d, -11614.6d, 9287.0d, -911.9d, 2754.5d, -2446.1d, 1207.3d, -2887.3d, $
    0.d, -8467.4d, -1383.8d, 5697.7d, -2056.3d, 3081.5d, -721.2d, 1352.5d, -210.1d, 1567.6d, $
    0.d, -4692.6d, 4445.8d, -2378.6d, -2204.3d, 164.1d, -1361.6d, -2031.5d, 1411.8d, -714.3d, 1676.5d, replicate(0.d, 200)]

knm=18
nm = 10
a = dindgen(18)*0.d
b = dindgen(18)*0.d
rec = dindgen(172)*0.d

for n=1, 18 do begin
  n2 = 2*n-1;
  n2 = n2*(n2-2);
  for m=1, n do begin
    mn=n*(n-1)/2+m;
    rec(mn)=double((n-m)*(n+m-2.d))/n2;
  endfor
endfor
s = 1.d;
for n=2, 18 do begin
  mn=n*(n-1)/2+1.d;
  s=s*double(2.d*n-3.d)/double(n-1);
  g(mn)=g(mn)*s;
  h(mn)=h(mn)*s;
  p=s;
  for m=2, n do begin
    aa = 1.d;
    if m eq 2 then aa = 2.d
    p=p*sqrt(aa*double(n-m+1)/double(n+m-2));
    mnn=mn+m-1;
    g(mnn)=g(mnn)*p;
    h(mnn)=h(mnn)*p;
  endfor
endfor

if knm ne nm then begin
  knm = nm
  k = knm+1
endif
pp = 1.0d/r_rj;
p = pp;
for n=1, k do begin
  p = p*pp;
  a(n) = p;
  b(n) = p*n;
endfor

p = 1.d
d = 0.d
bbr = 0.d
bbt = 0.d
bbf = 0.d
cos_phi = cos(elong_rads)
sin_phi = sin(elong_rads)
cos_theta = cos(colat_rads)
sin_theta = sin(colat_rads)
bk = sin_theta lt 1e-5

for m=1,k do begin
  bm = m eq 1
  if bm eq 1 then begin
    x = 0.d
    y = 1.d
  endif else begin
    mm = m-1
    w = x
    x = w*cos_phi + y*sin_phi
    y = y*cos_phi - w*sin_phi
  endelse
  q = p
  z = d
  bi = 0.d
  p2 = 0.d
  d2 = 0.d
  for n=m, k do begin
    an = a(n)
    mn = n*(n-1)/2.d + m
    e = g(mn)
    hh = h(mn)
    w = e*y+hh*x
    if abs(p2) lt 1e-38 then p2 = 0.d
    if abs(q) lt 1e-38 then q = 0.d
    bbr = bbr + b(n)*w*q
    bbt = bbt - an*w*z
    if (bm ne 1.d) then begin
      qq = q
      if bk eq 1.d then qq=z
      bi = bi+an*(e*x-hh*y)*qq
    endif
    xk=rec(mn)
    dp = cos_theta*z - sin_theta*q - xk*d2
    pm = cos_theta*q - xk*p2
    d2 = z
    p2 = q
    z = dp
    q = pm
  endfor
  d = sin_theta*d + cos_theta*p
  p = sin_theta*p
  if (bm ne 1.d) then begin
    bi = bi*mm
    bbf = bbf + bi
  endif
endfor

br = bbr
bt  = bbt
if bk eq 1 then begin
  if (cos_theta lt 0.d) then bbf = -bbf
  bf = bbf
endif
if bk ne 1 then bf = bbf/sin_theta

fieldarray = [[br],[bt],[bf]]
return,fieldarray



END
