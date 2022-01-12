FUNCTION jovian_jrm09_internal_xyz, x_rj, y_rj, z_rj
  ;% Code to calculate the JRM09 model of Jupiter's internal magnetic field model
  ;% Reference: Connerney et al., 2018, https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2018GL077312
  ;%
  ;% Required inputs (System III Spherical, right handed):
  ;%  x_rj       - Jupiter SYSIII right-handed position in x, in Rj.
  ;%  y_rj       - Jupiter SYSIII right-handed position in y, in Rj.
  ;%  z_rj       - Jupiter SYSIII right-handed position in z, in Rj.
  ;%
  ;% Outputs:
  ;%  B - Cartesian Magnetic field vector the JRM09 internal magnetic field model, [Bx, By, Bz], units of nT.
  ;%
  ;% Usage:
  ;% For internal field only: B = jovian_jrm09_internal_xyz(x_rj, y_rj, z_rj)
  ;%
  ;% This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu).
  ;% Last updated January 2022.
  ;% It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009.
  ;% Thanks to Masafumi Imai for providing code for his version of the JRM09 model, which was used to test and validate this code.
  ;% Thanks to Gabby Provan, Matt James, and Marty Brennan for helpful discussions.

  ON_ERROR, 2 ; % Exit code if an error in main, don't stop in code - no Matlab equivalent, just delete line in Matlab

  ;% check inputs are same size
  N_input = N_ELEMENTS(x_rj)
  scalar_input = (N_input EQ 1) ;% scalar or not

  IF (N_input NE N_ELEMENTS(y_rj)) THEN MESSAGE,'ERROR: First argument x_rj must be the same size as 2nd argument y_rj'
  IF (N_input NE N_ELEMENTS(z_rj)) THEN MESSAGE,'ERROR: First argument x_rj must be the same size as 3rd argument z_rj'
  IF (ISA(x_rj, NUMBER=1) EQ 0) OR (SIZE(x_rj, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: First  argument x_rj must be a scalar number or 1D column array of numbers'
  IF (ISA(y_rj, NUMBER=1) EQ 0) OR (SIZE(y_rj, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Second argument y_rj must be a scalar number or 1D column array of numbers'
  IF (ISA(z_rj, NUMBER=1) EQ 0) OR (SIZE(z_rj, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Third  argument z_rj must be a scalar number or 1D column array of numbers'

  ;% Covert to Doubles, and rename input variables so as not to over-write them (only an IDL issue)
  x_in = double(x_rj); % X in SYSIII, units Rj
  y_in = double(y_rj); % Y in SYSIII, units Rj
  z_in = double(z_rj); % Z in SYSIII, units Rj

  rho_rj_sq = x_in*x_in + y_in*y_in
  r_rj = sqrt(rho_rj_sq + z_in*z_in)

  colat_rads = acos(z_in/r_rj)
  elong_rads = atan(y_in,x_in)
  IF scalar_input THEN BEGIN
    IF (elong_rads LT 0)      THEN elong_rads      += 2d*!DPI
  ENDIF ELSE BEGIN
    ind = WHERE(elong_rads LT 0, NULL = 1)
    IF (N_ELEMENTS(ind) NE 0) THEN elong_rads[ind] += 2d*!DPI
  ENDELSE

  ;% Do calculation in spherical
  brtp = jovian_jrm09_internal_rtp(r_rj, colat_rads, elong_rads);

  ;% r_rj already set up earlier
  rho_rj = sqrt(rho_rj_sq)
  cos_phi   = x_in/rho_rj
  sin_phi   = y_in/rho_rj
  cos_theta = z_in  /r_rj
  sin_theta = rho_rj/r_rj


  ;% Convert to cartesian coordinates
  bxyz = [ $
    [brtp[*,0]*sin_theta*cos_phi + brtp[*,1]*cos_theta*cos_phi - brtp[*,2]*sin_phi], $ ;% bx
    [brtp[*,0]*sin_theta*sin_phi + brtp[*,1]*cos_theta*sin_phi + brtp[*,2]*cos_phi], $ ;% by
    [brtp[*,0]*cos_theta         - brtp[*,1]*sin_theta                            ]  $ ;% bz
    ];% size n x 3

  RETURN,bxyz
END