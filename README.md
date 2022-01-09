# jovian_jrm09_internal
This repository contains IDL code used to calculate the JRM09 model (<a href="https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1002/2018GL077312
">Connerney et al., 2018</a>) of Jupiter's internal magnetic field model. This code can be used with <a href="https://github.com/marissav06/con2020_idl/blob/main/con2020_model_rtp.pro">con2020_model_rtp.pro</a>, which gives the field produced by the Connerney et al. (2020) current sheet. However, note that this JRM09 code currently only takes scalar inputs (i.e. inputs r_rj, colat_rads, elong_rads should be individual positions) while con2020_model_rtp.pro takes vector inputs.

<h3>Required inputs (System III Spherical, right handed):</h3>
<ul>
  <li>r_rj       - radial distance, in Rj. </li>                   
  <li>colat_rads - colatitude, in radians. Value(s) should be 0 <= colat_rads <=  pi. </li>
  <li>elong_rads - East longitude, right handed, in radians. Value(s) should be 0 <= elong_rads <= 2pi. </li>
</ul>

<h3>Outputs:</h3>
Spherical Magnetic field vector from the JRM09 internal magnetic field model, [Br, Btheta, Bphi], units of nT.

<h3>Usage:</h3>
<ul>
  <li>For internal field only: B = jovian_jrm09_internal(r_rj, colat_rads, elong_rads)</li>
 <li>For full field model: B = jovian_jrm09_internal(r_rj, colat_rads, elong_rads) + con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads)</li>
</ul>


This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson and was last updated December 2021
It is based on a routine originally written by K. Khurana and translated into IDL by Marissa Vogt in 2009 (for the VIP4 internal field model). 
Thanks to Masafumi Imai for providing code for his version of the JRM09 model, which was used to test and validate this code
