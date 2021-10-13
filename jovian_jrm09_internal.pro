FUNCTION jovian_jrm09_internal, r_rj, colat_rads, elong_rads
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
  ;%
  ;% This code was written by Marissa Vogt (mvogt@bu.edu) and Rob Wilson (rob.wilson@lasp.colorado.edu)
  ;% Last updated October 2021
  ;% It is based on a routine originally written by K. Khurana, translated into IDL by Marissa Vogt in 2009
  ;% Thanks to Masafumi Imai for providing code for his version of the JRM09 model, which was used to test and validate this code

  ON_ERROR, 2 ; % Exit code if an error in main, don't stop in code - no Matlab equivalent, just delete line in Matlab

  N_input = N_ELEMENTS(r_rj)
  scalar_input = (N_input EQ 1) ;% scalar or not

  ;% Check inputs r_rj, colat_rads and elong_rads are all numbers, and same size (also scalar or 1D only)
  IF (ISA(r_rj      , NUMBER=1) EQ 0) OR (SIZE(r_rj      , N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: First  argument    r_rj    must be a scalar number or 1D array of numbers';
  IF (ISA(colat_rads, NUMBER=1) EQ 0) OR (SIZE(colat_rads, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Second argument colat_rads must be a scalar number or 1D array of numbers';
  IF (ISA(elong_rads, NUMBER=1) EQ 0) OR (SIZE(elong_rads, N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Third  argument elong_rads must be a scalar number or 1D array of numbers';
  IF (N_input NE N_ELEMENTS(colat_rads)) THEN MESSAGE,'ERROR: First argument r_rj must be the same size as 2nd argument colat_rads';
  IF (N_input NE N_ELEMENTS(elong_rads)) THEN MESSAGE,'ERROR: First argument r_rj must be the same size as 3rd argument elong_rads';
  IF scalar_input THEN BEGIN
    IF (    r_rj   LE 0d) OR (    r_rj   GE 200d   ) THEN MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'
    IF (colat_rads LT 0d) OR (colat_rads GT !DPI   ) THEN MESSAGE,'ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'
    IF (elong_rads LT 0d) OR (elong_rads GT !DPI*2d) THEN MESSAGE,'ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'
  ENDIF ELSE BEGIN
    min_x = MIN(r_rj, MAX=max_x)
    min_y = MIN(colat_rads, MAX=max_y)
    min_z = MIN(elong_rads, MAX=max_z)
    IF (min_x LE 0d) OR (max_x GE 200d   ) THEN MESSAGE,'ERROR: First  argument, Position    r_rj   , must be in units of Rj and >0 and <200 only, and not outside that range (did you use km instead?)'
    IF (min_x LT 0d) OR (max_y GT !DPI   ) THEN MESSAGE,'ERROR: Second argument, Position colat_rads, must be in units of radians and >= 0 and <=   Pi only, and not outside that range (did you use degrees instead?)'
    IF (min_z LT 0d) OR (max_z GT !DPI*2d) THEN MESSAGE,'ERROR: Third  argument, Position elong_rads, must be in units of radians and >= 0 and <= 2*Pi only, and not outside that range (did you use degress instead?)'
  ENDELSE

  ;============
  ;Begin hard-coding for JRM09
  ;Values from Connerney et al. 2020
  ;See supplemental online information Table S1, https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2018GL077312&file=grl57087-sup-0005-2018GL077312-ds01.txt
  ;============
  ;% Used function _jovian_jrm09_internal_rjw_setup (below in this file) to copy/paste rec, g and h
  ;% rec, g and h never change - so why calculate them each time?
  degree = 20
  order = 10
  rec = [$
    0d,                                0d,        0.33333333333333331482962d,                                0d, $
    0.26666666666666666296592d,        0.20000000000000001110223d,                                0d,        0.25714285714285711748062d, $
    0.22857142857142856429142d,        0.14285714285714284921269d,                                0d,        0.25396825396825395415590d, $
    0.23809523809523808202115d,        0.19047619047619046561692d,        0.11111111111111110494321d,                                0d, $
    0.25252525252525254151337d,        0.24242424242424243097105d,        0.21212121212121212709967d,        0.16161616161616162989922d, $
    0.09090909090909091161414d,                                0d,        0.25174825174825177231952d,        0.24475524475524476630817d, $
    0.22377622377622377602968d,        0.18881118881118880148406d,        0.13986013986013987042689d,        0.07692307692307692734701d, $
    0d,        0.25128205128205127749652d,        0.24615384615384616751044d,        0.23076923076923078204103d, $
    0.20512820512820512108831d,        0.16923076923076924016343d,        0.12307692307692308375522d,        0.06666666666666666574148d, $
    0d,        0.25098039215686274161499d,        0.24705882352941177515504d,        0.23529411764705882026405d, $
    0.21568627450980393245317d,        0.18823529411764705621124d,        0.15294117647058824704942d,        0.10980392156862744945656d, $
    0.05882352941176470506601d,                                0d,        0.25077399380804954454049d,        0.24767801857585139413409d, $
    0.23839009287925697067045d,        0.22291021671826624639401d,        0.20123839009287924906033d,        0.17337461300309597866942d, $
    0.13931888544891640746570d,        0.09907120743034056320475d,        0.05263157894736841813099d,                                0d, $
    0.25062656641604008633806d,        0.24812030075187968547468d,        0.24060150375939848288454d,        0.22807017543859647856763d, $
    0.21052631578947367252397d,        0.18796992481203006475354d,        0.16040100250626565525636d,        0.12781954887218044403241d, $
    0.09022556390977443108170d,        0.04761904761904761640423d,                                0d,        0.25051759834368531043580d, $
    0.24844720496894409644817d,        0.24223602484472050999642d,        0.23188405797101449556941d,        0.21739130434782608092270d, $
    0.19875776397515526605630d,        0.17598343685300207872579d,        0.14906832298136646342002d,        0.11801242236024844789455d, $
    0.08281573498964803214939d,        0.04347826086956521618454d,                                0d,        0.25043478260869567186830d, $
    0.24869565217391303990624d,        0.24347826086956522728677d,        0.23478260869565217849875d,        0.22260869565217392129775d, $
    0.20695652173913042792819d,        0.18782608695652172614565d,        0.16521739130434781595014d,        0.13913043478260869734164d, $
    0.10956521739130434256460d,        0.07652173913043477937457d,        0.04000000000000000083267d,                                0d, $
    0.25037037037037035425158d,        0.24888888888888888173412d,        0.24444444444444443642617d,        0.23703703703703704608330d, $
    0.22666666666666665519436d,        0.21333333333333334702608d,        0.19703703703703703831174d,        0.17777777777777778456247d, $
    0.15555555555555555802272d,        0.13037037037037035869247d,        0.10222222222222222820509d,        0.07111111111111111104943d, $
    0.03703703703703703498107d,                                0d,        0.25031928480204340692339d,        0.24904214559386972371868d, $
    0.24521072796934864634899d,        0.23882503192848020256989d,        0.22988505747126436462580d,        0.21839080459770116027229d, $
    0.20434227330779056175381d,        0.18773946360153256907033d,        0.16858237547892720997744d,        0.14687100893997445671957d, $
    0.12260536398467432317450d,        0.09578544061302682322001d,        0.06641123882503192910054d,        0.03448275862068965469387d, $
    0d,        0.25027808676307006230388d,        0.24916573971078975757720d,        0.24582869855394884339717d, $
    0.24026696329254726425262d,        0.23248053392658510341029d,        0.22246941045606227760345d,        0.21023359288097887009883d, $
    0.19577308120133482538527d,        0.17908787541713014346278d,        0.16017797552836485208694d,        0.13904338153503892350216d, $
    0.11568409343715238546402d,        0.09010011123470522409473d,        0.06229143492769743939430d,        0.03225806451612903136272d, $
    0d,        0.25024437927663734093642d,        0.24926686217008797719075d,        0.24633431085043988595373d, $
    0.24144672531769306722538d,        0.23460410557184752100568d,        0.22580645161290321953906d,        0.21505376344086021833668d, $
    0.20234604105571846188738d,        0.18768328445747800570231d,        0.17106549364613879427033d,        0.15249266862170088310258d, $
    0.13196480938416421668791d,        0.10948191593352883665968d,        0.08504398826979471526233d,        0.05865102639296188025142d, $
    0.03030303030303030387138d,                                0d,        0.25021645021645022577417d,        0.24935064935064935043307d, $
    0.24675324675324675216537d,        0.24242424242424243097105d,        0.23636363636363635909454d,        0.22857142857142856429142d, $
    0.21904761904761904656169d,        0.20779220779220780590535d,        0.19480519480519481456682d,        0.18008658008658010030167d, $
    0.16363636363636363535434d,        0.14545454545454544748040d,        0.12554112554112553667984d,        0.10389610389610390295267d, $
    0.08051948051948051854332d,        0.05541125541125541120735d,        0.02857142857142857053643d]
  g   = [$
    0d,                                0d,   410244.70000000001164153218269d,   -71498.30000000000291038304567d, $
    17505.59999999999854480847716d,   -98442.49328882319969125092030d,    42166.34389756242308067157865d,    10046.50000000000000000000000d, $
    -115711.13977311669441405683756d,    30841.14733335159326088614762d,    -2142.83839947159822258981876d,  -151573.62500000000000000000000d, $
    -45642.10215250826877309009433d,    -9415.35553115892798814456910d,   -23183.43100524596593459136784d,   -13190.78728838805909617803991d, $
    -141935.84999999997671693563461d,    47619.25007516491314163431525d,   124193.04328343032102566212416d,   -77191.29987306696421001106501d, $
    -5769.73075946518656564876437d,    -2568.20347420563030027551576d,  -300582.97499999997671693563461d,   188897.03523126235813833773136d, $
    176219.39807558275060728192329d,  -125279.49167958697944413870573d,    14568.18469249509325891267508d,     2590.20913175944133399752900d, $
    5094.72643062895076582208276d,    16044.59999999999854480847716d,   165497.62303578440332785248756d,  -188120.73349108296679332852364d, $
    -51533.85619684254197636619210d,   -79632.08074599411338567733765d,    11455.48572598611099238041788d,    -7006.09637449074671167181805d, $
    1921.06723268604082477395423d,   505710.56250000005820766091347d,   129665.25000000000000000000000d,  -375914.50046967255184426903725d, $
    6366.18842053944717918056995d,  -110265.51694787663291208446026d,   -12861.08441189933364512398839d,   -25680.06972210412277490831912d, $
    -1835.99981426163253672712017d,    -1524.90263109687475662212819d,   918443.19531249988358467817307d,  -388096.44079238711856305599213d, $
    28346.79585738653622684068978d,   171882.35647868152591399848461d,   187708.65848237561294808983803d,   -85005.63010931370081380009651d, $
    31091.84173000367445638403296d,    -8649.99331942896606051363051d,     3298.44757452672547515248880d,    -1203.96883845257525535998866d, $
    -414889.08398437500000000000000d,   488932.02253022126387804746628d,   448310.26802277535898610949516d,   578200.21502919984050095081329d, $
    346825.93308208207599818706512d,     1204.82452975374712877965067d,    74644.73375998696428723633289d,     -932.00811829940028019336751d, $
    23711.52712466476441477425396d,     1525.17384009177953885227907d,      771.06330156869501024630154d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d]
  h   = [$
    0d,                                0d,                                0d,    21330.50000000000000000000000d, $
    0d,   -72793.41890493947721552103758d,    16760.36284452099789632484317d,                                0d, $
    -100910.71037478504877071827650d,    81496.28381673301919363439083d,   -21775.60208180246991105377674d,                                0d, $
    177057.11204695011838339269161d,   108828.53403772377350833266973d,    -1937.07712643301852040167432d,      271.47411104698363715215237d, $
    0d,   461033.11139938322594389319420d,    -5756.22459277780399133916944d,    29500.89399184978901757858694d, $
    24092.34750470571452751755714d,     6740.87640657726115023251623d,                                0d,   274721.00218349619535729289055d, $
    -158302.75800766979227773845196d,     5664.86031229477885062806308d,    70239.09162317455047741532326d,    -9651.15831540766339458059520d, $
    2421.05129224630400130990893d,                                0d,  -270501.83727851061848923563957d,  -317074.53216031729243695735931d, $
    53925.73158082475129049271345d,    66612.60292471760476473718882d,   -37360.45546854781423462554812d,    -3695.70433387703633343335241d, $
    -3679.15443532142626281711273d,                                0d,  -161525.20312500000000000000000d,  -651374.26444599486421793699265d, $
    384663.57749869779217988252640d,   -24380.75866950407362310215831d,    40850.84987612628174247220159d,   -16793.03281485293700825423002d, $
    3026.49177465601997027988546d,    -1809.49012278727877855999395d,                                0d, -1078776.11541115446016192436218d, $
    -150349.92758701223647221922874d,   472811.32743136369390413165092d,  -115925.43081370404979679733515d,   103818.65529778851487208157778d, $
    -12547.38753045641169592272490d,    10189.09246170325241109821945d,     -542.89372143209163823485142d,      954.74582717435089307400631d, $
    0d, -1141644.22994741331785917282104d,   936694.13928736466914415359497d,  -393135.81781678373226895928383d, $
    -257618.41363149805692955851555d,    12129.55247439201775705441833d,   -56261.42789238762634340673685d,   -40717.73101774691895116120577d, $
    11552.18924515208345837891102d,    -1896.31274843787309691833798d,      995.21720307946509365137899d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d,                                0d, $
    0d,                                0d,                                0d]
  ;============
  ;End parts that are hard-coded for JRM09
  ;============

  ;% Changing inputs to Doubles, and not using input names (so as not to alter inputs, an IDL issue)
  r_rj_dbl       = DOUBLE(    r_rj  )
  colat_rads_dbl = DOUBLE(colat_rads)
  elong_rads_dbl = DOUBLE(elong_rads)

  ;% RJW: This code can be sped up by assuming k = 11 always, and replacing k with 11, Should we?
  IF degree NE order THEN BEGIN
    degree = order
    k = degree + 1
  ENDIF
  IF scalar_input THEN BEGIN
    a         = [0d,0d,0d,0d,0d,0d,0d,0d,0d, 0d, 0d, 0d] ;% = DBLARR(k+1)
    DINDGEN_k = [0d,1d,2d,3d,4d,5d,6d,7d,8d, 9d,10d,11d] ;% = DINDGEN(k+1), done manually for speed
  ENDIF ELSE BEGIN
    a           = DBLARR(N_input,k+1)
    DINDGEN_k   = a
    FOR i = 0,k DO DINDGEN_k[*,i] = i
  ENDELSE

  ; Instead of: a = (1d/r_rj_dbl)^DINDGEN_kp1, do the da to end of for loop
  da = 1d/r_rj_dbl
  IF scalar_input THEN BEGIN
    a[0] = da
    FOR i=1,k DO a[i] = a[i-1]*da
  ENDIF ELSE BEGIN ; the following vectorized 2 lines works for scalars, but is slower.
    a[*,0] = da
    FOR i=1,k DO a[*,i] = a[*,i-1]*da
  ENDELSE

  b = a * DINDGEN_k

  cos_phi   = cos(elong_rads_dbl)
  sin_phi   = sin(elong_rads_dbl)
  cos_theta = cos(colat_rads_dbl)
  sin_theta = sin(colat_rads_dbl)
  not_bk = (sin_theta GE 0.00001d); % = 1d-5 - also see bk both times below
  IF scalar_input THEN BEGIN
    ;%bk = (sin_theta LT 0.00001d); bk not needed for scalar
    zero_array = 0d
    p   = 1d
    d   = 0d
    bbr = 0d
    bbt = 0d
    bbf = 0d
    x = 0d
    y = 1d
  ENDIF ELSE BEGIN
    bk = (sin_theta LT 0.00001d)
    zero_array = DBLARR(N_input)
    p   = zero_array + 1d
    d   = zero_array
    bbr = zero_array
    bbt = zero_array
    bbf = zero_array
    x = zero_array;% 0s
    y = p;% 1s
  ENDELSE

  FOR m = 1, order DO BEGIN
    bm  = (m NE 1)
    IF bm THEN BEGIN
      m_minus_1 = DOUBLE(m - 1)
      w = x
      x = w*cos_phi + y*sin_phi
      y = y*cos_phi - w*sin_phi
    ENDIF
    q = p
    z = d
    bi = zero_array
    p2 = zero_array
    d2 = zero_array
    FOR n = m, k DO BEGIN
      mn = n*(n-1)/2 + m
      w  = g[mn]*y + h[mn]*x
      ;IF (abs(p2) LT 1d-38) THEN p2 = 0d ; RJW - Why have these lines?
      ;IF (abs(q)  LT 1d-38) THEN q  = 0d ; RJW - Why have these lines?
      IF scalar_input THEN BEGIN
        bbr += b[  n]*w*q
        bbt -= a[  n]*w*z
        IF bm THEN BEGIN
          IF not_bk THEN bi += a[n] * (g[mn]*x-h[mn]*y) * q  $
          ELSE           bi += a[n] * (g[mn]*x-h[mn]*y) * z
        ENDIF
      ENDIF ELSE BEGIN
        bbr += b[*,n]*w*q
        bbt -= a[*,n]*w*z
        IF bm THEN BEGIN
          qq = q
          ind = WHERE(bk, NULL = 1)
          IF (N_ELEMENTS(ind) NE 0) THEN qq[ind] = z[ind]
          bi += a[*,n] * (g[mn]*x-h[mn]*y) * qq
        ENDIF
      END
      xk = rec[mn] ; in IDL 8.4 it's faster to write this to xk, to use below twice.  In Matlab quicker to just use rec(nm)
      dp = cos_theta*z - sin_theta*q - xk*d2
      pm = cos_theta*q               - xk*p2
      d2 = z
      p2 = q
      z = dp
      q = pm
    ENDFOR
    d = sin_theta*d + cos_theta*p
    p = sin_theta*p
    IF bm THEN BEGIN
      bi  *= m_minus_1
      bbf += bi
    ENDIF
  ENDFOR

  ;br = bbr ; These don't change again
  ;bt = bbt ; These don't change again
  IF scalar_input THEN BEGIN
    IF not_bk THEN bf = bbf/sin_theta $
    ELSE BEGIN
      ;IF (cos_theta LT 0d) THEN bf = -bbf ELSE bf = bbf
      IF (cos_theta GE 0d) THEN bf = bbf ELSE bf = -bbf
    ENDELSE
  ENDIF ELSE BEGIN
    bf = bbf;% set size of array and do the 3rd case
    ind = WHERE(bk AND (cos_theta LT 0d), NULL = 1)
    IF (N_ELEMENTS(ind) NE 0) THEN bf[ind] = -bbf[ind]
    ind = WHERE(bk EQ 0, NULL = 1)
    IF (N_ELEMENTS(ind) NE 0) THEN bf[ind] =  bbf[ind]/sin_theta[ind]
  ENDELSE

  ;Brtp = [[br],[bt],[bf]]
  ;RETURN,Brtp
  RETURN,[[bbr],[bbt],[bf]]
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO _jovian_jrm09_internal_rjw_setup
  ; This is code to print values to the screen, which are then copied in to the above code.  We do not expect others to run this!
  ; This is mostly a copy/paste from earlier code of Marissa's - it's not optimized for speed - doesn't have to be in any way.
  ; We just want the output, which is independent of position.

  ON_ERROR, 2 ; % Exit code if an error in main, don't stop in code - no Matlab equivalent, just delete line in Matlab
  COMPILE_OPT HIDDEN


  ;============
  ;Begin hard-coding for JRM09
  ;Values from Connerney et al. 2020
  ;See supplemental online information Table S1, https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2018GL077312&file=grl57087-sup-0005-2018GL077312-ds01.txt
  ;============
  g = [0.d, 0.d, 410244.7d, -71498.3d, 11670.4d, -56835.8d, 48689.5d, $
    4018.6d, -37791.1d, 15926.3d, -2710.5d, $
    -34645.4d, -8247.6d, -2406.1d, -11083.8d, -17837.2d, $
    -18023.6d, 4683.9d, 16160.0d, -16402.0d, -2600.7d, -3660.7d, $
    -20819.6d, 9992.9d, 11791.8d, -12574.7d, 2669.7d, 1113.2d, 7584.9d, $
    598.4d, 4665.9d, -6495.7d, -2516.5d, -6448.5d, 1855.3d, -2892.9d, 2968.d, $
    10059.2d, 1934.4d, -6702.9d, 153.7d, -4124.2d, -867.2d, -3740.6d, -732.4d, -2433.2d, $
    9671.8d, -3046.2d, 260.9d, 2071.3d, 3329.6d, -2523.1d, 1787.1d, -1148.2d, 1276.5d, -1976.8d, $
    -2299.5d, 2009.7d, 2127.8d, 3498.3d, 2967.6d, 16.3d, 1806.5d, -46.5d, 2897.8d, 574.5d, 1298.9d, replicate(0.d, 200)] ; replicates go further than we need...

  h = [0.d, 0.d, 0.d, 21330.5d, 0.d, -42027.3d, 19353.2d, 0.d, -32957.3d, 42084.5d, -27544.2d, $
    0.d, 31994.5d, 27811.2d, -926.1d, 367.1d, 0.d, 45347.9d, -749.0d, 6268.5d, 10859.6d, 9608.4d,$
    0.d, 14533.1d, -10592.9d, 568.6d, 12871.7d, -4147.8d, 3604.4d, $
    0.d, -7626.3d, -10948.4d, 2633.3d, 5394.2d, -6050.8d, -1526.d, -5684.2d, $
    0.d, -2409.7d, -11614.6d, 9287.0d, -911.9d, 2754.5d, -2446.1d, 1207.3d, -2887.3d, $
    0.d, -8467.4d, -1383.8d, 5697.7d, -2056.3d, 3081.5d, -721.2d, 1352.5d, -210.1d, 1567.6d, $
    0.d, -4692.6d, 4445.8d, -2378.6d, -2204.3d, 164.1d, -1361.6d, -2031.5d, 1411.8d, -714.3d, 1676.5d, replicate(0.d, 200)] ; replicates go further than we need...
    
  degree = 20
  order = 10
  ;============
  ;End parts that are hard-coded for JRM09
  ;============
  
  a = dindgen(degree)*0.d
  b = dindgen(degree)*0.d
  rec = h*0.d

  for n=1, degree do begin
    n2 = 2*n-1;
    n2 = n2*(n2-2);
    for m=1, n do begin
      mn=n*(n-1)/2+m;
      rec(mn)=double((n-m)*(n+m-2.d))/n2;
    endfor
  endfor
  s = 1.d;
  for n=2, degree do begin
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

  ; Now Rob Prints this out to copy/paste in to earlier code

  HELP,rec,g,h
  ; Only print out for size of rec, although size g and h are longer
  FOR z=0,2 DO BEGIN
    CASE z OF
      0: PRINT,'rec = [$'
      1: PRINT,'g   = [$'
      2: PRINT,'h   = [$'
    ENDCASE
    CASE z OF
      0: x = rec
      1: x = g
      2: x = h
    ENDCASE
    TXT = ''
    i = -1L
    WHILE i LT N_ELEMENTS(rec)-2L DO BEGIN
      i += 1L
      if x[i] EQ 0d THEN BEGIN
        TXT += '                               0d, '
      endif else begin
        TXT += STRING(x[i],FORMAT='(F32.23)')+'d, '
      endelse
      IF (i+1) MOD 4L EQ 0 THEN BEGIN
        PRINT,TXT+ '$'
        TXT = ''
      ENDIF
    ENDWHILE
    TXT = STRMID(TXT,0,STRLEN(TXT)-2)+']'
    PRINT,TXT
  ENDFOR
  PRINT,'COPY/PASTE these in to source code'
  PRINT,'NOTE: For Matlab, remove the "d"s and also remove index 0, to only have 3 values on each first line'
END
