with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;

package Double_Double_Constants is

-- DESCRIPTION :
--   This package collects common definitions of double double constants
--   for use in double double mathematical functions.

  dd_eps : constant double_float := 4.93038065763132e-32;     -- 2^-104

-- PI and some multiples and fractions :

  pi : constant double_double
     := Create(3.141592653589793116e+00,1.224646799147353207e-16);
  twopi : constant double_double -- 2*pi
        := create(6.283185307179586232e+00,2.449293598294706414e-16);
  pi2 : constant double_double   -- pi/2
      := create(1.570796326794896558e+00,6.123233995736766036e-17);
  pi4 : constant double_double   -- pi/4
      := create(7.853981633974482790e-01,3.061616997868383018e-17);
  threepi4 : constant double_double -- 3*pi/4
           := create(2.356194490192344837e+00,9.1848509936051484375e-17);
  pi16 : constant double_double
       := create(1.963495408493620697e-01,7.654042494670957545e-18);

-- INVERSE FACTORIALS FOR TAYLOR EXPANSION :

  i_fac0  : constant double_double
          := Create(1.66666666666666657E-01, 9.25185853854297066E-18);
  i_fac1  : constant double_double
          := Create(4.16666666666666644E-02, 2.31296463463574266E-18);
  i_fac2  : constant double_double
          := Create(8.33333333333333322E-03, 1.15648231731787138E-19);
  i_fac3  : constant double_double
          := Create(1.38888888888888894E-03,-5.30054395437357706E-20);
  i_fac4  : constant double_double
          := Create(1.98412698412698413E-04, 1.72095582934207053E-22);
  i_fac5  : constant double_double
          := Create(2.48015873015873016E-05, 2.15119478667758816E-23);
  i_fac6  : constant double_double
          := Create(2.75573192239858925E-06,-1.85839327404647208E-22);
  i_fac7  : constant double_double
          := Create(2.75573192239858883E-07, 2.37677146222502973E-23);
  i_fac8  : constant double_double
          := Create(2.50521083854417202E-08,-1.44881407093591197E-24);
  i_fac9  : constant double_double
          := Create(2.08767569878681002E-09,-1.20734505911325997E-25);
  i_fac10 : constant double_double
          := Create(1.60590438368216133E-10, 1.25852945887520981E-26);
  i_fac11 : constant double_double
          := Create(1.14707455977297245E-11, 2.06555127528307454E-28);
  i_fac12 : constant double_double
          := Create(7.64716373181981641E-13, 7.03872877733453001E-30);
  i_fac13 : constant double_double
          := Create(4.77947733238738525E-14, 4.39920548583408126E-31);
  i_fac14 : constant double_double
          := Create(2.81145725434552060E-15, 1.65088427308614326E-31);

  n_inv_fact : constant natural := 15;
  i_fac : array(0..n_inv_fact-1) of double_double
        := (i_fac0,i_fac1,i_fac2,i_fac3,i_fac4,i_fac5,i_fac6,i_fac7,
            i_fac8,i_fac9,i_fac10,i_fac11,i_fac12,i_fac13,i_fac14);

-- TABLES of sin(k*pi/16) and cos(k*pi/16)

  sin_t0 : constant double_double
         := create(1.950903220161282758e-01,-7.991079068461731263e-18);
  sin_t1 : constant double_double
         := create(3.826834323650897818e-01,-1.005077269646158761e-17);
  sin_t2 : constant double_double
         := create(5.555702330196021776e-01, 4.709410940561676821e-17);
  sin_t3 : constant double_double
         := create(7.071067811865475727e-01,-4.833646656726456726e-17);

  cos_t0 : constant double_double
         := create(9.807852804032304306e-01, 1.854693999782500573e-17);
  cos_t1 : constant double_double
         := create(9.238795325112867385e-01, 1.764504708433667706e-17);
  cos_t2 : constant double_double
         := create(8.314696123025452357e-01, 1.407385698472802389e-18);
  cos_t3 : constant double_double
         := create(7.071067811865475727e-01,-4.833646656726456726e-17);

  sin_table : array(0..3) of double_double := (sin_t0,sin_t1,sin_t2,sin_t3);
  cos_table : array(0..3) of double_double := (cos_t0,cos_t1,cos_t2,cos_t3);

end Double_Double_Constants;
