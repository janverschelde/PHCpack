with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;

package Triple_Double_Constants is

-- DESCRIPTION :
--   This package collects common definitions of triple double constants
--   for use in triple double mathematical functions.

  td_eps : constant double_float := 5.473822126268817e-48; -- 2^(-52*3-1)

-- PI and multiples and factors :

  pi : constant triple_double
     := create( 3.141592653589793116e+00, 1.224646799147353207e-16,
               -2.994769809718339666e-33);
  twopi : constant triple_double -- 2*pi
        := create( 6.283185307179586232e+00, 2.449293598294706414e-16,
                  -5.989539619436679332e-33);
  pi2 : constant triple_double -- pi/2
      := create( 1.570796326794896558e+00, 6.123233995736766036e-17,
                -1.497384904859169833e-33);
  pi4 : constant triple_double -- pi/4
      := create( 7.853981633974482790e-01, 3.061616997868383018e-17,
                -7.486924524295849165e-34);
  threepi4 : constant triple_double -- 3*pi/4
           := create( 2.356194490192344837e+00, 9.1848509936051484375e-17,
                      3.9168984647504003225e-33);
  pi1024 : constant triple_double -- pi/1204
         := create( 3.067961575771282340e-03, 1.195944139792337116e-19,
                   -2.924579892303066080e-36);

end Triple_Double_Constants;
