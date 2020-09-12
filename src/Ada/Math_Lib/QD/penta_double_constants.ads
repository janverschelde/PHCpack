with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;

package Penta_Double_Constants is

-- DESCRIPTION :
--   This package collects common definitions of penta double constants
--   for use in penta double mathematical functions.

  pd_eps : constant double_float := 6.747006683667535e-80; -- 2^(-263)

-- PI and multiples and factors :

  pi : constant penta_double
     := create( 3.14159265358979312E+00, 1.22464679914735321E-16,
               -2.99476980971833967E-33, 1.11245422086336528E-49,
                5.67223197964031574E-66);

  twopi : constant penta_double
        := create( 6.28318530717958623E+00, 2.44929359829470641E-16,
                  -5.98953961943667933E-33, 2.22490844172673056E-49,
                   1.13444639592806315E-65);

  pi2 : constant penta_double
      := create( 1.57079632679489656E+00, 6.12323399573676604E-17,
                -1.49738490485916983E-33, 5.56227110431682641E-50,
                 2.83611598982015787E-66);

  pi4 : constant penta_double
      := create( 7.85398163397448279E-01, 3.06161699786838302E-17,
                -7.48692452429584916E-34, 2.78113555215841320E-50,
                 1.41805799491007894E-66);

  threepi4 : constant penta_double
           := create( 2.35619449019234484E+00, 9.18485099360514844E-17,
                      3.91689846475040032E-33, -2.58679816327048642E-49,
                     -4.93609888149662566E-67);

  pi1024 : constant penta_double
         := create( 3.06796157577128234E-03, 1.19594413979233712E-19,
                   -2.92457989230306608E-36, 1.08638107506188016E-52,
                    5.53928904261749584E-69);

end Penta_Double_Constants;
