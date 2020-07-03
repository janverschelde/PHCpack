with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Coefficient_Convolutions;

package Shift_Coefficient_Convolutions is

-- DESCRIPTION :
--   Provides operations to shift the coefficients of convolution circuits,
--   represented by vectors with splitted real and imaginary parts.

-- COMPUTING POWERS OF THE VALUE IN THE SHIFT :

  procedure Powers_of_Shift
              ( pwt : in Standard_Floating_Vectors.Link_to_Vector;
                t : in double_float );

  -- DESCRIPTION :
  --   Computes all powers of t, needed to shift a coefficient
  --   vector of a series of degree deg, deg = pwt'last.

  -- REQUIRED : pwt'range = 0..deg, where deg is the degree
  --   of a power series, and deg > 0.

  -- ON ENTRY :
  --   pwt      space allocated for deg powers of t;
  --   t        value used in the shift.

  -- ON RETURN :
  --   pwt      powers of t.

  procedure Powers_of_Shift
              ( rpwt,ipwt : in Standard_Floating_Vectors.Link_to_Vector;
                rpt,ipt : in double_float );

  -- DESCRIPTION :
  --   Computes all powers of t = rpt + ipt*i,
  --   needed to shift a coefficient vector of a series of degree deg,
  --   deg = rpwt'last = ipwt'last.

  -- REQUIRED : rpwt'range = ipwt'range = 0..deg,
  --   where deg is the degree of a power series, and deg > 0.

  -- ON ENTRY :
  --   rpwt     space allocated for real parts of deg powers of t;
  --   ipwt     space allocated for imaginary parts of deg powers of t;
  --   rpt      real part of the value used in the shift;
  --   ipt      imaginary part of the value used in the shift.

  -- ON RETURN :
  --   rpwt     real parts of the powers of t;
  --   ipwt     imaginary parts of the powers of t.

-- SHIFTING COEFFICIENTS OF POWER SERIES :

  procedure Shift ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shift the real and imaginary parts in rcf and icf
  --   for the real value in pwt(1), using rwk and iwk as work vectors.

  -- REQUIRED : rcf'range = icf'range = rwk'range = iwk'range = 0..deg,
  --   where deg is the degree of the series.

  -- ON ENTRY :
  --   rcf      real parts of the complex coefficients of a series;
  --   icf      imaginary parts of the complex coefficients of a series;
  --   rwk      space allocated for the same range as rcf;
  --   iwk      space allocated for the same range as icf;
  --   pwt      powers of the values used in the shift.

  -- ON RETURN :
  --   rcf      real parts of the coefficients shifted with pwt(1);
  --   icf      imaginary parts of the coefficients shifted with pwt(1).

  procedure Map ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                  rsh,ish : in Standard_Floating_Vectors.Link_to_Vector;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shift the real and imaginary parts in rcf and icf
  --   for the real value in pwt(1),
  --   mapping the shifted coefficients in rsh and ish.

  -- REQUIRED : rcf'range = icf'range = rsh'range = ish'range = 0..deg,
  --   where deg is the degree of the series.

  -- ON ENTRY :
  --   rcf      real parts of the complex coefficients of a series;
  --   icf      imaginary parts of the complex coefficients of a series;
  --   rsh      space allocated for the same range as rcf;
  --   ish      space allocated for the same range as icf;
  --   pwt      powers of the values used in the shift.

  -- ON RETURN :
  --   rsh      real parts of the coefficients shifted with pwt(1);
  --   ish      imaginary parts of the coefficients shifted with pwt(1).

  procedure Shift ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shift the real and imaginary parts in rcf and icf
  --   for the real part rpwt(1) and imaginary part in ipwt(1),
  --   using rwk and iwk as work vectors.

  -- REQUIRED : rcf'range = icf'range = rwk'range = iwk'range = 0..deg,
  --   where deg is the degree of the series, and deg > 0.

  -- ON ENTRY :
  --   rcf      real parts of the complex coefficients of a series;
  --   icf      imaginary parts of the complex coefficients of a series;
  --   rwk      space allocated for the same range as rcf;
  --   iwk      space allocated for the same range as icf;
  --   rpwt     real part of the powers of the value used in the shift;
  --   ipwt     imaginary part of the powers of the value used in the shift.

  -- ON RETURN :
  --   rcf      real parts of the shifted coefficients;
  --   icf      imaginary parts of the shifted coefficients.

  procedure Map ( rcf,icf : in Standard_Floating_Vectors.Link_to_Vector;
                  rsh,ish : in Standard_Floating_Vectors.Link_to_Vector;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shift the real and imaginary parts in rcf and icf
  --   for the real part rpwt(1) and imaginary part in ipwt(1),
  --   mapping the shift coefficients in rsh and ish.

  -- REQUIRED : rcf'range = icf'range = rsh'range = ish'range = 0..deg,
  --   where deg is the degree of the series, and deg > 0.

  -- ON ENTRY :
  --   rcf      real parts of the complex coefficients of a series;
  --   icf      imaginary parts of the complex coefficients of a series;
  --   rsh      space allocated for the same range as rcf;
  --   ish      space allocated for the same range as icf;
  --   rpwt     real part of the powers of the value used in the shift;
  --   ipwt     imaginary part of the powers of the value used in the shift.

  -- ON RETURN :
  --   rsh      real parts of the shifted coefficients;
  --   ish      imaginary parts of the shifted coefficients.

  procedure Shift ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the shift procedure to all (rcf(k), icf(k)) pairs,
  --   for the value in pwt(1), using rwk and iwk as work vectors.

  -- REQUIRED : rcf'range = icf'range, and for all k in rcf'range:
  --   rcf(k)'range = icf(k)'range = rwk'range = iwk'range = 0..deg,
  --   where deg is the degree of the series, and deg > 0.

  -- ON ENTRY :
  --   rcf      vector of series coefficients with the
  --            real parts of the complex coefficients;
  --   icf      vector of series coefficients with the
  --            imaginary parts of the complex coefficients;
  --   rwk      space allocated for the same range as rcf;
  --   iwk      space allocated for the same range as icf;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   rcf      real parts of the shifted coefficients;
  --   icf      imaginary parts of the shifted coefficients.

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  rsh,ish : in Standard_Floating_VecVecs.VecVec;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the Map procedure to all (rcf(k), icf(k)) pairs,
  --   for the value in pwt(1), mapping the results into the
  --   (rsh(k), ish(k)) pairs of coefficients.

  -- REQUIRED : rcf'range = icf'range, and for all k in rcf'range:
  --   rcf(k)'range = icf(k)'range = rsh(k)'range = ish(k)'range = 0..deg,
  --   where deg is the degree of the series, and deg > 0.

  -- ON ENTRY :
  --   rcf      vector of series coefficients with the
  --            real parts of the complex coefficients;
  --   icf      vector of series coefficients with the
  --            imaginary parts of the complex coefficients;
  --   rsh      space allocated for the same range as rcf;
  --   ish      space allocated for the same range as icf;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   rsh      real parts of the shifted coefficients;
  --   ish      imaginary parts of the shifted coefficients.

  procedure Shift ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the shift procedure to all (rcf(k), icf(k)) pairs,
  --   for the real part rpwt(1) and imaginary part in ipwt(1),
  --   using rwk and iwk as work vectors.

  -- REQUIRED : rcf'range = icf'range, and for all k in rcf'range:
  --   rcf(k)'range = icf(k)'range = rwk'range = iwk'range = 0..deg,
  --   where deg is the degree of the series, and deg > 0.

  -- ON ENTRY :
  --   rcf      vector of series coefficients with the
  --            real parts of the complex coefficients;
  --   icf      vector of series coefficients with the
  --            imaginary parts of the complex coefficients;
  --   rwk      space allocated for the same range as rcf;
  --   iwk      space allocated for the same range as icf;
  --   rpwt     real part of the powers of the value used in the shift;
  --   ipwt     imaginary part of the powers of the value used in the shift.

  -- ON RETURN :
  --   rcf      real parts of the shifted coefficients;
  --   icf      imaginary parts of the shifted coefficients.

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  rsh,ish : in Standard_Floating_VecVecs.VecVec;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Applies the shift procedure to all (rcf(k), icf(k)) pairs,
  --   for the real part rpwt(1) and imaginary part in ipwt(1),
  --   mapping the shifted coefficients into (rsh(k), ish(k)) pairs.

  -- REQUIRED : rcf'range = icf'range, and for all k in rcf'range:
  --   rcf(k)'range = icf(k)'range = rsh(k)'range = ish(k)'range = 0..deg,
  --   where deg is the degree of the series, and deg > 0.

  -- ON ENTRY :
  --   rcf      vector of series coefficients with the
  --            real parts of the complex coefficients;
  --   icf      vector of series coefficients with the
  --            imaginary parts of the complex coefficients;
  --   rsh      space allocated for the same range as rcf;
  --   ish      space allocated for the same range as icf;
  --   rpwt     real part of the powers of the value used in the shift;
  --   ipwt     imaginary part of the powers of the value used in the shift.

  -- ON RETURN :
  --   rsh      real parts of the shifted coefficients;
  --   ish      imaginary parts of the shifted coefficients.

-- SHIFTING CIRCUITS WITH REAL VALUE :

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Shift ( c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuits;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuit(s) c with the value in pwt(1),
  --   using rwk and iwk as work space vectors.

  -- REQUIRED : rwk'range = iwk'range = pwt'range = 0..deg,
  --   where deg is the degree of the series in c.

  -- ON ENTRY :
  --   c        circuit(s) with splitted coefficient vectors;
  --   rwk      space allocated for the same range as c.rcf;
  --   iwk      space allocated for the same range as c.icf;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   c        circuit(s) with shifted coefficients.

-- SHIFTING CIRCUITS WITH COMPLEX VALUE :

  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Shift ( c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Shift ( c : in Standard_Coefficient_Convolutions.Circuits;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuit(s) c with the value 
  --   in rpwt(1) + rpwt(1)*i, using rwk and iwk as work space vectors.

  -- REQUIRED : rwk'range = iwk'range = pwt'range = 0..deg,
  --   where deg is the degree of the series in c.

  -- ON ENTRY :
  --   c        circuit(s) with splitted coefficient vectors;
  --   rwk      space allocated for the same range as c.rcf;
  --   iwk      space allocated for the same range as c.icf;
  --   rpwt     real parts of the powers of the value used in the shift;
  --   ipwt     imaginary parts of the powers of the value used in the shift.

  -- ON RETURN :
  --   c        circuit(s) with shifted coefficients.

-- MAPPING COEFFICIENTS SHIFTED WITH REAL VALUE INTO A CIRCUIT :

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Circuit;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuit c with the value in pwt(1),
  --   mapping the shifts of rcf, icf into the coefficients of c.

  -- REQUIRED : pwt'range = 0..deg,
  --   where deg is the degree of the series in c.

  -- ON ENTRY :
  --   rcf      real parts of the series coefficients of circuits,
  --            rcf(0) is the real part of the series in the constant;
  --   icf      real parts of the series coefficients of circuits,
  --            icf(0) is the real part of the series in the constant;
  --   c        circuit with splitted coefficient vectors;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   c        circuit with shifted coefficients.

-- MAPPING COEFFICIENTS SHIFTED WITH COMPLEX VALUE INTO A CIRCUIT :

  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Circuit;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf,icf : in Standard_Floating_VecVecs.VecVec;
                  c : in Standard_Coefficient_Convolutions.Link_to_Circuit;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuit c with the value 
  --   in rpwt(1) + ipwt(1)*i,
  --   mapping the shifts of rcf, icf into the coefficients of c.

  -- REQUIRED : pwt'range = 0..deg,
  --   where deg is the degree of the series in c.

  -- ON ENTRY :
  --   rcf      real parts of the series coefficients of circuits,
  --            rcf(0) is the real part of the series in the constant;
  --   icf      real parts of the series coefficients of circuits,
  --            icf(0) is the real part of the series in the constant;
  --   c        circuit with splitted coefficient vectors;
  --   rpwt     real parts of the powers of the value used in the shift;
  --   ipwt     imaginary parts of the powers of the shift value.

  -- ON RETURN :
  --   c        circuit with shifted coefficients.

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuits c with the value in pwt(1),
  --   mapping the shifts of rcf, icf into the coefficients of c.

  -- REQUIRED : pwt'range = 0..deg,
  --   where deg is the degree of the series in c.

  -- ON ENTRY :
  --   rcf      real parts of the series coefficients of circuits,
  --            rcf(k)(0) is the real part of the series in the constant;
  --   icf      real parts of the series coefficients of circuits,
  --            icf(k)(0) is the real part of the series in the constant;
  --   c        circuits with splitted coefficient vectors;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   c        circuits with shifted coefficients.

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  c : in Standard_Coefficient_Convolutions.Circuits;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuit(s) c with the value 
  --   in rpwt(1) + ipwt(1)*i,
  --   mapping the shifts of rcf, icf into the coefficients of c.

  -- REQUIRED : pwt'range = 0..deg,
  --   where deg is the degree of the series in c.

  -- ON ENTRY :
  --   rcf      real parts of the series coefficients of circuits,
  --            rcf(k)(0) is the real part of the series in the constant;
  --   icf      real parts of the series coefficients of circuits,
  --            icf(k)(0) is the real part of the series in the constant;
  --   c        circuits with splitted coefficient vectors;
  --   rpwt     real parts of the powers of the value used in the shift;
  --   ipwt     imaginary parts of the powers of the shift value.

  -- ON RETURN :
  --   c        circuits with shifted coefficients.

-- SHIFTING SYSTEMS OF CIRCUITS :

  procedure Shift ( s : in Standard_Coefficient_Convolutions.System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Shift ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuits of s with the value in pwt(1),
  --   using rwk and iwk as work space vectors.

  -- REQUIRED : rwk'range = iwk'range = pwt'range = 0..deg,
  --   where deg is the degree of the series in s.crc.

  -- ON ENTRY :
  --   s        system of circuits with splitted coefficient vectors;
  --   rwk      space allocated for the same range as c.rcf;
  --   iwk      space allocated for the same range as c.icf;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   s        circuit with shifted coefficients.

  procedure Shift ( s : in Standard_Coefficient_Convolutions.System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Shift ( s : in Standard_Coefficient_Convolutions.Link_to_System;
                    rwk,iwk : in Standard_Floating_Vectors.Link_to_Vector;
                    rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                    ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Shifts all coefficients in the circuits of s with the value 
  --   in rpwt(1) + ipwt(1)*i, using rwk and iwk as work space vectors.

  -- REQUIRED : rwk'range = iwk'range = pwt'range = 0..deg,
  --   where deg is the degree of the series in s.crc.

  -- ON ENTRY :
  --   s        system of circuits with splitted coefficient vectors;
  --   rwk      space allocated for the same range as c.rcf;
  --   iwk      space allocated for the same range as c.icf;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   s        circuit with shifted coefficients.

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.System;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  pwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Maps the coefficients in rcf and icf shifted with the value
  --   in pwt(1) into the coefficients of s.

  -- REQUIRED : pwt'range = 0..deg,
  --   where deg is the degree of the series in s.crc.

  -- ON ENTRY :
  --   rcf      real parts of the series coefficients of circuits,
  --            rcf(k)(0) is the real part of the series in the constant;
  --   icf      real parts of the series coefficients of circuits,
  --            icf(k)(0) is the real part of the series in the constant;
  --   s        system of circuits with splitted coefficient vectors;
  --   pwt      the powers of the value used in the shift.

  -- ON RETURN :
  --   s        circuit with shifted coefficients.

  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.System;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf : in Standard_Floating_VecVecVecs.VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Map ( rcf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  icf : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                  s : in Standard_Coefficient_Convolutions.Link_to_System;
                  rpwt : in Standard_Floating_Vectors.Link_to_Vector;
                  ipwt : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Maps the coefficients in rcf and icf shifted with the value
  --   in rpwt(1) + ipwt(1)*i into the coefficients of s.

  -- REQUIRED : pwt'range = 0..deg,
  --   where deg is the degree of the series in s.crc.

  -- ON ENTRY :
  --   rcf      real parts of the series coefficients of circuits,
  --            rcf(k)(0) is the real part of the series in the constant;
  --   icf      real parts of the series coefficients of circuits,
  --            icf(k)(0) is the real part of the series in the constant;
  --   s        system of circuits with splitted coefficient vectors;
  --   rpwt     real parts of powers of the value used in the shift;
  --   ipwt     imaginary part of powers of the value used in the shift.

  -- ON RETURN :
  --   s        circuit with shifted coefficients.

end Shift_Coefficient_Convolutions;
