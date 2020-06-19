with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;

package Shift_Coefficient_Convolutions is

-- DESCRIPTION :
--   Provides operations to shift the coefficients of convolution circuits,
--   represented by vectors with splitted real and imaginary parts.

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

end Shift_Coefficient_Convolutions;
