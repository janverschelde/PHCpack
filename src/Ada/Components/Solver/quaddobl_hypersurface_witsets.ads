with text_io;                         use text_io;
with Standard_Natural_Numbers;        use Standard_Natural_Numbers;
with Standard_Integer_Numbers;        use Standard_Integer_Numbers;
with Quad_Double_Numbers;             use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;        use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;        use QuadDobl_Complex_Vectors;

package QuadDobl_Hypersurface_Witsets is

-- DESCRIPTION :
--   This package offers generic procedures to compute a witness set for
--   a hypersurface defined by one polynomial in several variables,
--   in quad double precision.

-- PART I : the primitives in the Durand-Kerner method (aka Weierstrass)

  procedure Divided_Differences ( x : in Vector; f : in out Vector );

  -- DESCRIPTION :
  --   Computes in f the divided differences using the points in x.

  function Roots_of_Unity ( d : natural32 ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector with the d complex roots of unity.
  --   These serve as start values for the univariate root finder.

  function Compute_q ( i : integer32; a : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Computes the quotient needed in the Durand-Kerner step.

  procedure Write ( file : in file_type; z,err,res : in Vector );

  -- DESCRIPTION :
  --   Writes the current approximations for the roots in z,
  --   along with their residuals in res on the file.
  --   This procedure is called in the reporting version
  --   of the root finder.

-- PART II : generic procedures for the root finding method

  generic

    with function f ( x : Vector ) return Complex_Number;
    -- returns the value of the polynomial at the point x

  procedure Silent_Root_Finder0
              ( degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out quad_double );

  generic

    with function f ( x : Vector ) return Complex_Number;
    -- returns the value of the polynomial at the point x

  procedure Silent_Root_Finder1
              ( degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out quad_double );

  -- DESCRIPTION :
  --   Returns as many t's as the degree: f(x+t*v) = 0, modulo eps.

  -- REQUIRED : t'range = 1..degree.
  
  -- ON ENTRY :
  --   degree   degree of the hypersurface;
  --   eps      requirement on the accuracy of the roots;
  --   max_it   maximal number of iterations allowed;
  --   b        offset vector for the affine line b + t*v;
  --   v        direction for the affine line b + t*v.

  -- ON RETURN :
  --   fail     false if the accuracy eps was reached within max_it,
  --            true otherwise;
  --   t        degree many points satisfying f(x+t*v) = 0 if not fail;
  --   err      last correction term to the components of t;
  --   res      residual vector for all components of t;
  --   nrm      maximum norm of the residual vector.

  generic

    with function f ( x : Vector ) return Complex_Number;
    -- returns the value of the polynomial at the point x

  procedure Reporting_Root_Finder0
              ( file : in file_type;
                degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out quad_double );

  generic

    with function f ( x : Vector ) return Complex_Number;
    -- returns the value of the polynomial at the point x

  procedure Reporting_Root_Finder1
              ( file : in file_type;
                degree : in natural32; eps : in quad_double;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out quad_double );

  -- DESCRIPTION :
  --   Returns as many t's as the degree: f(x+t*v) = 0, modulo eps.

  -- REQUIRED : t'range = 1..degree.
  
  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   degree   degree of the hypersurface;
  --   eps      requirement on the accuracy of the roots;
  --   max_it   maximal number of iterations allowed;
  --   b        offset vector for the affine line b + t*v;
  --   v        direction for the affine line b + t*v.

  -- ON RETURN :
  --   fail     false if the accuracy eps was reached within max_it,
  --            true otherwise;
  --   t        degree many points satisfying f(x+t*v) = 0 if not fail;
  --   err      last correction term to the components of t;
  --   res      residual vector for all components of t;
  --   nrm      maximum norm of the residual vector res.

end QuadDobl_Hypersurface_Witsets;
