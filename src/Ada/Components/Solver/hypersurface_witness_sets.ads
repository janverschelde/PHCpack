with text_io;                         use text_io;
with Standard_Natural_Numbers;        use Standard_Natural_Numbers;
with Standard_Floating_Numbers;       use Standard_Floating_Numbers;
with Standard_Complex_Numbers;        use Standard_Complex_Numbers;
with Standard_Complex_Vectors;        use Standard_Complex_Vectors;

package Hypersurface_Witness_Sets is

-- DESCRIPTION :
--   This package offers generic procedures to compute a witness set for
--   a hypersurface defined by one polynomial in several variables.

  generic

    with function f ( x : Vector ) return Complex_Number;
    -- returns the value of the polynomial at the point x

  procedure Silent_Root_Finder0
              ( degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out double_float );

  generic

    with function f ( x : Vector ) return Complex_Number;
    -- returns the value of the polynomial at the point x

  procedure Silent_Root_Finder1
              ( degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out double_float );

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
                degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t : out Vector; nrm : out double_float );

  generic

    with function f ( x : Vector ) return Complex_Number;
    -- returns the value of the polynomial at the point x

  procedure Reporting_Root_Finder1
              ( file : in file_type;
                degree : in natural32; eps : in double_float;
                max_it : in natural32; fail : out boolean;
                b,v : in Vector; t,err,res : out Vector;
                nrm : out double_float );

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

end Hypersurface_Witness_Sets;
