with text_io;                             use text_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;

package Certify_Factor_with_Trace is

-- DESCRIPTION :
--   To check whether we have sufficiently many generic points on a line
--   to represent a factor of a polynomial, we compute the linear trace,
--   in standard double, double double, or quad double precision,
--   without output or with output.

  procedure Certify_Factor
               ( p : in Standard_Complex_Polynomials.Poly;
                 b,v,w : in Standard_Complex_Vectors.Vector;
                 testres : out double_float );
  procedure Certify_Factor
               ( p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,w : in DoblDobl_Complex_Vectors.Vector;
                 testres : out double_float );
  procedure Certify_Factor
               ( p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,w : in QuadDobl_Complex_Vectors.Vector;
                 testres : out double_float );

  -- DESCRIPTION :
  --   Applies linear traces to certify one factor, without output.

  -- ON ENTRY :
  --   p         polynomial in n variables, contains a factor;
  --   b         offset vector for a random line b + t*v;
  --   v         direction of a random line b + t*v;
  --   w         witness points on a factor of p, on a random line.

  -- ON RETURN :
  --   testres   test residual: absolute value of the difference 
  --             between the value of the linear trace at the test points,
  --             and the sum of one component over all test points.

  procedure Certify_Factor
               ( file : in file_type;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,w : in Standard_Complex_Vectors.Vector;
                 testres : out double_float );
  procedure Certify_Factor
               ( file : in file_type;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,w : in DoblDobl_Complex_Vectors.Vector;
                 testres : out double_float );
  procedure Certify_Factor
               ( file : in file_type;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,w : in QuadDobl_Complex_Vectors.Vector;
                 testres : out double_float );

  -- DESCRIPTION :
  --   Applies linear traces to certify one factor,
  --   with intermediate output written to the file.

end Certify_Factor_with_Trace;
