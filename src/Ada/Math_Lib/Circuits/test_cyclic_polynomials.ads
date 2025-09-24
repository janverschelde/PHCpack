with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;

package Test_Cyclic_Polynomials is

-- DESCRIPTION :
--   The goal of this test program is to test the performance of the
--   evaluation and differentiation schemes in reverse mode on the
--   benchmark problem of cyclic n-roots.
--   Run the performance test on n = 10 and m = 50000 to see the improvement.

  function Indexed_Support
             ( s : Standard_Natural_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the indexed supports of s.

  function Indexed_Supports
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return Standard_Integer_VecVecs.Array_of_VecVecs;

  -- DESCRIPTION :
  --   Returns the indexed supports of s.

  function Gradient_of_Cyclic
             ( s : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector ) 
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its support in s at x, in standard double precision.

  function Gradient_of_Cyclic
             ( s : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector ) 
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its support in s at x, in double double precision.

  function Gradient_of_Cyclic
             ( s : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector ) 
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its support in s at x, in quad double precision.

  function Indexed_Gradient_of_Cyclic
             ( s : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector ) 
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its indexed support in s at x, in double precision.

  function Indexed_Gradient_of_Cyclic
             ( s : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector ) 
             return DoblDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its indexed support in s at x, in double double precision.

  function Indexed_Gradient_of_Cyclic
             ( s : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector ) 
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its indexed support in s at x, in quad double precision.

  procedure Evaluate_and_Differentiate
              ( s : in Standard_Natural_VecVecs.Array_of_VecVecs;
                x : in Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in standard double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Evaluate_and_Differentiate
              ( s : in Standard_Natural_VecVecs.Array_of_VecVecs;
                x : in DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in double double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Evaluate_and_Differentiate
              ( s : in Standard_Natural_VecVecs.Array_of_VecVecs;
                x : in QuadDobl_Complex_Vectors.Vector;
                y : out QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in quad double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Indexed_Evaluate_and_Differentiate
              ( s : in Standard_Integer_VecVecs.Array_of_VecVecs;
                x : in Standard_Complex_Vectors.Vector;
                z : in out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   using the indexed supports in s, in standard double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  -- REQUIRED :
  --   The vector z is used as work space and is of range 0..x'last.

  procedure Indexed_Evaluate_and_Differentiate
              ( s : in Standard_Integer_VecVecs.Array_of_VecVecs;
                x : in DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   using the indexed supports in s, in double double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Indexed_Evaluate_and_Differentiate
              ( s : in Standard_Integer_VecVecs.Array_of_VecVecs;
                x : in QuadDobl_Complex_Vectors.Vector;
                y : out QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   using the indexed supports in s, in quad double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Evaluate_and_Differentiate
              ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in standard double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Evaluate_and_Differentiate
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in double double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Evaluate_and_Differentiate
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in QuadDobl_Complex_Vectors.Vector;
                y : out QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix );

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in quad double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  procedure Standard_Random_Point_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation of the polynomials
  --   in the cyclic n-roots problem at a random point,
  --   in double precision.

  procedure DoblDobl_Random_Point_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation of the polynomials
  --   in the cyclic n-roots problem at a random point,
  --   in double double precision.

  procedure QuadDobl_Random_Point_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation of the polynomials
  --   in the cyclic n-roots problem at a random point,
  --   in quad double precision.

  procedure Standard_Performance_Test ( n : integer32 );

  -- DESCRIPTION :
  --   Does a performance test on the cyclic n-roots problem,
  --   in double precision.

  procedure DoblDobl_Performance_Test ( n : integer32 );

  -- DESCRIPTION :
  --   Does a performance test on the cyclic n-roots problem,
  --   in double double precision.

  procedure QuadDobl_Performance_Test ( n : integer32 );

  -- DESCRIPTION :
  --   Does a performance test on the cyclic n-roots problem,
  --   in quad double precision.

  procedure Double_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Generates the support for the cyclic n-roots system and
  --   the coefficients in double precision.

  procedure DoblDobl_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Generates the support for the cyclic n-roots system and
  --   the coefficients in double double precision.

  procedure QuadDobl_Test ( n : in integer32 );

  -- DESCRIPTION :
  --   Generates the support for the cyclic n-roots system and
  --   the coefficients in quad double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts the user for the dimension and the precision.
  --   Runs then the selected test.

end test_cyclic_polynomials;
