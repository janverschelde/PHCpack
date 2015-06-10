with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;

package VarbPrec_Gradient_Evaluations is

-- DESCRIPTION :
--   The reverse mode of algorithmic differentiation is applied
--   to evaluate and differentiate polynomials, given in distributed form, 
--   with the common factors splitted from the products of the variables.
--   The evaluation and differentiation is performed jointly with the 
--   computation of the inverse condition number.
--   Four different levels of precision are supported:
--   double, double double, quad double, and arbitrary multiprecision.
--   Condition numbers that turn out as larger than the inverse of the 
--   working precision must be recomputed at a higher precision.

-- PART I : ordinary polynomials

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.VecVec;
               ydx : out Standard_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out double_float;
               gxnrc,gxdrc,gxrco : out double_float );
  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.VecVec;
               ydx : out DoblDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out double_double;
               gxnrc,gxdrc,gxrco : out double_double );
  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.VecVec;
               ydx : out QuadDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out quad_double;
               gxnrc,gxdrc,gxrco : out quad_double );
  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out Floating_Number;
               gxnrc,gxdrc,gxrco : out Floating_Number );

  -- DESCRIPTION :
  --   Computes the value of the polynomial and its gradient at x,
  --   with exponents in f, b, and coefficients in c.
  --   The condition of the evaluation is computed in (fxnrc, fxdrc, fxrco)
  --   and (gxnrc, gxdrc, gxrco) defines the condition of the gradient at x.
  --   The condition of the gradient is taken as the smallest inverse
  --   condition number over all its components.
 
  -- REQUIRED :
  --   The range of the vector ydx must be 0..x'last with its 
  --   0-th component the function value and the i-th
  --   component the i-th derivative of the sum at x.
  --   Moreover: numrco'range = ydx'range.
  --   The argument wrk serves as work space and has been allocated,
  --   in particular wrk'range = b'range.

  -- ON ENTRY :
  --   f       defines the common factors in each monomial,
  --           these are the monomials with higher powers that
  --           appear both in the function values and the derivatives;
  --   b       bit vectors that define products of variables,
  --           b(k)(i) = 1 if the i-th variable occurs in monomial k,
  --           b(k)(i) = 0 otherwise;
  --   c       c(k) is the coefficient of the k-th monomial;
  --   x       values for the variables: where to evaluate at;
  --   wrk     serves as workspace for all evaluated monomials.

  -- ON RETURN :
  --   wrk     used workspace, filled with values;
  --   ydx     ydx(0) is the value of the polynomial at x,
  --           ydx(k) is the k-th derivative of the polynomial at x;
  --   fxnrc   numerator of the inverse condition number of the evalution,
  --           this is the absolute value of the evaluation number ydx(0);
  --   fxdrc   denominator of the inverse condition number of the evalution,
  --           as the sum of the absolute values of the evaluated terms;
  --   fxrco   inverse condition number of the evaluation problem at x;
  --   gxnrc   numerator of the inverse condition number of the gradient;
  --   gxdrc   denominator of the inverse condition number of the gradient;
  --   gxrco   the inverse condition of the gradient is the smallest over
  --           all the components of the gradient, equals gxnrc/gxdrc.

-- PART II : ordinary polynomial systems

  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in Standard_Complex_VecVecs.VecVec;
               x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.Array_of_VecVecs;
               ydx : in Standard_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out double_float;
               gxnrc,gxdrc,gxrco : out double_float );
  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in DoblDobl_Complex_VecVecs.VecVec;
               x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.Array_of_VecVecs;
               ydx : in DoblDobl_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out double_double;
               gxnrc,gxdrc,gxrco : out double_double );
  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in QuadDobl_Complex_VecVecs.VecVec;
               x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.Array_of_VecVecs;
               ydx : in QuadDobl_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out quad_double;
               gxnrc,gxdrc,gxrco : out quad_double );
  procedure Jacobian_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.Array_of_VecVecs;
               c : in Multprec_Complex_VecVecs.VecVec;
               x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
               ydx : in Multprec_Complex_VecVecs.VecVec;
               fxnrc,fxdrc,fxrco : out Floating_Number;
               gxnrc,gxdrc,gxrco : out Floating_Number );

  -- DESCRIPTION :
  --   Given in f, b, and c the common factors, the bit vectors, and the
  --   coefficients of a polynomial system, the polynomial system is
  --   evaluated and differentiated at x, along with the computation
  --   of the condition number of the evaluation and differentiation problem.
  --   The result of the evaluation and differentiation is in ydx.
  --   The condition of the evaluation is computed in (fxnrc, fxdrc, fxrco)
  --   and (maxng, mindg, rcogd) defines the condition of the Jacobian at x.

  -- REQUIRED :
  --   f'range = b'range = c'range = wrk'range = ydx'range, where f'length
  --   equals the number of polynomials in the system.

  -- ON ENTRY :
  --   f       common factors for each polynomial in the system;
  --   b       bit vectors that define the products of variables;
  --   c       coefficients of the terms in each polynomial in the system;
  --   x       point where to evaluate and differentiate at;
  --   wrk     work space for each polynomial.

  -- ON RETURN :
  --   wrk     values used in the evaluation and differentiation;
  --   ydx     ydx(k) is a vector range 0..x'last, its first position
  --           at ydx(k)(0) contains the value of the k-th polynomial
  --           evaluated at x and the other components of ydx(k) contain
  --           the components of the evaluated gradient of the k-th polynomial
  --           so: ydx(k)(1..x'last) contains the gradient
  --           which is the k-th column of the Jacobian matrix;
  --   fxnrc   numerator of the inverse condition number of the evalution,
  --           this is the maximum of the absolute values of the evaluated
  --           polynomials in the vector ydx(0);
  --   fxdrc   denominator of the inverse condition number of the evalution,
  --           as the maximum of the sums of the absolute values of the 
  --           evaluated terms;
  --   fxrco   inverse condition number of the evaluation problem at x;
  --   gxnrc   numerator of the smallest inverse condition number,
  --           computed over all gradients;
  --   gxdrc   denominator of the smallest inverse condition number,
  --           computed over all gradients;
  --   rcodg   the inverse condition number is gxnrc/gxdrc.

end VarbPrec_Gradient_Evaluations;
