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

  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Standard_Complex_Vectors.Vector;
               wrk : in out Standard_Complex_VecVecs.VecVec;
               ydx : out Standard_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out double_float;
               maxng,mindg,rcogd : out double_float );
  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in DoblDobl_Complex_Vectors.Vector;
               wrk : in out DoblDobl_Complex_VecVecs.VecVec;
               ydx : out DoblDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out double_double;
               maxng,mindg,rcogd : out double_double );
  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in QuadDobl_Complex_Vectors.Vector;
               wrk : in out QuadDobl_Complex_VecVecs.VecVec;
               ydx : out QuadDobl_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out quad_double;
               maxng,mindg,rcogd : out quad_double );
  procedure Gradient_with_Inverse_Condition
             ( f,b : in Standard_Natural_VecVecs.VecVec;
               c,x : in Multprec_Complex_Vectors.Vector;
               wrk : in out Multprec_Complex_VecVecs.VecVec;
               ydx : out Multprec_Complex_Vectors.Vector;
               fxnrc,fxdrc,fxrco : out Floating_Number;
               maxng,mindg,rcogd : out Floating_Number );

  -- DESCRIPTION :
  --   Computes the value of the polynomial and its gradient at x,
  --   with exponents in f, b, and coefficients in c.
  --   The condition of the evaluation is computed in (fxnrc, fxdrc, fxrco)
  --   and (maxng, mindg, rcogd) defines the condition of the gradient at x.
 
  -- REQUIRED :
  --   The range of the vector ydx must be 0..x'last with its 
  --   0-th component the function value and the i-th
  --   component the i-th derivative of the sum at x.
  --   Moreover: numrco'range = ydx'range.
  --   The y serves as work space and has been allocated,
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
  --   maxng   largest numerator of the condition of the gradient at x,
  --           taken over every component of the gradient as the absolute
  --           value of the evaluated partial derivative at x;
  --   mindg   smallest denonnarator of the condition of the gradient at x,
  --           taken over every component k as the absolute value of ydx(k);
  --   rcodg   the inverse condition number of the gradient is maxng/mindg. 

end VarbPrec_Gradient_Evaluations;
