with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Polynomials;
with Multprec_Complex_Vectors;
with Multprec_Complex_Polynomials;

package VarbPrec_Polynomial_Evaluations is

-- DESCRIPTION :
--   Based on the condition number of polynomial evaluation,
--   we can adjust the working precision to obtained the result
--   accurate to a wanted number of decimal places.

  function Inverse_Condition_Number
             ( f : Standard_Complex_Polynomials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float; 
  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Polynomials.Poly;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double; 
  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Polynomials.Poly;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double; 
  function Inverse_Condition_Number
             ( f : Multprec_Complex_Polynomials.Poly;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the inverse of the condition number of evaluating f at z.
  --   This number is the sum of the magnitudes of the coefficients
  --   times the magnitudes of the evaluated monomials, divided by
  --   the magnitude of the function value |f(z)|.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f).

  -- ON ENTRY :
  --   f       a polynomial in several variables;
  --   z       values for the variables that appear in f.

end VarbPrec_Polynomial_Evaluations;
