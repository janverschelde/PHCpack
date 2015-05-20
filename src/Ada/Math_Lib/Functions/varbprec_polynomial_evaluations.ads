with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Laur_Systems;

package VarbPrec_Polynomial_Evaluations is

-- DESCRIPTION :
--   Based on the condition number of polynomial evaluation,
--   we can adjust the working precision to obtained the result
--   accurate to a wanted number of decimal places.
--   The inverse condition number has in its numerator the absolute
--   value of the function evaluated at the point.  In root finding
--   problems, this numerator leads to an excessive condition number
--   and we are mostly interested in knowing the magnitude of the
--   residual (whether it is 10^(-4) or 10^(-16)) and not in the
--   number of decimal places that are accurate.
--   Therefore, procedures that compute the inverse condition numbers
--   returns as well the numerators and their denominators.
--   For root finding problems, one should better ignore the numerator
--   and use the denominator of the inverse condition number as
--   the condition number for the evaluation problem.

-- PART I : ordinary polynomials

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Polynomials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out double_float );
  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Polynomials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out double_double );
  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Polynomials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out quad_double );
  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Polynomials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out Floating_Number );

  -- DESCRIPTION :
  --   Computes the inverse of the condition number of evaluating f at z.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f).

  -- ON ENTRY :
  --   f       a polynomial in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   fz      the value of the polynomial f at z;
  --   absfz   the magnitude of the function value;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals absfz/denrco.

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

-- PART II : Laurent polynomials

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Laurentials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out double_float );
  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Laurentials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out double_double );
  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Laurentials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out quad_double );
  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Laurentials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Numbers.Complex_Number;
               absfz,denrco,rco : out Floating_Number );

  -- DESCRIPTION :
  --   Computes the inverse of the condition number of evaluating f at z.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f).

  -- ON ENTRY :
  --   f       a polynomial in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   fz      the value of the polynomial f at z;
  --   absfz   the magnitude of the function value fz;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals numrco/absfz.

  function Inverse_Condition_Number
             ( f : Standard_Complex_Laurentials.Poly;
               z : Standard_Complex_Vectors.Vector ) return double_float; 
  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Laurentials.Poly;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double; 
  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Laurentials.Poly;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double; 
  function Inverse_Condition_Number
             ( f : Multprec_Complex_Laurentials.Poly;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number; 

  -- DESCRIPTION :
  --   Returns the inverse of the condition number of evaluating f at z.
  --   This number is the sum of the magnitudes of the coefficients
  --   times the magnitudes of the evaluated monomials, divided by
  --   the magnitude of the function value |f(z)|.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f) and z(i) /= 0
  --   for the negative exponents of variables x(i) in f.

  -- ON ENTRY :
  --   f       a polynomial in several variables;
  --   z       values for the variables that appear in f.

-- PART III : ordinary polynomial systems

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float );
  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double );
  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double );
  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number );

  -- DESCRIPTION :
  --   Evaluates the polynomial system f at z with the computation
  --   of the inverse condition number.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f).

  -- ON ENTRY :
  --   f       a polynomial system in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   fz      the polynomial system f evaluated at z;
  --   absfz   the magnitude of the function value;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals absfz/denrco = Inverse_Condition_Number(f,z).

  function Inverse_Condition_Number
             ( f : Standard_Complex_Poly_Systems.Poly_Sys;
               z : Standard_Complex_Vectors.Vector ) return double_float; 
  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double;
  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double;
  function Inverse_Condition_Number
             ( f : Multprec_Complex_Poly_Systems.Poly_Sys;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the smallest inverse condition number of evaluating
  --   every polynomial f(i) at z.

-- PART IV : Laurent polynomial systems

  procedure Inverse_Condition_Number
             ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
               z : in Standard_Complex_Vectors.Vector;
               fz : out Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float );
  procedure Inverse_Condition_Number
             ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               fz : out DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double );
  procedure Inverse_Condition_Number
             ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               fz : out QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double );
  procedure Inverse_Condition_Number
             ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               fz : out Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number );

  -- DESCRIPTION :
  --   Evaluates the Laurent polynomial system f at z,
  --   with the computation of the inverse condition number.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f) and z(i) /= 0
  --   for those indices i for which the i-th variable occurs with
  --   a negative exponent.

  -- ON ENTRY :
  --   f       a Laurent polynomial system in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   absfz   the magnitude of the function value;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals absfz/denrco = Inverse_Condition_Number(f,z);
  --   fz      the Laurent polynomial system f evaluated at z.

  function Inverse_Condition_Number
             ( f : Standard_Complex_Laur_Systems.Laur_Sys;
               z : Standard_Complex_Vectors.Vector ) return double_float; 
  function Inverse_Condition_Number
             ( f : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               z : DoblDobl_Complex_Vectors.Vector ) return double_double;
  function Inverse_Condition_Number
             ( f : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               z : QuadDobl_Complex_Vectors.Vector ) return quad_double;
  function Inverse_Condition_Number
             ( f : Multprec_Complex_Laur_Systems.Laur_Sys;
               z : Multprec_Complex_Vectors.Vector ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the smallest inverse condition number of evaluating
  --   every polynomial f(i) at z.

-- PART V : the wrappers

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Polynomials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Numbers.Complex_Number );
  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Polynomials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Numbers.Complex_Number );
  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Polynomials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Numbers.Complex_Number );
  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Polynomials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Evaluates the polynomial f at z with the computation
  --   of the inverse condition number.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f).

  -- ON ENTRY :
  --   f       a polynomial in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   absfz   the magnitude of the function value;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals absfz/denrco = Inverse_Condition_Number(f,z);
  --   fz      the polynomial f evaluated at z.

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Laurentials.Poly;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Numbers.Complex_Number );
  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Laurentials.Poly;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Numbers.Complex_Number );
  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Laurentials.Poly;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Numbers.Complex_Number );
  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Laurentials.Poly;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Evaluates the polynomial f at z with the computation
  --   of the inverse condition number.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f) and z(i) /= 0
  --   for those i for which the variable x(i) appears with negative
  --   exponent in f.

  -- ON ENTRY :
  --   f       a Laurent polynomial in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   absfz   the magnitude of the function value;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals absfz/denrco = Inverse_Condition_Number(f,z);
  --   fz      the Laurent polynomial f evaluated at z.

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Vectors.Vector );
  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Vectors.Vector );
  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Vectors.Vector );
  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Poly_Systems.Poly_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Analogue versions for systems of polynomial equations.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f).

  -- ON ENTRY :
  --   f       a polynomial system in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   absfz   the magnitude of the function value;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals absfz/denrco = Inverse_Condition_Number(f,z);
  --   fz      the polynomial system f evaluated at z.

  procedure Evaluate_with_Inverse_Condition
             ( f : in Standard_Complex_Laur_Systems.Laur_Sys;
               z : in Standard_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_float;
               fz : out Standard_Complex_Vectors.Vector );
  procedure Evaluate_with_Inverse_Condition
             ( f : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
               z : in DoblDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out double_double;
               fz : out DoblDobl_Complex_Vectors.Vector );
  procedure Evaluate_with_Inverse_Condition
             ( f : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
               z : in QuadDobl_Complex_Vectors.Vector;
               absfz,denrco,rco : out quad_double;
               fz : out QuadDobl_Complex_Vectors.Vector );
  procedure Evaluate_with_Inverse_Condition
             ( f : in Multprec_Complex_Laur_Systems.Laur_Sys;
               z : in Multprec_Complex_Vectors.Vector;
               absfz,denrco,rco : out Floating_Number;
               fz : out Multprec_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Analogue versions for systems of Laurent polynomial equations.

  -- REQUIRED : z'range = 1..Number_of_Unknowns(f) and z(i) /= 0
  --   for those indices i for which the i-th variable occurs with
  --   a negative exponent.

  -- ON ENTRY :
  --   f       a Laurent polynomial system in several variables;
  --   z       values for the variables that appear in f.

  -- ON RETURN :
  --   absfz   the magnitude of the function value;
  --   denrco  denominator of the inverse condition number, equals the
  --           sum of the magnitudes of the coefficients times the
  --           magnitudes of the evaluated monomials;
  --   rco     equals absfz/denrco = Inverse_Condition_Number(f,z);
  --   fz      the Laurent polynomial system f evaluated at z.

end VarbPrec_Polynomial_Evaluations;
