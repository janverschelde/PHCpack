with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Polynomials;

package Test_Standard_Coeff_Convolutions is

-- DESCRIPTION :
--   Tests the evaluation of the gradient of a polynomial in many variables,
--   in a power series of some fixed degree,
--   also for coefficient convolutions in standard double precision.

  function Leading_Coefficients
             ( s : Standard_Complex_Series_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of leading coefficients of the vector s.

  function One_Coefficients
             ( nbr,deg : integer32 )
             return Standard_Complex_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns a vector of range 1..nbr with vectors of range 0..deg
  --   to represent the constant one as the coefficients.

  function Standard_Make_Polynomial
             ( dim,deg : integer32;
               idx : Standard_Integer_VecVecs.VecVec;
               xps : Standard_Integer_VecVecs.VecVec;
               cff : Standard_Complex_Series_Vectors.Vector;
               expone,cffone : boolean )
             return Standard_CSeries_Polynomials.Poly;

  -- DESCRIPTION :
  --   Wraps the construction of a polynomial with series coefficients,
  --   with flags for special cases if expone and cffone.

  -- ON ENTRY :
  --   dim     dimension of the exponent vectors, number of variables;
  --   deg     degree of the power series;
  --   idx     exponent indices;
  --   xps     exponent vectors;
  --   cff     coefficients of the monomials;
  --   expone  true if all exponents are equal to one;
  --   cffone  true if all coefficients are equal to one.

  procedure Standard_Test ( dim,deg,nbr,pwr : in integer32;
                            expone,cffone : in boolean );

  -- DESCRIPTION :
  --   Generates a sequence of random exponents and tests the
  --   evaluation and differentiation in standard double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.
  --   expone   true if all exponents are equal to one;
  --   cffone   true if all coefficients are equal to one.

  procedure Standard_System_Test ( dim,deg,nbr,pwr : in integer32 );

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in standard double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables.

  procedure Standard_Input_Test ( deg : in integer32 );

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system,
  --   makes convolution circuits, evaluates and differentiates
  --   at a vector of random series, in double precision.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the degree, the dimension, and
  --   the number of monomials.  Then runs the tests.

end Test_Standard_Coeff_Convolutions;
