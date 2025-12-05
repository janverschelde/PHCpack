with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Test_Leading_Evaluations is

-- DESCRIPTION :
--   Tests the evaluation of Laurent monomials at leading terms
--   of series with real powers.

  function Random_Monomial_Support
             ( dim,low,upp : integer32 )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of exponents of dimension dim,
  --   with values randomly generated between low and upp,
  --   as the support of a random monomial.

  function Random_Polynomial_Support
             ( nbr,dim,low,upp : integer32 )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns as many as nbr random exponents in a vector of range 1..nbr,
  --   of dimension dim, with values between low and upp,
  --   as the support of a random Laurent polynomial.

  procedure Random_Laurent_System
              ( nbp,dim,low,upp : in integer32;
                nbm : in Standard_Integer_Vectors.Vector;
                deg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : out Standard_Complex_VecVecs.VecVec;
                tpw : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Generates a random Laurent polynomial system,
  --   where the coefficients are leading powers of series.

  -- REQUIRED :
  --   nbm'range = deg'range = cff'range = ctp'range = 1..nbp.

  -- ON ENTRY :
  --   nbp      number of polynomials in the system;
  --   dim      number of variables in the polynomials;
  --   low      lower bound for the exponents in the monomials;
  --   upp      upper bound for the exponents in the monomials,
  --   nbm      nbm(i) equals the number of monomials in the i-th polynomial.

  -- ON RETURN :
  --   deg      supports of the monomials in the system;
  --   cff      coefficients of the monomials;
  --   tpw      powers of t in the coefficients of the homotopy.      

  function Random_Leading_Powers
             ( dim : integer32 ) return Standard_Floating_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..dim with random powers of series
  --   with real positive powers. 

  procedure Random_Power_Series
              ( dim : in integer32;
                nbt : in Standard_Integer_Vectors.Vector;
                cff : out Standard_Complex_VecVecs.VecVec;
                pwr : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Generates dim power series with random coefficients and real powers.

  -- REQUIRED :
  --   cff'range = nbt'range = pwr'range = 1..dim.

  -- ON ENTRY :
  --   dim      number of power series;
  --   nbt      nbt(i) equals the number of terms in the i-th series;
  --   cff      coefficients of the power series;
  --   pwr      powers of the series.

  procedure Random_Homotopy
              ( pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pcf : in Standard_Complex_VecVecs.VecVec;
                ptp : in Standard_Floating_VecVecs.VecVec;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                htp : out Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Given a Laurent polynomial system and a power series,
  --   makes a homotopy that has this power series as a solution.

  -- REQUIRED :
  --   All vectors have 1..dim as the same range,
  --   where dim is the number of polynomials and the number of variables.

  -- ON ENTRY :
  --   pdg      exponents of the monomials in the system;
  --   pcf      coefficients of the system;
  --   ptp      leading powers of the series coefficients of the system;
  --   scf      coefficients of the power series solution;
  --   spw      real powers of the series solution.

  -- ON RETURN :
  --   hdg      exponents of the monomials in the homotopy;
  --   hcf      coefficients of the homotopy;
  --   htp      leading powers of the series coefficients of the homotopy.

  function Evaluate_Series
             ( cff : Standard_Complex_VecVecs.VecVec;
               pwr : Standard_Floating_VecVecs.VecVec; tpt : double_float )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Given in cff and pwr are the coefficients and real powers of a series,
  --   and in tpt a value for t.  Returns the value of the series.

  function Evaluate_Homotopy
             ( deg : Standard_Integer_VecVecs.Array_of_VecVecs;
               cff : Standard_Complex_VecVecs.VecVec;
               tpw : Standard_Floating_VecVecs.VecVec;
               zpt : Standard_Complex_Vectors.Vector; tpt : double_float )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the evaluation of the homotopy defined by supports
  --   in deg, coefficients in cff, and leading powers in tpw,
  --   at the point zpt and parameter tpt.

  procedure Test_Random_Homotopy
              ( hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : in Standard_Complex_VecVecs.VecVec;
                htp : in Standard_Floating_VecVecs.VecVec;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates the series with coefficients in scf and powers in spw
  --   at a random value for t into the homotopy, with supports in hdg,
  --   coefficients in hcf, and leading powers in htp.

  procedure Test_Monomial_Derivative
              ( deg : in Standard_Integer_Vectors.Vector;
                pwr : in Standard_Floating_Vectors.Vector;
                cff : in Standard_Complex_Vectors.Vector;
                idx : in integer32; err : out double_float );

  -- DESCRIPTION :
  --   Tests the derivative of the monomial with exponents in deg
  --   at the series with leading powers in pwr and coefficients in cff,
  --   with respect to the variable index idx.
  --   On return in err is the magnitude of the error of a random point test.

  procedure Test_Monomial ( dim : in integer32 );
               
  -- DESCRIPTION :
  --   Tests monomial evaluation and differentation in dim many variables
  --   at the leading terms of a series with real positive powers.

  procedure Test_Polynomial ( nbr,dim : in integer32 );

  -- DESCRIPTION :
  --   Tests polynomial evaluation and differentiation in dim many variables,
  --   of a polynomial with nbr many terms, at the leading terms of a series
  --   with real positive powers.

  procedure Test_System
              ( nbp,dim : in integer32;
                nbm : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests polynomial evaluation in dim many variables,
  --   of a system of nbp polynomials, where the i-th polynomial
  --   has nbm(i) many terms, at the leading terms of a series with 
  --   real positive powers.

  procedure Test_Homotopy
              ( dim : in integer32;
                nbm,nbt : in Standard_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Tests generation of a Laurent homotopy of dim polynomials
  --   in dim variables, with number of monomials in nbm,
  --   with a generated power series solution, where the number
  --   of terms in each series is defined by nbt.

  procedure Main;

  -- DESCRIPTION :
  --   Prompts for the dimension and then launches a test.

end Test_Leading_Evaluations;
