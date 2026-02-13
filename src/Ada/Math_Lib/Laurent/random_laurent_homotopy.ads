with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;

package Random_Laurent_Homotopy is

-- DESCRIPTION :
--   Provides procedures to construct a homotopy of Laurent polynomials
--   where the coefficients contain leading terms of real power series,
--   with a prescribed solution series, for testing purposes.

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
                tpw : out Standard_Floating_VecVecs.VecVec;
                intpow : in boolean := false );

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
  --   nbm      nbm(i) equals the number of monomials in the i-th polynomial;
  --   intpow   if true, then converts the real powers of t to integers.

  -- ON RETURN :
  --   deg      supports of the monomials in the system;
  --   cff      coefficients of the monomials;
  --   tpw      powers of t in the coefficients of the homotopy.      

  procedure Random_Laurent_System
              ( nbp,dim,low,upp,size : in integer32;
                nbm : in Standard_Integer_Vectors.Vector;
                deg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                cff : out Standard_Complex_VecVecs.Array_of_VecVecs;
                tpw : out Standard_Floating_VecVecs.Array_of_VecVecs;
                intpow : in boolean := false );

  -- DESCRIPTION :
  --   Generates a random Laurent polynomial system,
  --   where the coefficients are real powered series of the given size.

  -- REQUIRED :
  --   nbm'range = deg'range = cff'range = ctp'range = 1..nbp.

  -- ON ENTRY :
  --   nbp      number of polynomials in the system;
  --   dim      number of variables in the polynomials;
  --   low      lower bound for the exponents in the monomials;
  --   upp      upper bound for the exponents in the monomials,
  --   nbm      nbm(i) equals the number of monomials in the i-th polynomial;
  --   size     size of the power series;
  --   intpow   if true, then converts the real powers of t to integers.

  -- ON RETURN :
  --   deg      supports of the monomials in the system;
  --   cff      coefficients of the monomials;
  --   tpw      powers of t in the coefficients of the homotopy.      

  procedure Product_Homotopy_Polynomial
              ( pdg : in Standard_Integer_VecVecs.VecVec;
                pcf : in Standard_Complex_Vectors.Vector;
                ptp : in Standard_Floating_Vectors.Vector;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec;
                idxfac : in integer32;
                hdg : out Standard_Integer_VecVecs.VecVec;
                hcf : out Standard_Complex_Vectors.Vector;
                htp : out Standard_Floating_Vectors.Vector;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given a Laurent polynomial and a power series,
  --   makes a homotopy polynomial that has the power series as a factor.

  -- REQUIRED :
  --   The size of the vectors hdg, hcf, and htp is the product
  --   of the number of monomials and two plus the number of terms 
  --   in the series.

  -- ON ENTRY :
  --   pdg      exponents of the monomials in the system;
  --   pcf      coefficients of the system;
  --   ptp      leading powers of the series coefficients of the system;
  --   scf      coefficients of the power series solution;
  --   spw      real powers of the series solution;
  --   idxfac   index i of the variable used for the factor x(i) - series;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   hdg      exponents of the monomials in the homotopy;
  --   hcf      coefficients of the homotopy;
  --   htp      leading powers of the series coefficients of the homotopy.

  procedure Product_Homotopy
              ( pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pcf : in Standard_Complex_VecVecs.VecVec;
                ptp : in Standard_Floating_VecVecs.VecVec;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                htp : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given a Laurent polynomial system and a power series,
  --   makes a homotopy that has this power series as a solution.
  --   The polynomials in the homotopy are the product of a factor
  --   of the form (x - series)*polynomial.

  -- REQUIRED :
  --   All vectors have 1..dim as the same range,
  --   where dim is the number of polynomials and the number of variables.

  -- ON ENTRY :
  --   pdg      exponents of the monomials in the system;
  --   pcf      coefficients of the system;
  --   ptp      leading powers of the series coefficients of the system;
  --   scf      coefficients of the power series solution;
  --   spw      real powers of the series solution;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   hdg      exponents of the monomials in the homotopy;
  --   hcf      coefficients of the homotopy;
  --   htp      leading powers of the series coefficients of the homotopy.

  procedure Canonical_Binomial
              ( idx,dim : in integer32;
                bdg : out Standard_Integer_VecVecs.VecVec;
                bcf : out Standard_Complex_Vectors.Vector;
                btp : out Standard_Floating_Vectors.Vector );

  -- DESCRIPTION :
  --    Defines a canonical binomial of the form x(idx) - 1,
  --    where idx is in the range 1..dim.

  -- REQUIRED : bdg'range = 1..2 = bcf'range = btp'range.

  -- ON ENTRY :
  --   idx      index of the current variable;
  --   dim      total number of variables.

  -- ON RETURN :
  --   bdg      degrees of the monomials;
  --   bcf      corresponding coefficients of the monomials;
  --   btp      powers of t are zero.

  procedure Binomial_Homotopy
              ( pdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                pcf : in Standard_Complex_VecVecs.VecVec;
                ptp : in Standard_Floating_VecVecs.VecVec;
                hdg : out Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : out Standard_Complex_VecVecs.VecVec;
                htp : out Standard_Floating_VecVecs.VecVec;
                vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Given a Laurent polynomial system and a power series,
  --   makes a homotopy that starts at canonical binomials.
  --   The polynomials in the homotopy are of the form
  --   of the form (x - 1) + polynomial.

  -- REQUIRED :
  --   All vectors have 1..dim as the same range,
  --   where dim is the number of polynomials and the number of variables.

  -- ON ENTRY :
  --   pdg      exponents of the monomials in the system;
  --   pcf      coefficients of the system;
  --   ptp      leading powers of the series coefficients of the system;
  --   vrblvl   is the verbose level.

  -- ON RETURN :
  --   hdg      exponents of the monomials in the homotopy;
  --   hcf      coefficients of the homotopy;
  --   htp      leading powers of the series coefficients of the homotopy.

  procedure Scale_Homotopy_Powers
              ( hct : in out Standard_Floating_Vectors.Vector );
  procedure Scale_Homotopy_Powers
              ( hct : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Subtracts the minimum power from all other powers in hct.

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

  procedure Test_Product_Homotopy
              ( hdg : in Standard_Integer_VecVecs.Array_of_VecVecs;
                hcf : in Standard_Complex_VecVecs.VecVec;
                htp : in Standard_Floating_VecVecs.VecVec;
                scf : in Standard_Complex_VecVecs.VecVec;
                spw : in Standard_Floating_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Evaluates the series with coefficients in scf and powers in spw
  --   at a random value for t into the homotopy, with supports in hdg,
  --   coefficients in hcf, and leading powers in htp.

end Random_Laurent_Homotopy;
