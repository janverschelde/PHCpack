with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;

package Standard_Lseries_Polynomials is

-- DESCRIPTION :
--   A Laurent series polynomial is a polynomial in several variables
--   where the coefficients are Laurent series, represented by
--   (1) a vector of leading exponents of the Laurent series coefficients,
--   (2) coefficient vectors of the Laurent series for each monomial;
--   (3) exponent vectors of the monomials in the polynomial.

  procedure Write ( plead : in Standard_Integer_Vectors.Vector;
                    pcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                    pmons : in Standard_Integer_VecVecs.VecVec;
                    s : in string := "p" );

  -- DESCRIPTION :
  --   Writes a Laurent series polynomial.

  -- ON ENTRY :
  --   plead    leading exponents for the Laurent series coeffficients;
  --   pcffs    coefficient vectors of the Laurent series for each monomial;
  --   pmons    exponent vectors of the monomials;
  --   s        string used as the name of the polynomial.

  procedure Make_Random_Polynomial
              ( dim,nbr,deg,pwr,low,upp : in integer32;
                lead : out Standard_Integer_Vectors.Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec;
                mons : out Standard_Integer_VecVecs.VecVec );

  -- DESCRIPTION :
  --   Makes a random polynomial with Laurent series as coefficients.

  -- ON ENTRY :
  --   dim      the dimension is the number of variables;
  --   nbr      number of monomials in the polynomial;
  --   deg      degree of the series;
  --   pwr      largest power for every variable;
  --   low      lower bound on leading exponents of the series;
  --   upp      upper bound on leading exponents of the series.

  -- ON RETURN :
  --   lead     an array of range 1..nbr with the leading exponents
  --            of the power series coefficients;
  --   cffs     coefficient vectors of the power series coefficients;
  --   mons     exponents of the monomials in the polynomial.

  procedure Eval ( deg,mlead : in integer32;
                   cff : in Standard_Complex_Vectors.Link_to_Vector;
                   mon : in Standard_Integer_Vectors.Link_to_Vector;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates a monomial at a Laurent series.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   mlead    leading exponent of the Laurent series coefficient
  --            of the monomial;
  --   cff      Laurent series coefficient of the monomial;
  --   mon      exponents of the monomial;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation.

  -- ON RETURN :
  --   ye       leading exponent of the result of the evaluation;
  --   yc       coefficient vector of the value of the monomial.

  procedure Eval ( deg : in integer32;
                   plead : in Standard_Integer_Vectors.Vector;
                   pcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   pmons : in Standard_Integer_VecVecs.VecVec;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a Laurent series.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   plead    leading exponents of the Laurent series coefficients;
  --   pcffs    coefficient vectors of the Laurent series coefficients;
  --   pmons    exponents of the monomials in the polynomial;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation.

  -- ON RETURN :
  --   ye       leading exponent of the result of the evaluation;
  --   yc       coefficient vector of the value of the polynomial.

  function Index_of_Degrees
             ( mons : Standard_Integer_VecVecs.VecVec;
               idx : integer32;
               degs : Standard_Integer_Vectors.Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the index where to place the degrees respective
  --   to the first exponents in the range 1..idx.
  --   This function is needed in case tdx /= 0.

  -- ON ENTRY :
  --   mons    exponents of monomials processed from 1 to idx-1;
  --   idx     current free index in mons;
  --   degs    exponents of a new monomial.

  -- ON RETURN :
  --   idx is returned if the exponents in degs do not occur in mons,
  --   otherwise returns the index in mons where degs exponents occur.

  procedure Make_Series_Polynomial
              ( p : in Poly; dim,nvr,tdx,deg : in integer32;
                lead : out Standard_Integer_Vectors.Link_to_Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec;
                mons : out Standard_Integer_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Given a Laurent polynomial p and the index for t,
  --   makes the data structures to represent p as polynomial
  --   with Laurent series coefficients in t.

  -- ON ENTRY :
  --   p       a Laurent polynomial with possibly negative exponents in t;
  --   dim     total number of variables in p, including t,
  --   nvr     number of variables without t;
  --   tdx     index of t as one of the dim variables in p,
  --           if tdx is zero, then dim must equal nvr,
  --           otherwise nvr = dim - 1.

  -- ON RETRUN :
  --   lead    leading exponents of the Laurent series coefficients;
  --   cffs    coefficients in the Laurent series for each monomials;
  --   mons    exponents of the monomials.

end Standard_Lseries_Polynomials;
