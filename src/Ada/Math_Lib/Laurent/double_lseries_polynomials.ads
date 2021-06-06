with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Integer_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecVecs;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_JacoMats;    use Standard_Complex_Laur_JacoMats;

package Double_Lseries_Polynomials is

-- DESCRIPTION :
--   A Laurent series polynomial is a polynomial in several variables
--   where the coefficients are Laurent series, represented by
--   (1) a vector of leading exponents of the Laurent series coefficients,
--   (2) coefficient vectors of the Laurent series for each monomial;
--   (3) exponent vectors of the monomials in the polynomial.

-- DATA STRUCTURES :

  type Table is record
    nbt : integer32; -- number of terms
    nvr : integer32; -- number of variables
    lead : Standard_Integer_Vectors.Link_to_Vector; -- leading exponents
    cffs : Standard_Complex_VecVecs.Link_to_VecVec; -- coefficient vectors
    mons : Standard_Integer_VecVecs.Link_to_VecVec; -- monomial exponents
  end record;

  type Table_Vector ( nbt : integer32 ) is record -- nbt is number of tables
    nvr : integer32; -- number of variables in each table
    lead : Standard_Integer_VecVecs.VecVec(1..nbt);   -- leading exponents
    cffs : Standard_Complex_VecVecs.Array_of_VecVecs(1..nbt); -- coefficients
    mons : Standard_Integer_VecVecs.Array_of_VecVecs(1..nbt); -- monomials
  end record;

  type Link_to_Table_Vector is access Table_Vector;

  type Table_Vector_Array is
    array ( integer32 range <> ) of Link_to_Table_Vector;

-- CONSTRUCTORS :

  function Random_Table ( dim,nbr,deg,pwr,low,upp : integer32 ) return Table;

  -- DESCRIPTION :
  --   Makes a random polynomial with Laurent series as coefficients.
  --   Wraps the Make_Random_polynomial procedure.

  -- ON ENTRY :
  --   dim      the dimension is the number of variables;
  --   nbr      number of monomials in the polynomial;
  --   deg      degree of the series;
  --   pwr      largest power for every variable;
  --   low      lower bound on leading exponents of the series;
  --   upp      upper bound on leading exponents of the series.

  function Make_Table ( p : Poly; dim,nvr,tdx,deg : integer32 ) return Table;

  -- DESCRIPTION :
  --   Given a Laurent polynomial p and the index for t,
  --   makes the data structures to represent p as a polynomial
  --   with Laurent series coefficients in t.
  --   Wraps the Make_Series_Polynomial procedure.

  -- ON ENTRY :
  --   p       a Laurent polynomial with possibly negative exponents in t;
  --   dim     total number of variables in p, including t,
  --   nvr     number of variables without t;
  --   tdx     index of t as one of the dim variables in p,
  --           if tdx is zero, then dim must equal nvr,
  --           otherwise nvr = dim - 1.

  function Make_Table_Vector
             ( p : Laur_Sys; dim,nvr,tdx,deg : integer32 )
             return Table_Vector;

  -- DESCRIPTION :
  --   Given a Laurent polynomial system p, returns the triplet of data
  --   to represent p as a system with Laurent series coefficients.
  --   Wraps the Make_Series_System procedure.

  -- ON ENTRY :
  --   p       a Laurent system with possibly negative exponents in t;
  --   dim     total number of variables in p, including t,
  --   nvr     number of variables without t;
  --   tdx     index of t as one of the dim variables in p,
  --           if tdx is zero, then dim must equal nvr,
  --           otherwise nvr = dim - 1.

  function Make_Table_Vector_Array
             ( jp : Jaco_Mat; tdx,deg : integer32 )
             return Table_Vector_Array;

  -- DESCRIPTION :
  --   Returns an array of table vectors for the symbolic Jacobian matrix,
  --   skipping the column with index tdx.
  --   The degree of the series equals deg.

-- EVALUATORS :

  procedure Eval ( deg : in integer32; tab : in Table;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ye : out integer32;
                   yc : out Standard_Complex_Vectors.Vector;
                   verbose : in boolean := true );

  -- DESCRIPTION :
  --   Evaluates a Laurent series polynomial at a Laurent series.
  --   Wraps the Eval procedure.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   tab      data for a Laurent series polynomial;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation;
  --   verbose  is the verbose flag.

  -- ON RETURN :
  --   ye       leading exponent of the result of the evaluation;
  --   yc       coefficient vector of the value of the polynomial.

  procedure Eval ( deg : in integer32; tab : in Table_Vector;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   ylead : out Standard_Integer_Vectors.Vector;
                   ycffs : out Standard_Complex_VecVecs.VecVec;
                   verbose : in boolean := true );

  -- DESCRIPTION :
  --   Evaluates a Laurent series polynomial vector at a Laurent series.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   tab      data for a Laurent series polynomial;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation;
  --   ylead    must include the range for the number of tables;
  --   ycffs    must also include the range for the number of tables;
  --   verbose  is the verbose flag.

  -- ON RETURN :
  --   ylead    leading exponents of the result of the evaluation;
  --   ycffs    coefficient vectors of the evaluation result;
  --            allocations will be made as needed.

  procedure Eval ( deg : in integer32; tab : in Table_Vector_Array;
                   xlead : in Standard_Integer_Vectors.Vector;
                   xcffs : in Standard_Complex_VecVecs.Link_to_VecVec;
                   Alead : in out Standard_Integer_Matrices.Matrix;
                   Acffs : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                   verbose : in boolean := true );

  -- DESCRIPTION :
  --   Evaluates a Laurent series Jacobian matrix at a Laurent series.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   tab      data for a Laurent series Jacobian matrix;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation;
  --   Acffs    allocated space for all rows and columns of Laurent
  --            series of degree deg for the evaluated matrix;
  --   verbose  is the verbose flag.

  -- ON RETURN :
  --   Alead    leading exponents of the evaluated Jacobian matrix;
  --   Acffs    Jacobian matrix evaluated at (xlead, xcffs).

-- OUTPUT :

  procedure Write ( tab : in Table; s : in string := "p" );
  procedure Write ( tab : in Table_Vector; s : in string := "p" );

  -- DESCRIPTION :
  --   Writes the table or table vector tab.

  procedure Write ( tva : in Table_Vector_Array; s : in string := "p" );

  -- DESCRIPTION :
  --   Writes the table vector array.

-- BASIC OPERATIONS :

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
                   yc : out Standard_Complex_Vectors.Vector;
                   verbose : in boolean := true );

  -- DESCRIPTION :
  --   Evaluates a polynomial at a Laurent series.

  -- ON ENTRY :
  --   deg      only coefficients in the range 0..deg are considered;
  --   plead    leading exponents of the Laurent series coefficients;
  --   pcffs    coefficient vectors of the Laurent series coefficients;
  --   pmons    exponents of the monomials in the polynomial;
  --   xlead    leading exponents of the argument for the evaluation;
  --   xcffs    coefficient vectors of the argument for the evaluation;
  --   verbose  if true, then the value of every monomial is written,
  --            and also the value after each update.

  -- ON RETURN :
  --   ye       leading exponent of the result of the evaluation;
  --   yc       coefficient vector of the value of the polynomial.

  function tsymbol_index return integer32;

  -- DESCRIPTION :
  --   Returns the index of the symbol t used in the Laurent series
  --   in the input Laurent polynomials from the symbol table contents.
  --   Returns 0 if the symbol table does not contain the symbol t.

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

  function Minimum_Laurent_Series_Degree
             ( p : Poly; tdx : integer32 ) return integer32;
  function Minimum_Laurent_Series_Degree
             ( p : Laur_Sys; tdx : integer32 ) return integer32;

  -- DESCRIPTION :
  --   Given a nonzero index tdx for the position of t in p,
  --   returns the smallest value for the degree in the range of coefficients
  --   so no coefficient of p is lost.

  -- REQUIRED : tdx /= 0.

  procedure Make_Series_Polynomial
              ( p : in Poly; dim,nvr,tdx,deg : in integer32;
                lead : out Standard_Integer_Vectors.Link_to_Vector;
                cffs : out Standard_Complex_VecVecs.Link_to_VecVec;
                mons : out Standard_Integer_VecVecs.Link_to_VecVec );

  -- DESCRIPTION :
  --   Given a Laurent polynomial p and the index for t,
  --   makes the data structures to represent p as a polynomial
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
  --   cffs    coefficients in the Laurent series for each monomial;
  --   mons    exponents of the monomials.

  procedure Make_Series_System
              ( p : in Laur_Sys; dim,nvr,tdx,deg : in integer32;
                lead : out Standard_Integer_VecVecs.VecVec;
                cffs : out Standard_Complex_VecVecs.Array_of_VecVecs;
                mons : out Standard_Integer_VecVecs.Array_of_VecVecs );

  -- DESCRIPTION :
  --   Given a Laurent polynomial system p, returns the triplet of data
  --   to represent p as a system with Laurent series coefficients.

  -- ON ENTRY :
  --   p       a Laurent system with possibly negative exponents in t;
  --   dim     total number of variables in p, including t,
  --   nvr     number of variables without t;
  --   tdx     index of t as one of the dim variables in p,
  --           if tdx is zero, then dim must equal nvr,
  --           otherwise nvr = dim - 1.

  -- ON RETRUN :
  --   lead    leading exponents of the Laurent series coefficients;
  --   cffs    coefficients in the Laurent series for each monomial;
  --   mons    exponents of the monomials.

end Double_Lseries_Polynomials;
