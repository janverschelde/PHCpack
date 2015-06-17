with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;            use Standard_Natural_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;

package QuadDobl_Monomial_Evaluations is

-- DESCRIPTION :
--   Offers routines to evaluate sequences of monomials efficiently.
--   A monomial is defined as a vector of natural numbers.

  function Eval ( e : Standard_Natural_Vectors.Vector;
                  x : QuadDobl_Complex_Vectors.Vector )
                return Complex_Number;

  -- DESCRIPTION :
  --   Returns the product of x(i)^e(i), for i in e'range.

  -- REQUIRED : x'range = e'range.

  function Eval ( e : Standard_Natural_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRITPION :
  --   Returns the values of the monomials defined by e at x.

  function Power_Table
                ( n,m : integer32;
                  d : Standard_Natural_Vectors.Vector;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns an n-by-m matrix with consecutive powers of x.
  --   The (i,j)-th entry on return contains x(i)^j for all j <= d(i).

  -- REQUIRED : m >= max(d).

  function Eval ( e : Standard_Natural_Vectors.Vector;
                  p : QuadDobl_Complex_Matrices.Matrix )
                return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the monomial defined by e at x,
  --   using the powers of x in the table in p.
 
  function Eval ( d : Standard_Natural_Vectors.Vector;
                  e : Standard_Natural_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the monomials defined by e at x,
  --   using d = Largest_Degrees(x'last,e).

  function Eval_with_Power_Table 
                ( e : Standard_Natural_VecVecs.VecVec;
                  x : QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the sequence of monomials defined by e at x,
  --   using a table of consecutive powers of x.

end QuadDobl_Monomial_Evaluations;
