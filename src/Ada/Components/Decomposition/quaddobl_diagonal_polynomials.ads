with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with QuadDobl_Complex_VecVecs;          use QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Systems;     use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;     use QuadDobl_Complex_Laur_Systems;
with Permutations;                      use Permutations;

package QuadDobl_Diagonal_Polynomials is

-- DESCRIPTION :
--   The functions in this package help in the construction of diagonal
--   homotopies to intersect pure dimensional varieties, in the extrinic
--   version, with double double precision arithmetic.

  function Create ( n,i : integer32 ) return QuadDobl_Complex_Polynomials.Term;
  function Create ( n,i : integer32 ) return QuadDobl_Complex_Laurentials.Term;
  function Create ( n,i : integer32 ) return QuadDobl_Complex_Polynomials.Poly;
  function Create ( n,i : integer32 ) return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Returns the i-th variable as a term or polynomial in n variables.
  --   This is used to make a slack variable in the homotopy vanish.

  function Insert_Variables
             ( n : integer32; t : QuadDobl_Complex_Polynomials.Term )
             return QuadDobl_Complex_Polynomials.Term;
  function Insert_Variables
             ( n : integer32; t : QuadDobl_Complex_Laurentials.Term )
             return QuadDobl_Complex_Laurentials.Term;
  function Insert_Variables
             ( n : integer32; p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly;
  function Insert_Variables
             ( n : integer32; p : QuadDobl_Complex_Laurentials.Poly )
             return QuadDobl_Complex_Laurentials.Poly;
  function Insert_Variables ( n : integer32; p : Poly_Sys ) return Poly_Sys;
  function Insert_Variables ( n : integer32; p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Inserts n variables to the term t, or to every term in p.

  function Append_Variables
             ( n : integer32; t : QuadDobl_Complex_Polynomials.Term )
             return QuadDobl_Complex_Polynomials.Term;
  function Append_Variables
             ( n : integer32; t : QuadDobl_Complex_Laurentials.Term )
             return QuadDobl_Complex_Laurentials.Term;
  function Append_Variables
             ( n : integer32; p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly;
  function Append_Variables
             ( n : integer32; p : QuadDobl_Complex_Laurentials.Poly )
             return QuadDobl_Complex_Laurentials.Poly;
  function Append_Variables ( n : integer32; p : Poly_Sys ) return Poly_Sys;
  function Append_Variables ( n : integer32; p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Appends n variables to the term t, or to every term in p.

  function Truncate
             ( t : QuadDobl_Complex_Polynomials.Term; n : integer32 )
             return QuadDobl_Complex_Polynomials.Term;
  function Truncate
             ( t : QuadDobl_Complex_Laurentials.Term; n : integer32 )
             return QuadDobl_Complex_Laurentials.Term;
  function Truncate
             ( p : QuadDobl_Complex_Polynomials.Poly; n : integer32 )
             return QuadDobl_Complex_Polynomials.Poly;
  function Truncate
             ( p : QuadDobl_Complex_Laurentials.Poly; n : integer32 )
             return QuadDobl_Complex_Laurentials.Poly;

  -- DESCRIPTION :
  --   Truncates to the first n variables in term or poly.
  --   Terms with variables of indices > n are omitted in the polynomial.

  function Collapse
             ( t : QuadDobl_Complex_Polynomials.Term; n : integer32 )
             return QuadDobl_Complex_Polynomials.Term;
  function Collapse
             ( t : QuadDobl_Complex_Laurentials.Term; n : integer32 )
             return QuadDobl_Complex_Laurentials.Term;
  function Collapse
             ( p : QuadDobl_Complex_Polynomials.Poly; n : integer32 )
             return QuadDobl_Complex_Polynomials.Poly;
  function Collapse
             ( p : QuadDobl_Complex_Laurentials.Poly; n : integer32 )
             return QuadDobl_Complex_Laurentials.Poly;
  function Collapse ( p : Poly_Sys; n : integer32 ) return Poly_Sys;
  function Collapse ( p : Laur_Sys; n : integer32 ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the term or polynomial after elimination of the diagonal.
  --   Terms with slack variables in them are omitted in the polynomial.

  function Collapse
             ( t : QuadDobl_Complex_Polynomials.Term;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Polynomials.Term;
  function Collapse
             ( t : QuadDobl_Complex_Laurentials.Term;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Laurentials.Term;
  function Collapse
             ( p : QuadDobl_Complex_Polynomials.Poly;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Polynomials.Poly;
  function Collapse
             ( p : QuadDobl_Complex_Laurentials.Poly;
               n : integer32; q : Permutation )
             return QuadDobl_Complex_Laurentials.Poly;
  function Collapse ( p : Poly_Sys; n : integer32; q : Permutation )
                    return Poly_Sys;
  function Collapse ( p : Laur_Sys; n : integer32; q : Permutation )
                    return Laur_Sys;

  -- DESCRIPTION :
  --   Returns a term, polynomial, or system in n variables, omitting the
  --   slack variables, and collapsing the second group of variables onto
  --   the first one, using the permutation q.

  function Diagonal ( n : integer32 ) return Poly_Sys;
  function Diagonal ( n : integer32 ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns the equations of the diagonal system x(i) - y(i) = 0,
  --   for i from 1 to n, so Diagonal(n) has 2*n variables.

  function Product ( n1,n2 : integer32; p1,p2 : Poly_Sys ) return Poly_Sys;
  function Product ( n1,n2 : integer32; p1,p2 : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   The system on return has n1+n2 variables and consists of the
  --   equations in p1 (with n2 variables added to each polynomial) and
  --   the equations in p2 (with n1 variables inserted to each polyomial).

  function Product ( n,k : integer32; hyp1,hyp2 : VecVec ) return VecVec;

  -- DESCRIPTION :
  --   Returns k vectors, taken from the first n variables from hyp1
  --   and hyp2, and as constant the sum of both hyperplanes.

end QuadDobl_Diagonal_Polynomials;
