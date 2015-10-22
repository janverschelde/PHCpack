with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Matrices;
with Brackets;                          use Brackets;
with Standard_Bracket_Polynomials;
with DoblDobl_Bracket_Polynomials;
with QuadDobl_Bracket_Polynomials;
with Remember_Numeric_Minors;           use Remember_Numeric_Minors;
with Remember_Symbolic_Minors;          use Remember_Symbolic_Minors;

package Numeric_Schubert_Conditions is

-- DESCRIPTION :
--   Given the symbolic bracket systems that represent Schubert conditions
--   as minors of a matrix of type [ X | F ],
--   this package formulates the corresponding polynomial systems.

  function Degree ( b : Bracket; k : natural32 ) return integer32;

  -- DESCRIPTION :
  --   Returns the number of entries in b that are <= k.
  --   This number determines the degree of the polynomial equations
  --   resulting from the intersection conditions det([ X | F ]) = 0,
  --   if b is the selection of columns of [ X | F ].

  function Permute ( b,p : Bracket ) return Bracket;

  -- DESCRIPTION :
  --   Returns p(b), the entries in b are filtered through p.

  function Substitute
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial;
               t : Standard_Numeric_Minors )
             return Standard_Bracket_Polynomials.Bracket_Polynomial;
  function Substitute
             ( p : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : DoblDobl_Numeric_Minors )
             return DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
  function Substitute
             ( p : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : QuadDobl_Numeric_Minors )
             return QuadDobl_Bracket_Polynomials.Bracket_Polynomial;

  -- DESCRIPTION :
  --   Substitutes the second bracket of each monomial in p
  --   by the corresponding value in the numeric minor table, computed in
  --   standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p        bracket polynomial encodes a Laplace expansion,
  --            with two brackets in each monomial, respectively
  --            the symbolic and numeric minor;
  --   t        remember table for numerical minors.

  function Substitute
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial;
               t : Standard_Numeric_Minors; rows : Bracket )
             return Standard_Bracket_Polynomials.Bracket_Polynomial;
  function Substitute
             ( p : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : DoblDobl_Numeric_Minors; rows : Bracket )
             return DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
  function Substitute
             ( p : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : QuadDobl_Numeric_Minors; rows : Bracket )
             return QuadDobl_Bracket_Polynomials.Bracket_Polynomial;

  -- DESCRIPTION :
  --   Substitutes the second bracket of each monomial in p
  --   by the corresponding value in the numeric minor table, computed in
  --   standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   p       bracket polynomial encodes a Laplace expansion,
  --           with two brackets in each monomial, respectively
  --           the symbolic and numeric minor;
  --   t       remember table for numerical minors;
  --   rows    optional displacement of the rows selected by the brackets
  --           in the Laplace expansion.

  function Substitute
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial;
               t : Standard_Symbolic_Minors )
             return Standard_Complex_Polynomials.Poly;
  function Substitute
             ( p : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : DoblDobl_Symbolic_Minors )
             return DoblDobl_Complex_Polynomials.Poly;
  function Substitute
             ( p : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : QuadDobl_Symbolic_Minors )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Substitutes the minors in p by the corresponding polynomials of t,
  --   in standard double, double double, or quad double precision.

  function Substitute
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial;
               t : Standard_Symbolic_Minors; rows : Bracket )
             return Standard_Complex_Polynomials.Poly;
  function Substitute
             ( p : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : DoblDobl_Symbolic_Minors; rows : Bracket )
             return DoblDobl_Complex_Polynomials.Poly;
  function Substitute
             ( p : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               t : QuadDobl_Symbolic_Minors; rows : Bracket )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Substitutes the minors in p by the corresponding polynomials of t,
  --   in standard double, double double, or quad double precision,
  --   taking the displacement with the rows into account.

  function Substitute
             ( p : Standard_Bracket_Polynomials.Bracket_Polynomial;
               nt,st : Standard_Symbolic_Minors; rows : Bracket )
             return Standard_Complex_Polynomials.Poly;
  function Substitute
             ( p : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               nt,st : DoblDobl_Symbolic_Minors; rows : Bracket )
             return DoblDobl_Complex_Polynomials.Poly;
  function Substitute
             ( p : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               nt,st : QuadDobl_Symbolic_Minors; rows : Bracket )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Substitutes minors in p by the corresponding polynomials in nt,
  --   computed in standard double, double double, or quad double precision,
  --   for the first brackets and st for the second brackets of every
  --   term in p.

  function Select_Columns
             ( A : Standard_Complex_Matrices.Matrix;
               col : Bracket; d,k : integer32 )
             return Standard_Complex_Matrices.Matrix;
  function Select_Columns
             ( A : DoblDobl_Complex_Matrices.Matrix;
               col : Bracket; d,k : integer32 )
             return DoblDobl_Complex_Matrices.Matrix;
  function Select_Columns
             ( A : QuadDobl_Complex_Matrices.Matrix;
               col : Bracket; d,k : integer32 )
             return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Selects as many as d columns of A, a matrix with coefficients
  --   in standard double, double double, or quad double precision,
  --   where d = col'last - Degree(b,k), for some k-bracket b.

  function Select_Columns
             ( A : Standard_Complex_Poly_Matrices.Matrix;
               col : Bracket; d,k : integer32 )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Select_Columns
             ( A : DoblDobl_Complex_Poly_Matrices.Matrix;
               col : Bracket; d,k : integer32 )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Select_Columns
             ( A : QuadDobl_Complex_Poly_Matrices.Matrix;
               col : Bracket; d,k : integer32 )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Selects as many as d columns of A, a matrix of polynomials
  --   with coefficients in standard double, double double,
  --   or quad double precision, where d = col'last - Degree(b,k),
  --   for some k-bracket b.

  function Select_Columns
              ( x : Standard_Complex_Poly_Matrices.Matrix;
                col : Bracket; d : integer32 )
              return Standard_Complex_Poly_Matrices.Matrix;
  function Select_Columns
              ( x : DoblDobl_Complex_Poly_Matrices.Matrix;
                col : Bracket; d : integer32 )
              return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Select_Columns
              ( x : QuadDobl_Complex_Poly_Matrices.Matrix;
                col : Bracket; d : integer32 )
              return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Selects as many as d columns of x, a matrix of polynomials with
  --   coefficients in standard double, double double, or quad double
  --   precision, using the first d indices from col.

  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               A : Standard_Complex_Matrices.Matrix ) 
             return Standard_Bracket_Polynomials.Bracket_Polynomial;
  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               A : DoblDobl_Complex_Matrices.Matrix ) 
             return DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               A : QuadDobl_Complex_Matrices.Matrix ) 
             return QuadDobl_Bracket_Polynomials.Bracket_Polynomial;

  -- DESCRIPTION :
  --   Applies Laplace expansion to one minor defined by row and col
  --   to express that the intersection of a k-plane with an f-plane,
  --   in standard double, double double, or quad double precision,
  --   generated by the columns of A, in n-space has dimension i.

  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               X : Standard_Complex_Poly_Matrices.Matrix; 
               A : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Polynomials.Poly;
  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               X : DoblDobl_Complex_Poly_Matrices.Matrix; 
               A : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Polynomials.Poly;
  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               X : QuadDobl_Complex_Poly_Matrices.Matrix; 
               A : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Applies Laplace expansion to one minor defined by row and col
  --   to express that the intersection of a k-plane with an f-plane,
  --   generated by the columns of A, in n-space has dimension i,
  --   in standard double, double double, or quad double precision.
  --   The matrix X is a matrix of polynomials with coefficients in
  --   standard double, double double, or quad double precision.

  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               X,A : Standard_Complex_Poly_Matrices.Matrix)
             return Standard_Complex_Polynomials.Poly;
  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               X,A : DoblDobl_Complex_Poly_Matrices.Matrix)
             return DoblDobl_Complex_Polynomials.Poly;
  function Laplace_One_Minor
             ( n,k : integer32; row,col : Bracket;
               X,A : QuadDobl_Complex_Poly_Matrices.Matrix)
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Applies Laplace expansion to one minor defined by row and col
  --   to express that the intersection of a k-plane with an f-plane,
  --   generated by the columns of A, in n-space has dimension i.
  --   Both X and A are represented as polynomials with coefficients
  --   in standard double, double double, or quad double precision.

  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : Standard_Bracket_Polynomials.Bracket_Polynomial;
               A : Standard_Complex_Matrices.Matrix )
             return Standard_Bracket_Polynomials.Bracket_Polynomial;
  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               A : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               A : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Bracket_Polynomials.Bracket_Polynomial;

  -- DESCRIPTION :
  --   Retrieves the row and column from the one monomial in fm
  --   and computes the condition that in n-space, a k-plane meets
  --   an f-plane, generated by the columns of A, in an i-space,
  --   where the matrices that represent the f-plane are given in
  --   standard double, double double, or quad double precision.

  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : Standard_Bracket_Polynomials.Bracket_Polynomial;
               X : Standard_Complex_Poly_Matrices.Matrix;
               A : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Polynomials.Poly;
  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               X : DoblDobl_Complex_Poly_Matrices.Matrix;
               A : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Polynomials.Poly;
  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               X : QuadDobl_Complex_Poly_Matrices.Matrix;
               A : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Retrieves the row and column from the one monomial in fm
  --   and computes the condition that in n-space, a k-plane meets
  --   an f-plane, generated by the columns of A, in an i-space,
  --   where the matrices that represent the f-plane are given in
  --   standard double, double double, or quad double precision.
  --   The matrices X are polynomials with coefficients in
  --   standard double, double double, or quad double precision.

  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : Standard_Bracket_Polynomials.Bracket_Polynomial;
               X,A : Standard_Complex_Poly_Matrices.Matrix )
             return Standard_Complex_Polynomials.Poly;
  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
               X,A : DoblDobl_Complex_Poly_Matrices.Matrix )
             return DoblDobl_Complex_Polynomials.Poly;
  function Elaborate_One_Flag_Minor
             ( n,k,f,i : integer32;
               fm : QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
               X,A : QuadDobl_Complex_Poly_Matrices.Matrix )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Retrieves the row and column from the one monomial in fm
  --   and computes the condition that in n-space, a k-plane meets
  --   an f-plane, generated by the columns of A, in an i-space.
  --   The matrices X and A are polynomials with coefficients in
  --   standard double, double double, or quad double precision.

  function Expand ( n,k,nq : integer32; lambda : Bracket;
                    X : Standard_Complex_Poly_Matrices.Matrix;
                    flag : Standard_Complex_Matrices.Matrix )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Expand ( n,k,nq : integer32; lambda : Bracket;
                    X : DoblDobl_Complex_Poly_Matrices.Matrix;
                    flag : DoblDobl_Complex_Matrices.Matrix )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Expand ( n,k,nq : integer32; lambda : Bracket;
                    X : QuadDobl_Complex_Poly_Matrices.Matrix;
                    flag : QuadDobl_Complex_Matrices.Matrix )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DECRIPTION :
  --   Expands all symbolic and numeric minors to return a polynomial
  --   system of nq equations for the Schubert conditions in lambda on
  --   a k-plane X in n-space meeting the given flag,
  --   in standard double, double double, or quad double precision.

  function Expand ( n,k,nq : integer32; lambda : Bracket;
                    X,F : Standard_Complex_Poly_Matrices.Matrix )
                  return Standard_Complex_Poly_Systems.Poly_Sys;
  function Expand ( n,k,nq : integer32; lambda : Bracket;
                    X,F : DoblDobl_Complex_Poly_Matrices.Matrix )
                  return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Expand ( n,k,nq : integer32; lambda : Bracket;
                    X,F : QuadDobl_Complex_Poly_Matrices.Matrix )
                  return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   All minors, both in X and F are now polynomials for expansion
  --   into a polynomial system for Schubert conditions on a k-plane in X
  --   meeting a flag F, with computations for the expansion done in
  --   standard_double, double double, or quad double precision.

  -- REQUIRED :
  --   The number of unknowns in all polynomial matrices is the same for all.

-- MINIMAL REPRESENTATION of Schubert Problems :

  function Minimal_Expand
             ( n,k,nq : integer32; lambda : Bracket;
               X : Standard_Complex_Poly_Matrices.Matrix;
               flag : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Minimal_Expand
             ( n,k,nq : integer32; lambda : Bracket;
               X : DoblDobl_Complex_Poly_Matrices.Matrix;
               flag : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Minimal_Expand
             ( n,k,nq : integer32; lambda : Bracket;
               X : QuadDobl_Complex_Poly_Matrices.Matrix;
               flag : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DECRIPTION :
  --   Expands all symbolic and numeric minors to return a polynomial
  --   system of nq equations for the Schubert conditions in lambda on
  --   a k-plane X in n-space meeting the given flag.
  --   This uses a more efficient formulation of the Schubert problem.

  -- ON ENTRY :
  --   n       ambient dimension;
  --   k       dimension of the solution plane;
  --   nq      total number of equations, which must be computed with
  --           Symbolic_Schubert_Conditions.Number_of_NotAbove(k,lambda);
  --   lambda  conditions imposed by the flag;
  --   X       matrix of indeterminates along some localization pattern;
  --   flag    numerical values for the coordinates in the flag.

-- WRAPPER FUNCTIONS :

  function Expanded_Polynomial_Equations 
             ( n,k : integer32; cond : Bracket;
               flag : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Poly_Systems.Poly_Sys;
  function Expanded_Polynomial_Equations 
             ( n,k : integer32; cond : Bracket;
               flag : DoblDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys;
  function Expanded_Polynomial_Equations 
             ( n,k : integer32; cond : Bracket;
               flag : QuadDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys;

  -- DESCRIPTION :
  --   Returns the expanded polynomial equation for the condition on
  --   a k-plane in n-space with the localization map in locmap,
  --   imposed by the conditions in cond and with respect to the given flag,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n        ambient dimension;
  --   k        dimension of the solution plane;
  --   cond     conditions imposed by Schubert conditions;
  --   flag     used for formulating the polynomial equations.

  -- ON RETURN :
  --   Expanded polynomial equations that define the Schubert conditions
  --   via all expanded minors.

end Numeric_Schubert_Conditions;
