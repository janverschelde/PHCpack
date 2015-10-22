with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Floating_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;

package Numeric_Minor_Equations is

-- DESCRIPTION :
--   This package evaluates the symbolic equations in the Pieri homotopies.

-- EXPANDING ACCORDING A BRACKET MONOMIAL :

  function Expanded_Minors
               ( cffmat : Standard_Floating_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly;

  function Expanded_Minors
               ( cffmat : Standard_Complex_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly;

  function Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly;

  function Lifted_Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bm : Bracket_Monomial ) return Poly;

  -- DESCRIPTION :
  --   Expansion of coefficient and polynomial minors along the Laplace
  --   expansion formula in bm creates a polynomial equation.
  --   With the prefix Lifted_, the polynomials in polmat are extended
  --   with a zero lifting.

  -- ON ENTRY :
  --   cffmat    coefficient matrix, represents m-plane;
  --   cntmat    polynomial matrix, represents moving m-plane,
  --             the continuation parameter is the last variable;
  --   polmat    polynomial matrix, contains the pattern of the p-plane;
  --   bm        quadratic bracket monomial, the first bracket is a coefficient
  --             minor and has zero as its first entry, the second bracket is
  --             a polynomial minor.

-- EXPANDING ACCORDING A BRACKET POLYNOMIAL :

  function Expanded_Minors
               ( cffmat : Standard_Floating_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly;

  function Expanded_Minors
               ( cffmat : Standard_Complex_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly;

  function Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly;

  function Lifted_Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bp : Bracket_Polynomial ) return Poly;

  -- DESCRIPTION :
  --   Expansion of coefficient and polynomial minors along the Laplace
  --   expansion formula in bp creates a polynomial equation.
  --   With the prefix Lifted_, the polynomials in polmat are extended
  --   with a zero lifting.

  -- ON ENTRY :
  --   cffmat    coefficient matrix, represents m-plane, m = n-p;
  --   cntmat    polynomial matrix, represents moving m-plane, m = n-p,
  --             the continuation parameter is the last variable;
  --   polmat    polynomial matrix, contains the pattern of the p-plane;
  --   bp        Laplace expansion of one minor, the coefficient minors come
  --             first and have a zero as first element.

-- EXPANDING TO CONSTRUCT POLYNOMIAL SYSTEMS :

  function Expanded_Minors
               ( cffmat : Standard_Floating_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys;

  function Expanded_Minors
               ( cffmat : Standard_Complex_Matrices.Matrix;
                 polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys;

  function Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys;

  function Lifted_Expanded_Minors
               ( cntmat,polmat : Standard_Complex_Poly_Matrices.Matrix;
                 bs : Bracket_System ) return Poly_Sys;

  -- DESCRIPTION :
  --   Expansion of coefficient and polynomial minors along the Laplace
  --   expansion formulas in bs creates a polynomial system.
  --   With the prefix Lifted_, the polynomials in polmat are extended
  --   with zero lifting.

  -- ON ENTRY :
  --   cffmat    coefficient matrix, represents m-plane;
  --   cntmat    polynomial matrix, represents moving m-plane,
  --             the continuation parameter is the last variable;
  --   polmat    polynomial matrix, contains the pattern of the p-plane;
  --   bs        Laplace expansion of all minors, the first equation is
  --             the generic one and should not count in the range of
  --             the resulting polynomial system.

  function Evaluate ( p : Poly; x : Standard_Complex_Matrices.Matrix )
                    return Complex_Number;

  -- DESCRIPTION :
  --   Evaluates the polynomial p at the matrix x, where x is a value
  --   for the polynomial matrix used above to define p.

  function Evaluate ( p : Poly_Sys; x : Standard_Complex_Matrices.Matrix )
                    return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Evaluates the polynomial system p at the matrix x, where x is a value
  --   for the polynomial matrix used above to define p.

  procedure Embed ( t : in out Term );
  procedure Embed ( p : in out Poly );
  procedure Embed ( p : in out Poly_Sys );
  procedure Embed ( m : in out Standard_Complex_Poly_Matrices.Matrix );

  -- DESCRIPTION :
  --   Augments the number of variables with one, as is required to embed
  --   the polynomials in a homotopy.

  function Linear_Homotopy ( target,start : Poly ) return Poly;

  -- DESCRIPTION :
  --   Returns (1-t)*start + t*target, with t an additional last variable.

  function Linear_Interpolation
              ( target,start : Poly; k : integer32 ) return Poly;

  -- DESCRIPTION :
  --   Returns (1-t)*start + t*target, with t the k-th variable.

  procedure Divide_Common_Factor ( p : in out Poly; k : in integer32 );

  -- DESCRIPTION :
  --   If the k-th variable occurs everywhere in p with a positive power,
  --   then it will be divided out.

end Numeric_Minor_Equations;
