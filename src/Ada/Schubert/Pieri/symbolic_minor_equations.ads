with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Matrices;
with Brackets;                           use Brackets;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Systems;           use Standard_Bracket_Systems;

package Symbolic_Minor_Equations is

-- DESCRIPTION :
--   This package generates the equations that arise from the intersection
--   conditions along the deformations determined by the Pieri trees.
--   The equations are symbolic: the coefficients are brackets.

-- LOCALIZATION PATTERNS :

  function Schubert_Pattern
             ( n : natural32; b1,b2 : Bracket )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Schubert_Pattern
             ( n : natural32; b1,b2 : Bracket )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Schubert_Pattern
             ( n : natural32; b1,b2 : Bracket )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTON :
  --   Returns the representation of the pattern of the p-plane that satisfies
  --   the Schubert conditions as a polynomial matrix.
  --   This definition is used in the original Pieri implementation.

  function Localization_Pattern
             ( n : natural32; top,bottom : Bracket )
             return Standard_Complex_Poly_Matrices.Matrix;
  function Localization_Pattern
             ( n : natural32; top,bottom : Bracket )
             return DoblDobl_Complex_Poly_Matrices.Matrix;
  function Localization_Pattern
             ( n : natural32; top,bottom : Bracket )
             return QuadDobl_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the matrix of indeterminates for the top and bottom pivots.
  --   The dimension of the working space equals n.
  --   This definition is used in the second Pieri implementation.

-- SYMBOLIC REPRESENTATIONS OF THE EQUATIONS :

  function Number_of_Maximal_Minors ( n,m : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of maximal minors of an n-by-m matrix,
  --   where n >= m, i.e.: returns n!/(m!*(n-m)!).

  function Number_of_Minors ( n,m,s : natural32 ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of minors of size s of an n-by-m matrix.

  -- REQUIRED : s <= min(n,m).

  function Maximal_Minors ( n,m : natural32 ) return Bracket_Monomial;

  -- DESCRIPTION :
  --   Generates all maximal minors of an n-by-m matrix, n >= m.
  --   The collection of all these m-by-m minors is represented as the
  --   product of m-brackets.  Every bracket determines a selection of
  --   m rows from the n-by-m matrix.

  -- REQUIRED : n >= m.

  function Minors ( n,m,s : natural32 ) return Bracket_Polynomial;

  -- DESCRIPTION :
  --   The bracket polynomial on return encodes all size s minors of
  --   a n-by-m matrix.  Every monomial in the returned polynomial
  --   starts with a row selection, the brackets following the row
  --   selection represent the corresponding selected columns.

  -- REQUIRED : s <= min(n,m).

  function Minor_Equations
             ( m,d : natural32; bm : Bracket_Monomial ) return Bracket_System;

  -- DESCRIPTION :
  --   Returns the bracket polynomials that arise from expanding the
  --   above maximal m-minors into blocks of respective sizes m-d and d.
  --   Equation 0 in the result is the generic bracket representation of
  --   the Laplace expansion.

  function Expanded_Minor
             ( m : Standard_Complex_Poly_Matrices.Matrix; b : Bracket )
             return Standard_Complex_Polynomials.Poly;
  function Expanded_Minor
             ( m : DoblDobl_Complex_Poly_Matrices.Matrix; b : Bracket )
             return DoblDobl_Complex_Polynomials.Poly;
  function Expanded_Minor
             ( m : QuadDobl_Complex_Poly_Matrices.Matrix; b : Bracket )
             return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   The minor that selects rows from m, according to b is expanded.

  function Extend_Zero_Lifting
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Extends every term with a new variable t^0.

end Symbolic_Minor_Equations;
