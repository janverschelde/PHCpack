with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;
with Brackets;                           use Brackets;

package Curves_into_Grassmannian_io is

-- DESCRIPTION :
--   This package facilitates output of curves of degree q that produce
--   p-planes in n-dimensional space, where n = m+p.

  procedure One_Set_up_Symbol_Table
              ( m,p,q : in natural32; top,bottom : in Bracket );
  procedure Two_Set_up_Symbol_Table
              ( m,p,q : in natural32; top,bottom : in Bracket );

  -- DESCRIPTION :
  --   Set up of the symbol table with enough symbols to represent the
  --   matrix of polynomials as prescribed in the localization pattern.
  --   The polynomials in the matrix are polynomials in "s" and "t" which are
  --   the first two variables in the symbol table.  The other variables are
  --   added columnwise in the pattern defined by top and bottom brackets.
  --   The x_ij's are represented by "xijsk": k-th block at (i,j)-th place.
  --   If the symbol table is not empty on calling this procedure,
  --   then its contents will be destroyed and replaced by the x_ij's.
  --   One_* : top or bottom type; Two_* : mixed type of node.

  procedure Reduce_Symbols ( top,bottom : in Bracket; locmap : in Matrix );

  -- DESCRIPTION :
  --   Removes the symbols that correspond to the ones in the localization
  --   map prescribed by top and bottom pivots.

end Curves_into_Grassmannian_io;
