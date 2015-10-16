with Symbol_Table;                       use Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;

package Matrix_Indeterminates is

-- DESCRIPTION :
--   This package provides a facility to manipulate symbols for the
--   indeterminates xij of a matrix, where i and j are allowed to run
--   in the range 1..F, along the usual hexadecimal coding.

  procedure Initialize_Symbols ( n,d : in natural32 );

  -- DESCRIPTION :
  --   Initializes the symbol table with the variable order
  --     x11 > x12 > .. > x1d > x21 > x22 > .. > x2d > .. > xn1 > .. > xnd.

  function Dimension ( locmap : Matrix ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of free variables in the localization map,
  --   counted by the number of "2" values in the matrix.

  procedure Initialize_Symbols
               ( d : in natural32; locmap : in Matrix );
  procedure Initialize_Symbols ( locmap : in Matrix );

  -- DESCRIPTION :
  --   Initializes the symbol table with those variables xij
  --   as indicated by the "2" positions in the localization map.
  --   The d on entry equals Dimension(locmap).

  function X_ij ( i,j : natural32 ) return Symbol;

  -- DESCRIPTION : 
  --   Returns the symbol that represents the variable xij.

  function Monomial ( n,d,i,j : natural32 )
                    return Standard_Complex_Polynomials.Poly;
  function Monomial ( n,d,i,j : natural32 )
                    return DoblDobl_Complex_Polynomials.Poly;
  function Monomial ( n,d,i,j : natural32 )
                    return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns the representation of X_ij as a polynomial.

  procedure Reduce_Symbols ( locmap : in Matrix );

  -- DESCRIPTION :
  --   Reduces the number of symbols in the symbol table, removing all 
  --   symbols that correspond to zeros and ones in the localization map.

  procedure Clear_Symbols;

  -- DESCRIPTION :
  --   Destruction of the symbol table.

end Matrix_Indeterminates;
