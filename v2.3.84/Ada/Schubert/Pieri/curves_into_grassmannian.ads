with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Brackets;                           use Brackets;

package Curves_into_Grassmannian is

-- DESCRIPTION :
--   Provides tools to deal with curves of degree q that produce
--   p-planes in n-dimensional space, where n = m+p.
--   The curves are represented symbolically by n-by-p matrices
--   of homogeneous polynomials in s and t.  This symbolic
--   representation corresponds to a localization pattern.
--   The degree q equals the degree of the maximal minors.
--   The numerical coefficients are stored in (m+p)*(q+1)-by-p
--   matrices with top and bottom pivots of the embedding.
--   Localization maps are matrices of the same format that
--   indicate where to put the ones to fix a chart to compute in.

-- CREATOR :

  function Symbolic_Create ( m,p,q : natural32; top,bottom : Bracket )
                           return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Creates the symbolic form of a q-curve that produces p-planes
  --   into (m+p)-dimensional space, as prescribed by the localization
  --   pattern defined by top and bottom pivots.
  --   The variable s is the next-to-last and t is the last variable
  --   in the polynomials of the matrix on return.

-- SELECTORS :

  function Number_of_Variables ( top,bottom : Bracket ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of x_ij-variables needed to represent the matrix
  --   of polynomials as prescribed in the localization pattern defined
  --   by the top and bottom pivots.  The x_ij's are the indeterminate
  --   coefficients of the polynomials in the matrix.
  --   This function determines the degrees of freedom in the poset.
  --   Note that the dimension of the corresponding space is p less because
  --   we can scale the columns independently dividing by a nonzero constant.

  function Standard_Coordinate_Frame
             ( m,p,q : natural32; top,bottom : Bracket;
               coeff : Standard_Complex_Matrices.Matrix )
             return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the localization map for the curve as (m+p)*(q+1)-by-q matrix, 
  --   taking the first nonzero occuring coefficients as places for the ones.

  function Eval ( c : Standard_Complex_Poly_Matrices.Matrix;
                  s,t : Complex_Number )
                return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Evaluates the curve at the point (s,t), returns a matrix
  --   of indeterminate coefficients.
  --   The number of variables in the polynomials remains the same.

  function Elim ( c : Standard_Complex_Poly_Matrices.Matrix;
                  s,t : Complex_Number )
                return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Evaluates the curve at the point (s,t), returns a matrix
  --   of indeterminate coefficients and eliminates s and t.

  function Substitute ( c : Standard_Complex_Poly_Matrices.Matrix;
                        v : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Substitutes the indeterminate coefficients in the matrix c
  --   by the values in the vector v and returns a polynomial matrix
  --   in the variables that are left free.

  -- REQUIRED :
  --   The number of unknowns in c is larger than or equal to v'length.

  function Convert ( c : Standard_Complex_Poly_Matrices.Matrix )
                   return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Converts a matrix of constant polynomials into a regular matrix
  --   of standard complex numbers.

  function Column_Localize ( top,bottom : Bracket;
                             locmap : Standard_Natural_Matrices.Matrix;
                             pm : Standard_Complex_Poly_Matrices.Matrix )
                           return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Applies the localization map to the matrix pm.

  function Column_Localize ( top,bottom : Bracket;
                             locmap : Standard_Natural_Matrices.Matrix;
                             p : Poly_Sys ) return Poly_Sys; 

  -- DESCRIPTION :
  --   Applies the localization map to the polynomial system p,
  --   assuming that the variables are ordered columnwise.

  function Column_Vector_Rep
             ( top,bottom : Bracket;
               cffmat : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector representation of the coefficient matrix,
  --   assuming the variables are stored columnwise.
  --   This is for nonlocalized coefficient matrices.

  function Column_Vector_Rep ( locmap : Standard_Natural_Matrices.Matrix;
                               cffmat : Standard_Complex_Matrices.Matrix )
                             return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector representation of a localized coefficient matrix.

  function Column_Matrix_Rep
             ( locmap : Standard_Natural_Matrices.Matrix;
               cffvec : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Converts from vector to matrix representation of the coefficients,
  --   assuming the variables are stored columnwise in the map.

-- MODIFIER :

  procedure Swap ( c : in out Standard_Complex_Poly_Matrices.Matrix;
                   k,l : in integer32 );

  -- DESCRIPTION :
  --   Swaps the variables k and l in the polynomial matrix.

  function Insert ( c : Standard_Complex_Poly_Matrices.Matrix; k : integer32 )
                  return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   The map on return has a k-th variable inserted with zero power.

end Curves_into_Grassmannian;
