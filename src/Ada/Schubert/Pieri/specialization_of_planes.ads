with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Brackets;                           use Brackets;
with Standard_Complex_Matrices;
with Standard_Complex_Poly_Matrices;

package Specialization_of_Planes is

-- DESCRIPTION :
--   Set up of the moving cycles in the Pieri homotopy algorithm.
--   The U-matrix generates the special m-plane.

  function Random_Upper_Triangular
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns an n*n random upper triangular matrix that has 1's on its
  --   anti-diagonal, to be used for the first Pieri tree.

  function Random_Lower_Triangular
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns an n*n random lower triangular matrix that has 1's on its
  --   diagonal, to be used for the second Pieri tree.

  function U_Matrix ( F : Standard_Complex_Matrices.Matrix; b : Bracket )
                    return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Determines the matrix U in the moving cycles for the Pieri tree.

  -- ON ENTRY :
  --   F         general triangular matrix, output of either
  --             Random_Upper_Triangular or Random_Lower_Triangular;
  --   b         bracket of first jumping-branching node down the tree.

  function Special_Plane ( m : integer32; b : Bracket )
                         return Standard_Complex_Matrices.Matrix;   

  -- DESCRIPTION :
  --   Generates the special m-plane that every child of the node with
  --   the pivots in the bracket b intersects.
  --   The plane on return is spanned by the standard basis vectors
  --   except those that are indexed by the entries in the bracket b.

  function Special_Bottom_Plane ( m : integer32; b : Bracket )
                                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Same as special m-plane, but with random numbers above the diagonal.

  function Special_Top_Plane ( m : integer32; b : Bracket )
                             return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Same as special m-plane, but with random numbers below the diagonal.

  function Special_Plane
                ( n,m,k : integer32; b : Bracket;
                  special : in Standard_Complex_Matrices.Matrix )
                return Standard_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a special (m+1-k)-plane, by random combinations from the
  --   n-dimensional basis vectors of special indexed by the pivots in b.
  --   The basis vectors have random numbers above the diagonal.
  --   Note that here not the complement of the indices in b is used!

  function Moving_U_Matrix
               ( n : integer32; U,L : Standard_Complex_Matrices.Matrix ) 
               return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Determines the moving U-matrix for the Pieri tree, when i=0.
  --   Returns a polynomial matrix in the continuation parameter t,
  --   which is the last variable in the polynomials on return.

  -- ON ENTRY :
  --   n         number of variables of the polynomials in the matrix;
  --   U         output of U_Matrix function listed above, start m-plane;
  --   L         target m-plane that will be folded in during the deformation.

  function Moving_U_Matrix
               ( U : Standard_Complex_Matrices.Matrix;
                 i,r : integer32; b : bracket ) 
               return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   Determines the moving U-matrix for the Pieri tree, when i>0.
  --   Returns a polynomial matrix in the continuation parameter t.

  -- ON ENTRY :
  --   U         output of U_Matrix function listed above;
  --   i         counts the number of nodes till 1st jumping node down;
  --   b         bracket of first jumping-branching node down the tree.

  function Lower_Section
             ( M : Standard_Complex_Poly_Matrices.Matrix;
               row : integer32 ) return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   The columns in M that have nonzero entries in the rows strictly
  --   larger than the given row index will be removed.  This corresponds
  --   to intersecting the column space of M with <e_1,..,e_row>.

  function Upper_Section
             ( M : Standard_Complex_Poly_Matrices.Matrix;
               row : integer32 ) return Standard_Complex_Poly_Matrices.Matrix;

  -- DESCRIPTION :
  --   The columns in M that have nonzero entries in the rows strictly
  --   lower than the given row index will be removed.  This corresponds
  --   to intersecting the column space of M with <e_row,..,e_n>.

end Specialization_of_Planes;
