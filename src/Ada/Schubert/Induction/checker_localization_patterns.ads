with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Matrices;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;

package Checker_Localization_Patterns is

-- DESCRIPTION :
--   This package defines localization patterns,
--   based on the relative positions of black and white checkers.

-- Part I: patterns of the moving flag, determined by black checkers

  function Moving_Flag ( p : Standard_Natural_Vectors.Vector )
                       return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the pattern of the matrix of generators in the moving flag
  --   defined by the permutation p, using the (0,1,2) convention:
  --   2 can be any number while 0 and 1 are what they are.

  function Transformation ( n,p : integer32 )
                          return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the pattern of the n-by-n transformation matrix
  --   using the pivot p.  The pivot is the row of the descending or
  --   falling checker, obtained from the current or next permutation
  --   depending respectively whether specializing or generalizing.
  --   In the (0,1,2) convention '2' is the homotopy parameter t.

  -- REQUIRED : 1 <= p < n.

-- Part II: patterns of the k-plane, determined by white checkers

  procedure Sort_Rows ( r,c : in out Standard_Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Sorts the position of the white checkers along the row indices in r
  --   in increasing order and permutes the corresponding columns in c.
  --   This routine is auxiliary to building localization patterns.

  function Column_Pattern
              ( n,k : integer32; p,r,c : Standard_Natural_Vectors.Vector )
              return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the pattern for the n-by-k coordinate matrix for a k-plane
  --   in n-space whose shape is determined by black and white checkers.
  --   Generators of the k-plane are in the columns of a n-by-k matrix.

  -- ON ENTRY :
  --   n        dimension of the ambient space;
  --   k        dimension of the plane;
  --   p        permutation determines the position of the black checkers;
  --   r        rows of the white checkers;
  --   c        columns of the white checkers (start index at end of c).

  -- ON RETURN :
  --   An integer n-by-k matrix with entries 0, 1, or 2.
  --   The '2' in the matrix on return indicates a variable entry,
  --   entries with 0 or 1 are elements fixed to 0 or 1 respectively.

  function Row_Pattern
              ( n,k : integer32; p,r,c : Standard_Natural_Vectors.Vector )
              return Standard_Natural_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the Transpose of Column_Pattern(n,k,p,r,c).
  --   Storing the generators of a k-plane in the rows of an k-by-n matrix
  --   usually saves space when printing the pattern.

  function Degree_of_Freedom
              ( lp : Standard_Natural_Matrices.Matrix ) return natural32;

  -- DESCRIPTION :
  --   Returns the number of entries filled with 2 in the
  --   localization pattern lp.

  function Degree_of_Freedom
              ( p,r,c : Standard_Natural_Vectors.Vector ) return natural32;

  -- DESCRIPTION :
  --   Returns Degree_of_Freedom(Localization_Pattern(p,r,c)).

  function Rank ( lp : Standard_Natural_Matrices.Matrix;
                  i,j : integer32 ) return integer32;

  -- DESCRIPTION : 
  --   The rank of the variable x(i,j) is its order in the matrix,
  --   when flattened to a vector in row wise fashion.

  function Permute_Index
              ( p,q : Standard_Natural_Matrices.Matrix ) return integer32;

  -- DESCRIPTION :
  --   On input in p and q are two consecutive row patterns,
  --   respectively for the previous and current node.
  --   On return is 0 if no permutations are necessary,
  --   otherwise the index on return indicates the variable
  --   that should be permuted with the next one.

-- PART III : mapping patterns to complex matrices

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : Standard_Complex_Vectors.Vector )
               return Standard_Complex_Matrices.Matrix;
  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : DoblDobl_Complex_Vectors.Vector )
               return DoblDobl_Complex_Matrices.Matrix;
  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : QuadDobl_Complex_Vectors.Vector )
               return QuadDobl_Complex_Matrices.Matrix;

  -- DESCRIPTION :
  --   Maps the values in the vector x into a complex matrix,
  --   in standard double, double double, or quad double precision,
  --   using the localization pattern m.

  -- REQUIRED : x'length >= Degree_of_Freedom(m).

  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : Standard_Complex_Matrices.Matrix )
               return Standard_Complex_Vectors.Vector;
  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : DoblDobl_Complex_Matrices.Matrix )
               return DoblDobl_Complex_Vectors.Vector;
  function Map ( m : Standard_Natural_Matrices.Matrix;
                 x : QuadDobl_Complex_Matrices.Matrix )
               return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Maps the values in the matrix x into a complex vector,
  --   in standard double, double double, or quad double precision,
  --   using the localization pattern m.

end Checker_Localization_Patterns;
