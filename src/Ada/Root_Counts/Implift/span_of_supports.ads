with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;

package Span_of_Supports is

-- DESCRIPTION :
--   If the span of the support sets of a polynomial does not equal
--   the dimension of the ambient space, then the mixed volume will
--   equal zero --- if it is applicable and the system is square ---
--   and in all cases, a coordinate transformation needs to be applied
--   for an efficient polyhedral solver.  This package provides tools
--   to compute the span and in case of not full rank, a basis of
--   normal vectors for use in a coordinate transformation.

  procedure Triangulate_Span
              ( s : in Lists_of_Integer_Vectors.List;
                n : in integer32; rank : out natural32;
                ipvt : out Standard_Integer_Vectors.Vector;
                span : out Standard_Integer_Matrices.Matrix );
  procedure Triangulate_Span
              ( s : in Lists_of_Integer_Vectors.List;
                n : in integer32; rank : out natural32;
                ipvt : out Standard_Integer_Vectors.Vector;
                span : out Standard_Integer64_Matrices.Matrix );

  -- DESCRIPTION :
  --   Computes the rank of the point configuration given in s
  --   and triangulates the vector basis of its span.
  --   If span is a 64-bit integer matrix, then 64-bit arithmetic
  --   will be applied to triangulate the matrix.

  -- ON ENTRY :
  --   s        list of n-dimensional points;
  --   n        ambient dimension, all points in s are of range 1..n.

  -- ON RETURN :
  --   rank     rank of the span of the point configuration;
  --   ipvt     records coordinate swaps for zeroes;
  --   span     upper triangular matrix with nonzero elements on
  --            the diagonal for rows 1 to rank.

  function Normal_Vectors
              ( U : in Standard_Integer_Matrices.Matrix;
                n,r : in integer32 )
              return Standard_Integer_Matrices.Matrix;
  function Normal_Vectors
              ( U : in Standard_Integer64_Matrices.Matrix;
                n,r : in integer32 )
              return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Given an upper triangular matrix of rank r in n-space,
  --   returns a basis of vectors in the kernel of the matrix U.
  --   If U is a 64-bit matrix, then 64-bit arithmetic is used.

  -- ON ENTRY :
  --   U        an upper triangular matrix, with nonzero elements
  --            on the diagonal for all rows from 1 to r;
  --   n        number of rows and columns in the matrix U;
  --   r        rank of the matrix U.

  -- ON RETURN :
  --   a  matrix n-r linearly independent columns so that U*V = 0.

  function Apply_Pivots
              ( A : Standard_Integer_Matrices.Matrix;
                ipvt : Standard_Integer_Vectors.Vector )
              return Standard_Integer_Matrices.Matrix;
  function Apply_Pivots
              ( A : Standard_Integer64_Matrices.Matrix;
                ipvt : Standard_Integer_Vectors.Vector )
              return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Applies the pivoting information in ipvt to the matrix in A.

  -- ON ENTRY :
  --   A        a matrix with vectors in its columns;
  --   ipvt     pivoting information on the rows of A,
  --            that is: A'range(1) = ipvt'range.

  procedure Rank_of_Support
              ( s : in Lists_of_Integer_Vectors.List; n : in integer32;
                debug : in boolean; rank : out natural32;
                normals : out Standard_Integer_Matrices.Link_to_Matrix );
  procedure Rank_of_Support
              ( s : in Lists_of_Integer_Vectors.List; n : in integer32;
                debug : in boolean; rank : out natural32;
                normals : out Standard_Integer64_Matrices.Link_to_Matrix );

  -- DESCRIPTION :
  --   Computes the rank of a point configuration
  --   and also the normals perpendicular to it span if not full rank.
  --   If normals is a 64-bit matrix, then 64-bit arithmetic is used.

  -- ON ENTRY :
  --   s        a list of points in n-space;
  --   n        ambient dimension, all points in s have range 1..n;
  --   debug    if true, then additional output is written.

  -- ON RETURN :
  --   rank     the rank of the span of the point configuration in s;
  --   normals  if rank < n, then normals will have n - rank columns
  --            containing the coordinates of the vectors that are
  --            perpendicular to the span.

  function Rank32_of_Support
              ( s : Lists_of_Integer_Vectors.List; n : integer32;
                debug : boolean ) return natural32;
  function Rank64_of_Support
              ( s : Lists_of_Integer_Vectors.List; n : integer32;
                debug : boolean ) return natural32;

  -- DESCRIPTION :
  --   Returns the rank of the n-vectors obtained by differences of
  --   all points with the first point in s.
  --   If debug is true, then output will be written to screen.
  --   The Rank32_ uses 32-bit and Rank64_ works with 64-bit arithmetic.
 
  -- REQUIRED : Length_Of(s) >= 2.

  function Cayley_Embedding
              ( point : Standard_Integer_Vectors.Vector; k,r : integer32 )
              return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the Cayley embedding of the given point,
  --   extended with the coordinates of the k-th point
  --   of the r-dimensional standard simplex.

  function Cayley_Embedding
              ( supports : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
              return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Given a sequence of r distinct support sets,
  --   this function adds to every point of the k-th support 
  --   the coordinates corresponding the k-th corner of the 
  --   standard (r-1)-dimensional simplex.
 
  -- REQUIRED : supports'range = 1..r for some r >= 1.

  function Remove_Cayley_Embedding
             ( transfo : Standard_Integer_Matrices.Matrix;
               dim : in integer32 )
             return Standard_Integer_Matrices.Matrix;

  -- DECRIPTION :
  --   Returns the first dim rows and columns of the matrix transfo.

  function Remove_Cayley_Rows
             ( mat : Standard_Integer_Matrices.Matrix;
               dim : in integer32 )
             return Standard_Integer_Matrices.Matrix;

  -- DECRIPTION :
  --   Returns the first dim rows and all columns of the matrix mat.

-- FOR TESTING PURPOSES :

  function Support_Function
              ( s : Lists_of_Integer_Vectors.List;
                normals : Standard_Integer_Matrices.Matrix )
              return Standard_Integer_Matrices.Matrix;
  function Support_Function
              ( s : Lists_of_Integer_Vectors.List;
                normals : Standard_Integer64_Matrices.Matrix )
              return Standard_Integer64_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns the matrix of values of the support function for the points
  --   in the list s with the normal vectors in the columns of normals.
  --   The matrix on return has as many rows as there are points in s
  --   and as many columns as there are normals.

  function Random_Lower
              ( dim,low,upp : integer32 )
              return Standard_Integer_Matrices.Matrix;

  -- DESCRIPTION :
  --   Returns a random square matrix of dimension dim with
  --   ones on the diagonal and random numbers in the range low..upp
  --   for elements below the diagonal.

  function Random_Support
              ( dim,nbp,low,upp,rnk : integer32 )
              return Lists_of_Integer_Vectors.List;

  -- DESCRIPTION :
  --   Returns a list of points, randomly generated and with an upper
  --   bound on the rank of the configuration.

  -- REQUIRED : 0 < rnk <= dim

  -- ON ENTRY :
  --   dim      ambient dimension, range of generated points is 1..dim;
  --   nbp      number of points in the list on return;
  --   low      lower bound on the numbers in the random number generator;
  --   upp      upper bound on the numbers in the random number generator;
  --   rnk      bound on the rank, the rank could be lower if low and upp
  --            are very close to each other and/or nbp is small.

  function Random_Tuple_of_Supports
              ( dim,mix,low,upp,rnk : integer32;
                nbp : Standard_Integer_Vectors.Vector )
              return Arrays_of_Integer_Vector_Lists.Array_of_Lists;

  -- DESCRIPTION :
  --   Generates a tuple of point configuration.

  -- ON ENTRY :
  --   dim      ambient dimension;
  --   mix      number of distinct supports;
  --   low      lower bound on the coordinates in the supports;
  --   upp      upper bound on the coordinates in the supports;
  --   rnk      bound on the rank of the point configuration;
  --   nbp      vector of range 1..mix,
  --            nbp(k) equals the number of points in the k-th support.

  -- ON RETURN :
  --   A tuple of supports of range 1..mix of the given characteristics.

end Span_of_Supports;
