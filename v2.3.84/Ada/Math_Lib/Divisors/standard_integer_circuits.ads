with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;

package Standard_Integer_Circuits is

-- DESCRIPTION :
--   A circuit in a point configuration is a smallest affinely dependent
--   subset of points.  Given a point configuration with coordinates of
--   the points in the columns of a matrix, this package offers tools to
--   to compute the affine dependencies.

  function Submatrix ( A : Matrix; columns : Vector ) return Matrix;

  -- DESCRIPTION :
  --   Returns the submatrix of A defined by selecting the columns with
  --   indices in the vector columns.

  function Rank ( A : Matrix; columns : Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the rank of the submatrix of A selected by the indices
  --   in the columns vector.  If the rank of the submatrix of A, as defined
  --   by the columns is full, then the submatrix contains a basis for the
  --   circuit computation.

  function Circuit ( file : file_type;
                     A : Matrix; basis : Vector; k : integer32 ) return Vector;
  function Circuit ( A : Matrix; basis : Vector; k : integer32 ) return Vector;

  -- DESCRIPTION :
  --   Returns the coefficients of the affine dependencies that span the
  --   circuit defined by the points indexed by the basis and k.
  --   Writes extra output to file when file is provided.
  --   The vector on return has range basis'first..basis'last+1.
  --   The coefficients of the vector on return are coefficients with
  --   the points with corresponding indices in the basis.
  --   The last coefficient of the vector on return is the coefficient
  --   of point k.

  -- REQUIRED : rank(A,basis) = basis'last and k is not in basis. 

  function Kernel_Vector ( n,k : integer32; b,c : Vector ) return Vector;

  -- DESCRIPTION :
  --   Given a circuit defined by the basis b, index k, and coefficients in c,
  --   the corresponding n vector in the kernel of the matrix A is returned.

  function Is_In ( b : Vector; k : integer32 ) return boolean;

  -- DESCRIPTION :
  --   Returns true if b(i) = k for some i in b'range and false otherwise.

  function Is_Zero ( v : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if every element in v is zero, returns false otherwise.

  generic
    with procedure Report ( c : in Vector; k : integer32;
                            continue : out boolean );
  procedure Enumerate_All_Circuits ( A : in Matrix; basis : in Vector );

  -- DESCRIPTION :
  --   Enumerates all affine dependencies of the points in the columns
  --   of the matrix A that contain the given basis and some k.
  --   Each time a new circuit is found, the procedure Report is called.
  --   The enumeration stops when continue is set to false.

  generic
    with procedure Report ( b : in Vector; continue : out boolean );
  procedure Enumerate_All_Bases ( A : in Matrix; d : in integer32 ); 

  -- DESCRIPTION :
  --   Enumerates all choices of d columns in A with rank equal to d.
  --   After each choice of columns with rank d, the procedure Report
  --   is called.  The enumeration stops when continue turns false.

  generic
    with procedure Report ( b : in Vector; continue : out boolean );
  procedure Enumerate_Circuits ( A : in Matrix; d : in integer32 ); 

  -- DESCRIPTION :
  --   Enumerates all candidate circuits of dimension d, defined by d+1
  --   columns of A.  While the first d columns span a basis, it could
  --   still be that all coefficients of the circuit are zero.
  --   After each new circuit is computed, Report is called.
  --   The enumeration stops when continue turns false.

  function Cardinality ( n,d : integer32 ) return integer32;

  -- DESCRIPTION :
  --   The number of circuits of d+1 points of a point configuration
  --   of n points is of the order (n choose d) times (n-d).
  --   This function returns this number and is used in the matrix
  --   representation of all kernel vectors.

  function Circuits ( A : Matrix; d : integer32 ) return List;

  -- DESCRIPTION :
  --   Returns a list of all kernel vectors of d+1 circuits.

  function List_to_Matrix ( L : List ) return Matrix;

  -- DESCRIPTION :
  --   Returns matrix with in its columns the elements of the list L.

  -- REQUIRED : Length_Of(L) > 0.

  function Circuits ( A : Matrix; d : integer32 ) return Matrix;

  -- DESCRIPTION :
  --   Returns all the circuits of d+1 points of A in a matrix of
  --   kernel vectors in the columns.  If B = Circuits(A,d), then A*B = 0.

end Standard_Integer_Circuits;
