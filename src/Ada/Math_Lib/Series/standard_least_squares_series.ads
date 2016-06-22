with Standard_Integer_Vectors;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_Matrices;

package Standard_Least_Squares_Series is

-- DESCRIPTION :
--   With QR decomposition on an overdetermined linear system
--   leads to a solution of the system in the least squares sense.

  procedure QRD ( x : in out Standard_Dense_Series_Matrices.Matrix;
                  qraux : in out Standard_Dense_Series_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean );

  -- DESCRIPTION :
  --   Uses Householder transformations to compute the QR decomposition of x.
  --   Column pivoting based on 2-norms of reduced columns is optional.

  -- REQUIRED : jpvt'range = x'range(2) = qraux'range,
  --            x'length(1) >= x'length(2)

  -- ON ENTRY :
  --   x          matrix whose decomposition is to be computed;
  --   jpvt       controls the selection of the pivot columns.
  --              The k-th column x(k) of x is placed in one of three classes
  --              according to the value of jpvt(k):
  --                if jpvt(k) > 0, then x(k) is an initial column.
  --                if jpvt(k) = 0, then x(k) is a free column.
  --                if jpvt(k) < 0, then x(k) is a final column.
  --              before the decomposition is computed, initial columns are
  --              moved to the beginning of the array x and final columns to
  --              the end.  Both initial and final columns are frozen in place
  --              during the computation and only free columns are moved.
  --              At the k-th stage of the reduction, if x(k) is occupied by a
  --              free column it is interchanged with the free column of
  --              largest reduced norm.  jpvt is not referenced if not piv.
  --   piv        column pivoting is performed if and only if piv is true.

  -- ON RETURN :
  --   x          x contains in its upper triangle the upper triangular matrix
  --              R of the QR factorization.  Below its diagonal x contains
  --              information from which the orthogonal part of the
  --              decomposition can be recovered.  Note that if pivoting has
  --              been requested, the decomposition is not that of the original
  --              matrix x but that of x with its columns permuted as described
  --              by jpvt.
  --   qraux      qraux contains further information required to recover
  --              the orthogonal part of the decomposition.
  --   jpvt       jpvt(k) contains the index of the column of the original
  --              matrix that has been interchanged into the k-th column,
  --              if pivoting was requested.
  
  -- ACKNOWLEDGMENT :
  --   This Ada Version is a translation of the LINPACK version, dated 08/14/78,
  --   written by G.W. Stewart, University of Maryland, Argonne National Lab.

end Standard_Least_Squares_Series;
