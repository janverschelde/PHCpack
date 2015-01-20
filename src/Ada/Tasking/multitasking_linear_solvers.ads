with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;

package multitasking_linear_solvers is

-- DESCRIPTION :
--   This package provides multitasked solvers of linear systems,
--   based on Ada translations of linpack's lufac and lusolve.

-- LU FACTORIZATION :

  procedure silent_lufac
              ( t : in integer32;
                a : in out Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure reporting_lufac
              ( t : in integer32; 
                lu : in Standard_Complex_Matrices.Matrix;
                a : in out Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure silent_lufac
              ( t : in integer32;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure reporting_lufac
              ( t : in integer32;
                a : in out DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure silent_lufac
              ( t : in integer32;
                a : in out QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );
  procedure reporting_lufac
              ( t : in integer32;
                a : in out QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Computes an LU factorization with t tasks.  The reporting version
  --   writes output to allow to monitor the progress of the computations.

  -- REQUIRED : a'range(1) = a'range(2) = ipvt'range = 1..n.

  -- ON ENTRY :
  --   t        number of tasks;
  --   lu       for debugging purposes in the reporting versions,
  --            lu is the result of the LU decomposition;
  --   a        a complex n-by-n matrix for factorization;
  --   n        dimension of the matrix a.

  -- ON RETURN :
  --   a        contains L and U so a = L*U where L is lower
  --            and upper triangular matrix (modulo permutation);
  --   ipvt     pivot indices define the permutation matrix;
  --   info     0 is normal value,
  --            k if u(k,k) = 0.0 for singular matrices.

-- FORWARD and BACK SUBSTITUTION :

  procedure silent_lusolve
              ( t : in integer32;
                a : in Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out Standard_Complex_Vectors.Vector );
  procedure reporting_lusolve
              ( t : in integer32;
                a : in Standard_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out Standard_Complex_Vectors.Vector );
  procedure silent_lusolve
              ( t : in integer32;
                a : in DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out DoblDobl_Complex_Vectors.Vector );
  procedure reporting_lusolve
              ( t : in integer32;
                a : in DoblDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out DoblDobl_Complex_Vectors.Vector );
  procedure silent_lusolve
              ( t : in integer32;
                a : in QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out QuadDobl_Complex_Vectors.Vector );
  procedure reporting_lusolve
              ( t : in integer32;
                a : in QuadDobl_Complex_Matrices.Matrix;
                n : in integer32;
                ipvt : in Standard_Integer_Vectors.Vector;
                b : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Solves two triangular linear systems by t tasks.  The reporting
  --   version writes intermediate output for monitoring the computations.

  -- REQUIRED : a'range(1) = a'range(2) = b'range = 1..n.

  -- ON ENTRY :
  --   t        number of tasks;
  --   a        contains the factorization computed by lufac or lufco,
  --            a = L*U, modulo the permutation matrix defined by ipvt;
  --   n        dimension of the matrix, number of rows and columns;
  --   ipvt     pivoting information of lufac or lufco;
  --   b        right hand side vector the linear system.

  -- ON RETURN :
  --   b        solution to the linear system.

end multitasking_linear_solvers;
