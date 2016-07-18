with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Dense_Vector_Series;
with DoblDobl_Dense_Matrix_Series;

package DoblDobl_Matrix_Series_Solvers is

-- DESCRIPTION :
--   Given a linear system of truncated power series, linearization with
--   matrix series and vector series solves the linear system.
--   The coefficients and computations are in double double precision.

  procedure Solve_Lead_by_lufac
              ( A : in DoblDobl_Dense_Matrix_Series.Matrix;
                b : in DoblDobl_Dense_Vector_Series.Vector;
                a0lu : out DoblDobl_Complex_Matrices.Matrix;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                x : out DoblDobl_Dense_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies LU factorization and back substitution to compute the
  --   constant coefficient of the solution series to A*x = b.

  -- REQUIRED :
  --   A.deg >= 0 and b.deg >= 0.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series.

  -- ON RETURN :
  --   a0lu     LU factorization of A.cff(0);
  --   ipvt     pivoting information on the LU factorization of A.cff(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac;
  --   x        x.cff(0) is the constant coefficient of the solution
  --            and x.deg = 0, provided info = 0,
  --            if info /= 0, then x.deg = -1 and x.cff(0) is undefined.

  procedure Solve_Next_by_lusolve
              ( A : in DoblDobl_Dense_Matrix_Series.Matrix;
                b : in DoblDobl_Dense_Vector_Series.Vector;
                a0lu : in DoblDobl_Complex_Matrices.Matrix;
                ipvt : in Standard_Integer_Vectors.Vector;
                x : in out DoblDobl_Dense_Vector_Series.Vector );

  -- DESCRIPTION :
  --   Applies LU back substitution, calling lusolve,
  --   on the linear system with a modified right hand side vector,
  --   using previously computed coefficient series vector of x.

  -- REQUIRED :
  --   All coefficients up to x.deg are defined.

  -- ON ENTRY :
  --   A        the coefficient matrix as a matrix series;
  --   b        the right hand side as a vector series;
  --   a0lu     output of Solve_Lead_by_lufac, contains the LU
  --            factorization of the leading coefficient of A;
  --   ipvt     output of Solve_Lead_by_lufac, contains the pivoting
  --            information in the LU factorization of A.cff(0);
  --   x        previously computed coefficients of the solution,
  --            at the very least, x.cff(0) must be defined.

  -- ON RETURN :
  --   x        computed coefficients up to b.deg.

end DoblDobl_Matrix_Series_Solvers;
