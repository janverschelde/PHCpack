with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;

package Standard_Inlined_Linearization is

-- DESCRIPTION :
--   Test the performance on solving linear systems of power series with
--   linearization, flat data structures, and inlined linear solvers.

  procedure Row_Matrix_Multiply
              ( rArows,iArows : in Standard_Floating_VecVecs.Link_to_VecVec;
                rx,ix,ry,iy : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Computes the matrix-vector product y = A*x,
  --   with A given as the real and imaginary parts of its rows,
  --   and x and y as vectors of real and imaginary parts.

  -- REQUIRED :
  --   ry'range = iy'range = rArows'range = iArows'range, and 
  --   rx'range = ix'range = rArows(k)'range = iArows(k)'range,
  --   for all k in rArows'range.

  -- ON ENTRY :
  --   rArows   real parts of the complex numbers on the rows of A;
  --   iArows   imaginary parts of the complex numbers on the rows of A;
  --   rx       real parts of the numbers of the complex vector x;
  --   ix       imaginary parts of the numbers of the complex vector x.
  --   ry       space allocated for the real parts of y;
  --   iy       space allocated for the imaginary parts of y.

  -- ON RETURN :
  --   ry       real parts of the complex vector y = A*x;
  --   iy       imaginary parts of the complex vector y = A*x.

  procedure Inlined_Solve_by_lufac
              ( dim : in integer32; 
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, without condition number estimate,
  --   where the matrices A(k) are given as vector 
  --   using an inlined linear system solver,
  --   with allocated real work space vectors.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   dim      dimension of the vectors;
  --   b        the right hand side coefficients of a vector series;
  --   rv       all real pars of all A(k), for k in 1..degree;
  --   iv       all imaginary pars of all A(k), for k in 1..degree;
  --   rc       real parts of the columns of A(0);
  --   ic       imaginary parts of the columns of A(0);
  --   rb       allocated work space vector for all real parts
  --            of the solution vectors;
  --   ib       allocated work space vector for all imaginary parts
  --            of the solution vectors;
  --   ry       allocated work space vector of range 1..dim;
  --   iy       allocated work space vector of range 1..dim.

  -- ON RETURN :
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   b        b contains the coefficients of the solution series x,
  --            provided info = 0, otherwise b is unchanged;
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k);
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

  procedure Inlined_Solve_by_lufac
              ( A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, without condition number estimate,
  --   using an inlined linear system solver.
  --   Allocates and deallocates all real work space vectors.

  -- REQUIRED :
  --   A'range = b'range = 0..deg, for deg >= 0.
  --   Moreover, all matrices in A are square.

  -- ON ENTRY :
  --   A        the coefficient matrices in the matrix series;
  --   b        the right hand side coefficients of a vector series.

  -- ON RETURN :
  --   A        A(0) contains the output of lufac on A(0);
  --   b        b contains the coefficients of the solution series x,
  --            provided info = 0, otherwise b is unchanged;
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

end Standard_Inlined_Linearization;
