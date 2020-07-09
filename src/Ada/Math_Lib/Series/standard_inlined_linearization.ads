with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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

  procedure Inlined_Solve_Head_by_lufac
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rb,ib : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32 );

  -- DESCRIPTION :
  --   Solves the leading term of a linear series system,
  --   computing the LU factorization of the matrix with the real parts
  --   of its columns stored in rc and its imaginary parts in ic.

  -- REQUIRED :
  --   rc'range = ic'range = 1..dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..dim.
  --   rb'range = ib'range = 1..dim and ipvt'range = 1..dim.

  -- ON ENTRY :
  --   dim      dimension of the vectors;
  --   rc       real parts of the columns of A(0);
  --   ic       imaginary parts of the columns of A(0);
  --   rb       real parts of the leading right hand coefficients
  --            of the linear series system;
  --   ib       imaginary parts of the leading right hand coefficients
  --            of the linear series system.

  -- ON RETURN :
  --   rc       real parts of the output of the inlined lufac;
  --   ic       imaginary parts of the output of the inlined lufac;
  --   rb       rb stores the real parts of the leading coefficient
  --            vector of the solution series, if info = 0;
  --   ib       ib stores the imaginary parts of the leading coefficient
  --            vector of the solution series, if info = 0;
  --   ipvt     pivoting information on the LU factorization;
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

  procedure Inlined_Solve_Head_by_lufco
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rb,ib : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the leading term of a linear series system,
  --   computing the LU factorization of the matrix with the real parts
  --   of its columns stored in rc and its imaginary parts in ic.
  --   Estimates the inverse of the condition number.

  -- REQUIRED :
  --   rc'range = ic'range = 1..dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..dim.
  --   rb'range = ib'range = 1..dim and ipvt'range = 1..dim.

  -- ON ENTRY :
  --   dim      dimension of the vectors;
  --   rc       real parts of the columns of A(0);
  --   ic       imaginary parts of the columns of A(0);
  --   rb       real parts of the leading right hand coefficients
  --            of the linear series system;
  --   ib       imaginary parts of the leading right hand coefficients
  --            of the linear series system;
  --   ry       allocated work space vector of range 1..dim;
  --   iy       allocated work space vector of range 1..dim.

  -- ON RETURN :
  --   rc       real parts of the output of the inlined lufac;
  --   ic       imaginary parts of the output of the inlined lufac;
  --   rb       rb stores the real parts of the leading coefficient
  --            vector of the solution series, if rcond + 1.0 /= 1.0;
  --   ib       ib stores the imaginary parts of the leading coefficient
  --            vector of the solution series, if rcond + 1.0 /= 1.0;
  --   ipvt     pivoting information on the LU factorization;
  --   rcond    returned by lufco, estimate for the inverse condition
  --            number of the lead coefficient matrix,
  --            if zero, then deemed singular by lufco.

  procedure Inlined_Solve_Tail_by_lusolve
              ( deg,dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : in Standard_Integer_Vectors.Vector;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector;
                idx : in integer32 := 1 );

  -- DESCRIPTION :
  --   Computes the tail of the solution of a linear series system,
  --   starting at index idx and ending at the degree deg.

  -- REQUIRED : the info of the lufac was zero, and
  --   rc'range = ic'range = 1..dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..dim.
  --   rv'range = iv'range = 1..deg and for all k in rv'range:
  --   rv(k)'range = iv(k)'range = 1..dim and for all i in rv(k)'range:
  --   rv(k)(i)'range = iv(k)(i)'range = 1..dim.
  --   rb'range = ib'range = 0..deg and for all k in rb'range:
  --   rb(k)'range = ib(k)'range = 1..dim.

  -- ON ENTRY :
  --   deg      degree of the series, if omitted, then rb'last = deg,
  --            if provided, then only coefficients up to the index deg,
  --            with deg included will be considered, deg <= rb'last;
  --   dim      dimension of the vectors;
  --   rc       real parts of the columns of the head matrix;
  --   ic       imaginary parts of the columns of the head matrix;
  --   rv       all real parts of all A(k), for k in 1..degree;
  --   iv       all imaginary parts of all A(k), for k in 1..degree;
  --   rb       rb(k) stores the coefficient vectors of the solution,
  --            for all k < idx and for k in idx..deg, rb(k) stores
  --            the real parts of the right hand coefficients;
  --   ib       ib(k) stores the coefficient vectors of the solution,
  --            for all k < idx and for k in idx..deg, ib(k) stores
  --            the real parts of the right hand coefficients;
  --   ipvt     pivoting information computed by lufac;
  --   ry       allocated work space vector of range 1..dim;
  --   iy       allocated work space vector of range 1..dim;
  --   idx      start index at the 

  -- ON RETURN :
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   rb       rb(k) stores the real parts of the solution series,
  --            for k in range idx..deg;
  --   ib       ib(k) stores the imaginary parts of the solution series,
  --            for k in range idx..deg;
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

  procedure Inlined_Solve_by_lufac
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Inlined_Solve_by_lufac
              ( deg,dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                info : out integer32;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the linear system A(t)*x(t) = b(t) of power series,
  --   using LU factorization on the leading coefficient matrix of A(t),
  --   without condition number estimate.
  --   With linearization, A(t) is a power series with matrix coefficients
  --   and b(t) is a power series with vector coefficients.
  --   The matrix coefficients are represented as vectors of floating-point
  --   vectors, for a better performing inlined linear system solver.
  --   The leading matrix coefficient must be given in column format,
  --   while the other matrix coefficients must be given in row format.

  -- REQUIRED :
  --   rc'range = ic'range = 1..dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..dim.
  --   rv'range = iv'range = 1..deg and for all k in rv'range:
  --   rv(k)'range = iv(k)'range = 1..dim and for all i in rv(k)'range:
  --   rv(k)(i)'range = iv(k)(i)'range = 1..dim.
  --   rb'range = ib'range = 0..deg and for all k in rb'range:
  --   rb(k)'range = ib(k)'range = 1..dim.

  -- ON ENTRY :
  --   deg      degree of the series, if omitted, then rb'last = deg,
  --            if provided, then only coefficients up to the index deg,
  --            with deg included will be considered, deg <= rb'last;
  --   dim      dimension of the vectors;
  --   rc       real parts of the columns of A(0);
  --   ic       imaginary parts of the columns of A(0);
  --   rv       all real parts of all A(k), for k in 1..degree;
  --   iv       all imaginary parts of all A(k), for k in 1..degree;
  --   rb       real parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ib       imaginary parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ry       allocated work space vector of range 1..dim;
  --   iy       allocated work space vector of range 1..dim.

  -- ON RETURN :
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k);
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   info     returned by lufac, if nonzero, then the lead coefficient
  --            matrix was deemed singular by lufac.

  procedure Inlined_Solve_by_lufco
              ( dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector );
  procedure Inlined_Solve_by_lufco
              ( deg,dim : in integer32; 
                rc,ic : in Standard_Floating_VecVecs.Link_to_VecVec;
                rv,iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                rb,ib : in Standard_Floating_VecVecs.Link_to_VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float;
                ry,iy : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Solves the linear system A(t)*x(t) = b(t) of power series,
  --   using LU factorization on the leading coefficient matrix of A(t),
  --   with a condition number estimate.
  --   With linearization, A(t) is a power series with matrix coefficients
  --   and b(t) is a power series with vector coefficients.
  --   The matrix coefficients are represented as vectors of floating-point
  --   vectors, for a better performing inlined linear system solver.
  --   The leading matrix coefficient must be given in column format,
  --   while the other matrix coefficients must be given in row format.

  -- REQUIRED :
  --   rc'range = ic'range = 1..dim and for all k in rc'range:
  --   rc(k)'range = ic(k)'range = 1..dim.
  --   rv'range = iv'range = 1..deg and for all k in rv'range:
  --   rv(k)'range = iv(k)'range = 1..dim and for all i in rv(k)'range:
  --   rv(k)(i)'range = iv(k)(i)'range = 1..dim.
  --   rb'range = ib'range = 0..deg and for all k in rb'range:
  --   rb(k)'range = ib(k)'range = 1..dim.

  -- ON ENTRY :
  --   deg      degree of the series, if omitted, then rb'last = deg,
  --            if provided, then only coefficients up to the index deg,
  --            with deg included will be considered, deg <= rb'last;
  --   dim      dimension of the vectors;
  --   rc       real parts of the columns of A(0);
  --   ic       imaginary parts of the columns of A(0);
  --   rv       all real parts of all A(k), for k in 1..degree;
  --   iv       all imaginary parts of all A(k), for k in 1..degree;
  --   rb       real parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ib       imaginary parts of the right hand coefficients of b(t),
  --            where b(t) is a series with vector coefficients;
  --   ry       allocated work space vector of range 1..dim;
  --   iy       allocated work space vector of range 1..dim.

  -- ON RETURN :
  --   rc       real parts of the output of lufac on A(0);
  --   ic       imaginary parts of the output of lufac on A(0);
  --   rb       rb(k) stores the real parts of the solution x(k);
  --   ib       ib(k) stores the imaginary parts of the solution x(k);
  --   ipvt     pivoting information on the LU factorization of A(0);
  --   rcond    returned by lufco, estimate for the inverse condition
  --            number of the lead coefficient matrix,
  --            if zero, then deemed singular by lufco.

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
  --   This wrapper procedure is for testing purposes,
  --   as allocation and deallocation does not make it thread safe.

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

  procedure Inlined_Solve_by_lufco
              ( A : in Standard_Complex_VecMats.VecMat;
                b : in Standard_Complex_VecVecs.VecVec;
                ipvt : out Standard_Integer_Vectors.Vector;
                rcond : out double_float );

  -- DESCRIPTION :
  --   Solves the linear system A*x = b, using LU factorization on the
  --   leading coefficient matrix of A, with condition number estimate,
  --   using an inlined linear system solver.
  --   Allocates and deallocates all real work space vectors.
  --   This wrapper procedure is for testing purposes,
  --   as allocation and deallocation does not make it thread safe.

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
  --   rcond    returned by lufco, estimate for the inverse condition
  --            number of the lead coefficient matrix,
  --            if zero, then deemed singular by lufco.

end Standard_Inlined_Linearization;
