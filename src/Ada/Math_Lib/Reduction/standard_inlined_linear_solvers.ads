with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;

package Standard_Inlined_Linear_Solvers is

-- DESCRIPTION :
--   An inlined solver of a linear system with complex coefficients
--   operates on columns of real and imaginary parts.
--   The complex arithmetic is inlined in the loops over the columns.
--   The design of the procedures is taken from the Fortran77 LINPACK.

  procedure lufac ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec; 
                    dim : in integer32;
                    ipvt : in out Standard_Integer_Vectors.Vector;
                    info : out integer32 );

  -- DESCRIPTION :
  --   Computes the LU factorization of the complex matrix stored as
  --   a pair of real and imaginary parts of the column coefficients.
  --   Applies partial row pivoting.

  -- REQUIRED : rcols'range = icols'range = ipvt'range = 1..dim,
  --   and for k in 1..dim: rcols(k)'range = icols(k)'range.

  -- ON ENTRY :
  --   rcols    real parts of the columns of a complex matrix;
  --   icols    imaginary parts of the columns of a complex matrix;
  --   dim      dimension of the matrix;
  --   ipvt     work space for the pivoting information.

  -- ON RETURN :
  --   rcols    real parts of the LU factorization of the matrix;
  --   icols    imaginary parts of the LU factorization of the matrix;
  --   ipvt     pivoting information;
  --   info     0 if all went well, otherwise info indicates the first
  --            column where a zero pivot occurred.

  procedure lusolve ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec;
                      icols : in Standard_Floating_VecVecs.Link_to_VecVec;
                      dim : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      rb : in Standard_Floating_Vectors.Link_to_Vector;
                      ib : in Standard_Floating_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Uses the output of lufac to solve a linear system.

  -- REQUIRED : rcols'range = icols'range = ipvt'range = 1..dim,
  --   and for k in 1..dim: rcols(k)'range = icols(k)'range;
  --   and rb'range = ib'range = 1..dim.

  -- ON ENTRY :
  --   rcols    real parts of the LU factorization of the matrix;
  --   icols    imaginary parts of the LU factorization of the matrix;
  --   dim      dimension of the matrix;
  --   ipvt     pivoting information, as computed by lufac;
  --   rb       real parts of the right hand side vector;
  --   ib       imaginary parts of the right hand side vector.

  -- ON RETURN :
  --   rb       real parts of the solution vector;
  --   ib       imaginary parts of the solution vector.

  function Norm1 ( rcols : Standard_Floating_VecVecs.Link_to_VecVec;
                   icols : Standard_Floating_VecVecs.Link_to_VecVec )
                 return double_float;

  -- DESCRIPTION :
  --   Returns the 1-norm of the complex matrix with real parts of 
  --   its columns in rcols and the imaginary parts in icols.

  -- REQUIRED : rcols'range = icols'range, and for all k in rcols'range:
  --   rcols(k)'range = icols(k)'range.

  procedure estco ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    dim : in integer32;
                    ipvt : in Standard_Integer_Vectors.Vector;
                    anorm : in double_float;
                    rz : in Standard_Floating_Vectors.Link_to_Vector;
                    iz : in Standard_Floating_Vectors.Link_to_Vector;
                    rcond : out double_float );

  -- DESCRIPTION :
  --   Estimates the condition number, given a LU factorization.

  -- REQUIRED : rcols'range = icols'range = ipvt'range = 1..dim,
  --   and for k in 1..dim: rcols(k)'range = icols(k)'range;
  --   and rz'range = iz'range = 1..dim.

  -- ON ENTRY :
  --   rcols    columns with the real parts of the outcome of lufac;
  --   icols    columns with the imaginary parts of the outcome of lufac;
  --   dim      dimension of the matrix;
  --   ipvt     pivoting information, as computed by lufac;
  --   anorm    the 1-norm of the matrix, computed by Norm1,
  --            on the original matrix before the lufac;
  --   rz       work space vector of 1..dim for real parts;
  --   iz       work space vector of 1..dim for imaginary parts.

  -- ON RETURN :
  --   rcond    estimate for the inverse of the condition number.

  procedure lufco ( rcols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    icols : in Standard_Floating_VecVecs.Link_to_VecVec;
                    dim : in integer32;
                    ipvt : in out Standard_Integer_Vectors.Vector;
                    rz : in Standard_Floating_Vectors.Link_to_Vector;
                    iz : in Standard_Floating_Vectors.Link_to_Vector;
                    rcond : out double_float );

  -- DESCRIPTION :
  --   Computes a LU factorization
  --   Estimates the condition number, given a LU factorization.

  -- REQUIRED : rcols'range = icols'range = ipvt'range = 1..dim,
  --   and for k in 1..dim: rcols(k)'range = icols(k)'range;
  --   and rz'range = iz'range = 1..dim.

  -- ON ENTRY :
  --   rcols    real parts of the columns of a complex matrix;
  --   icols    imaginary parts of the columns of a complex matrix;
  --   dim      dimension of the matrix;
  --   ipvt     work space for the pivoting information.
  --   anorm    the 1-norm of the matrix, computed by Norm1,
  --            on the original matrix before the lufac;
  --   rz       work space vector of 1..dim for real parts;
  --   iz       work space vector of 1..dim for imaginary parts.

  -- ON RETURN :
  --   rcols    columns with the real parts of the outcome of lufac;
  --   icols    columns with the imaginary parts of the outcome of lufac;
  --   ipvt     pivoting information, as computed by lufac;
  --   rcond    estimate for the inverse of the condition number.

end Standard_Inlined_Linear_Solvers;
