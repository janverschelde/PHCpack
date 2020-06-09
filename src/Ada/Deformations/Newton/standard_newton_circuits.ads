with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Coefficient_Circuits;

package Standard_Newton_Circuits is

-- DESCRIPTION :
--   Newton's method is apply to approximate roots of polynomial systems,
--   represented by coefficient circuits in standard double precision.

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                info : out integer32; res,err : out double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Does one step with Newton's method on s, starting at v,
  --   using LU factorization to solve the linear system.
  --   No estimate for the condition number is computed.

  -- REQUIRED : s.neq = s.dim, i.e.: the system is square,
  --   and xr'range = xi'range = v'range = ipvt'range.

  -- ON ENTRY :
  --   s        coefficient system of circuits;
  --   v        vector with approximate values for a solution to s;
  --   xr       work space allocated for the real parts of v;
  --   xi       work space allocated for the imaginary parts of v;
  --   ipvt     pivoting vector for the LU factorization;
  --   verbose  flag for extra output to screen.

  -- ON RETURN :
  --   s        s.fx contains the update value to v, if info = 0,
  --            otherwise s.fx contains the function value of s at v,
  --            s.jm contains the LU factorization of the Jacobian
  --            at the given v;
  --   xr       real parts of the complex numbers in v;
  --   xi       imaginary parts of the complex numbers in v;
  --   ipvt     pivoting information of the LU factorization;
  --   info     if nonzero, then the Jacobian may be singular;
  --   res      residual, max norm of the function value at v;
  --   err      max norm of the update vector to v.

  procedure LU_Newton_Step
              ( s : in Standard_Coefficient_Circuits.Link_to_System;
                v : in out Standard_Complex_Vectors.Vector;
                xr,xi : in Standard_Floating_Vectors.Link_to_Vector;
                ipvt : in out Standard_Integer_Vectors.Vector;
                res,rco,err : out double_float;
                verbose : in boolean := true );

  -- DESCRIPTION :
  --   Does one step with Newton's method on s, starting at v,
  --   using LU factorization to solve the linear system.
  --   Provides an estimate for the inverse of the condition number 
  --   of the Jacobian matrix at v.

  -- REQUIRED : s.neq = s.dim, i.e.: the system is square,
  --   and xr'range = xi'range = v'range = ipvt'range.

  -- ON ENTRY :
  --   s        coefficient system of circuits;
  --   v        vector with approximate values for a solution to s;
  --   xr       work space allocated for the real parts of v;
  --   xi       work space allocated for the imaginary parts of v;
  --   ipvt     pivoting vector for the LU factorization;
  --   verbose  flag for extra output to screen.

  -- ON RETURN :
  --   s        s.fx contains the update value to v, if info = 0,
  --            otherwise s.fx contains the function value of s at v,
  --            s.jm contains the LU factorization of the Jacobian
  --            at the given v;
  --   xr       real parts of the complex numbers in v;
  --   xi       imaginary parts of the complex numbers in v;
  --   ipvt     pivoting information of the LU factorization;
  --   rco      estimate for the inverse of the condition number,
  --            if zero, then the Jacobian matrix may be singular;
  --   res      residual, max norm of the function value at v;
  --   err      max norm of the update vector to v.

end Standard_Newton_Circuits;
