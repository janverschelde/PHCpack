with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package PHCpack is

-- DESCRIPTION :
--   PHC is a numerical general-purpose solver for polynomial systems by
--   homotopy continuation.  This package collects all major features of
--   the package as organized according to the four stages of the solver.

-- 1. PRE-PROCESSING : SCALING AND REDUCTION

  procedure Equation_Scaling
                 ( file : in file_type; p : in Poly_Sys; s : out Poly_Sys );

  -- DESCRIPTION : equation scaling by dividing average coefficient.

  -- ON ENTRY :
  --   file        to write results on, must be opened for output;
  --   p           a polynomial system.

  -- ON RETURN :
  --   s           a scaled polynomial system.

  procedure Linear_Reduction
                 ( file : in file_type; p : in Poly_Sys; r : out Poly_Sys );

  -- DESCRIPTION : linear reduction of the coefficient matrix of p.

  -- ON ENTRY :
  --   file        to write results on, must be opened for output;
  --   p           a polynomial system.

  -- ON RETURN :
  --   r           a polynomial system with reduced coefficient matrix.

-- 2. ROOT COUNTING AND CONSTRUCTION OF START SYSTEM

  procedure Total_Degree
                 ( file : in file_type; p : in Poly_Sys; d : out natural32 );

  -- DESCRIPTION : computation of the total degree.

  -- ON ENTRY :
  --   file        to write results on, must be opened for output;
  --   p           a polynomial system.

  -- ON RETURN :
  --   d           total degree of the system p.

  procedure Total_Degree
                 ( file : in file_type; p : in Poly_Sys; d : out natural32;
                   q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION : construction of start system based on total degree.

  -- ON ENTRY :
  --   file        to write results on, must be opened for output;
  --   p           a polynomial system.

  -- ON RETURN :
  --   d           total degree of the system p.
  --   q           start system with same total degree as p;
  --   qsols       solutions of q.

  function Set_Structure_Bound ( p : Poly_Sys ) return natural32;

  -- DESCRIPTION :
  --   Returns the generalized Bezout bound, based on set structures.

  procedure Implicit_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural32 );

  -- DESCRIPTION : computation of mixed volume by implicit lifting.

  procedure Implicit_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural32;
                   q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION : construction of start system by implicit polyhedral homotopy.

  procedure Static_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural32 );

  -- DESCRIPTION : computation of mixed volume by static lifting.

  procedure Static_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural32;
                   q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION : construction of start system by static polyhedral homotopy.

  procedure Count_Roots
                 ( p : in out Poly_Sys; rc : out natural32;
                   q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Blackbox root counter, returns a start system q corresponding to the
  --   smallest root count rc.  The original equations in p may have been
  --   permuted to take advantage of equal supports.

-- 3. POLYNOMIAL CONTINUATION

  procedure Artificial_Parameter_Continuation
                 ( p,q : in Poly_Sys; sols : in out Solution_List;
                   k : in natural32 := 2;
                   a : in Complex_Number := Create(1.0);
                   target : in Complex_Number := Create(1.0) );
  procedure Artificial_Parameter_Continuation
                 ( file : in file_type; p,q : in Poly_Sys;
                   sols : in out Solution_List;
                   k : in natural32 := 2;
                   a : in Complex_Number := Create(1.0);
                   target : in Complex_Number := Create(1.0) );

  -- DESCRIPTION : continuation with the artificial-parameter homotopy
  --   h(x,t) = a*(1-t)^k*q(x) + t^k*p(x) = 0, for t going to the target.

  -- ON ENTRY :
  --   file        to write results on, must be opened for output;
  --   p           target system;
  --   q           start system;
  --   sols        start solutions;
  --   k           smoothing parameter to simulate small steps;
  --   a           random complex number to ensure regularity;
  --   target      target value for continuation parameter.

  -- ON RETURN :
  --   sols        solutions for t = target.

  procedure Natural_Parameter_Continuation
                 ( file : in file_type; h : in Poly_Sys; k : in natural32;
                   t0,t1 : in Complex_Number; sols : in out Solution_List );

  -- DESCRIPTION : continuation with a natural-parameter homotopy.

-- 4. POST-PROCESSING : VALIDATION

  procedure Refine_Roots 
                 ( p : in Poly_Sys; sols : in out Solution_List );
  procedure Refine_Roots 
                 ( file : in file_type; p : in Poly_Sys;
                   sols : in out Solution_List );

  -- DESCRIPTION : refines the roots and puts a report on file.

  -- ON ENTRY :
  --   file        to write results on, must be opened for output;
  --   p           a polynomial system system;
  --   sols        approximate solutions.

  -- ON RETURN :
  --   sols        refined solutions.

end PHCpack;
