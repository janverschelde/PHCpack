with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Quad_Double_Poly_Systems;
with Quad_Double_Poly_SysFun;
with Quad_Double_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_Quad_Sweepers is

-- DESCRIPTION :
--   A sweeper is a path tracker which monitors the orientation of
--   the tangent for quadratic turning points and the determinant
--   for general types of singularities, in quad double precision.

  function Real_Part ( x : QuadDobl_Complex_Vectors.Vector )
                     return Quad_Double_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the real parts of all components of x.

  procedure Interactive_Real_Sweep
              ( nq,nv : in natural32; eigval : in boolean;
                p : in Quad_Double_Poly_Systems.Poly_Sys;
                sols : in Solution_List );
  procedure Interactive_Real_Sweep
              ( nq,nv : in natural32; eigval : in boolean;
                p : in Quad_Double_Poly_Systems.Poly_Sys;
                sols : in Solution_List; repsols : in out Solution_List );
  procedure Interactive_Complex_Sweep
              ( nq,nv : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for an initial point and then starts
  --   the real continuation using the homotopy in p.
  --   In this interactive version, the user can monitor the progress
  --   of the corrector and the shooting method.

  -- ON ENTRY :
  --   nq       number of equations in p;
  --   nv       number of variables, parameter t included;
  --   eigval   eigenvalues will be computed if and only if eigval is true;
  --   p        a real system with nq equations in nv variables;
  --   sols     solutions to start the sweep at.

  -- ON RETURN :
  --   repsols  solutions along path recorded by user.

  procedure Silent_Real_Sweep
              ( eigval : in boolean;
                nq,nv : in natural32; t : in quad_double;
                f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jf : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Quad_Double_Vectors.Vector;
                dx : out Quad_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs a sweep starting at the given solution, till either
  --   a singularity is encountered, or till the target value for the
  --   parameter is reached, or till the number of steps is exhausted,
  --   without intermediate output to file or screen.

  -- ON ENTRY :
  --   output   for intermediate diagnostics during continuation;
  --   eigval   flag for eigenvalue computations;
  --   nq       number of equations in p;
  --   nv       number of variables, parameter t included;
  --   t        target value for x(x'last);
  --   f        a real system with nq equations in nv variables;
  --   jf       Jacobian matrix of f, also in evaluable format;
  --   x        current regular solution along a path.

  -- ON RETURN :
  --   x        solution at end is either a singularity, or
  --            corresponding to t, or somewhere along a path
  --            that required max number of steps;
  --   dx       tangent vector corresponding to the solution.
  --
  procedure Start_Real_Sweep
              ( eigval : in boolean;
                nq,nv : in natural32; t : in quad_double;
                f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jf : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Quad_Double_Vectors.Vector;
                dx : out Quad_Double_Vectors.Vector );
  procedure Start_Real_Sweep
              ( file : in file_type; output,eigval : in boolean;
                nq,nv : in natural32; t : in quad_double;
                f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jf : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out Quad_Double_Vectors.Vector;
                dx : out Quad_Double_Vectors.Vector );
  procedure Start_Complex_Sweep
              ( nq,nv : in natural32; t : in quad_double;
                f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector );
  procedure Start_Complex_Sweep
              ( file : in file_type; output : in boolean;
                nq,nv : in natural32; t : in quad_double;
                f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in out QuadDobl_Complex_Vectors.Vector;
                dx : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Performs a sweep starting at the given solution, till either
  --   a singularity is encountered, or till the target value for the
  --   parameter is reached, or till the number of steps is exhausted.

  -- ON ENTRY :
  --   file     must be opened for output;
  --   output   for intermediate diagnostics during continuation;
  --   eigval   flag for eigenvalue computations;
  --   nq       number of equations in p;
  --   nv       number of variables, parameter t included;
  --   t        target value for x(x'last);
  --   f        a real system with nq equations in nv variables;
  --   jf       Jacobian matrix of f, also in evaluable format;
  --   x        current regular solution along a path.

  -- ON RETURN :
  --   x        solution at end is either a singularity, or
  --            corresponding to t, or somewhere along a path
  --            that required max number of steps;
  --   dx       tangent vector corresponding to the solution.

-- DRIVER ROUTINES :

  procedure Ask_Sweep_Parameters ( target : out quad_double );

  -- DESCRIPTION :
  --   Prompts the user for a target value of the continuation variable.

  procedure Ask_Sweep_Parameters
              ( target : out quad_double; output : out boolean );

  -- DESCRIPTION :
  --   Prompts the user for the sweep parameters:
  --   (1) target is the target value for the continuation variable;
  --   (2) output is a flag for output of the corrector.

  procedure Ask_Sweep_Parameters
              ( target : out quad_double; output,eigval : out boolean );

  -- DESCRIPTION :
  --   Prompts the user for the sweep parameters:
  --   (1) target is the target value for the continuation variable;
  --   (2) output flags intermediate output of the corrector;
  --   (3) eigval is a flag to monitor the eigenvalues of the Jacobian
  --       matrix along a path.

  procedure Run_Sweep
              ( nq,nv : in natural32;
                p : in Quad_Double_Poly_Systems.Poly_Sys;
                sols : in Solution_List );

  -- DESCRIPTION :
  --   Runs the sweep on a polynomial system with real coefficients,
  --   using the solutions in sols to start at.

  procedure Run_Sweep
              ( nq,nv : in natural32;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in Solution_List );

  -- DESCRIPTION :
  --   Runs the sweep on a polynomial system with complex coefficients.

end QuadDobl_Quad_Sweepers;
