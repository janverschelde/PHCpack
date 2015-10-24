with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with Standard_Floating_Poly_Systems;
with Standard_Floating_Poly_SysFun;
with Standard_Floating_Jaco_Matrices;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;

package Parameter_Homotopy_Continuation is

-- DESCRIPTION :
--   Parameter continuation assumes a system with strictly fewer equations
--   than unknowns.  The definition of start and target values for the
--   parameters lead to naturally to interpolation between the two.
--   This interpolation may be seen as the default function to instantiate
--   the general parameter continuation routine with.

  function Define_Start ( v : Standard_Complex_Vectors.Vector;
                          ip : Standard_Integer_Vectors.Vector )
                        return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the elements in v indexed by the parameters in ip.

  function Define_Complex_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;
  function Define_Real_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Given the indices to the parameters, this function prompts the
  --   user for complex or real target values of the parameters.

  procedure Write_Complex_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in Standard_Complex_Vectors.Vector );
  procedure Write_Real_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the values for the parameters in a nice format to file.
  --   The _Real_ omits the imaginary parts.

  procedure Determine_Parameter_Values
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List;
                par : in Standard_Integer_Vectors.Vector;
                isreal : in out boolean;
                start,target : out Standard_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Determines the values for the parameters.

  -- REQUIRED : all solutions have the same value for the parameters.

  -- ON ENTRY :
  --   file     for writing start and target values to output file;
  --   sols     solution list with all the same values for the parameters;
  --   par      indices to the parameters in the solutions;
  --   isreal   true if the polynomial system has real coefficients.

  -- ON RETURN :
  --   isreal   true if start values for the parameters were all real,
  --            false otherwise or if isreal was already false;
  --   start    start values for the parameters;
  --   target   target values for the parameters, will be real if isreal.

  function Interpolate ( a,b : Standard_Complex_Vectors.Vector;
                         t : Standard_Complex_Numbers.Complex_Number )
                       return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns (1-t)*a + t*b, for t between 0 and 1.  Is the 1st
  --   default to use in the instantiation of Parameter_Continuation.

  function Circulate ( a,b : Standard_Complex_Vectors.Vector;
                       gamma,t : Standard_Complex_Numbers.Complex_Number )
                     return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns (1-s)*a + s*b, for s = t + gamma*t*(1-t).  This move
  --   through the parameter space ensures complex arithmetic.

  function Differentiate ( a,b : Standard_Complex_Vectors.Vector )
                         return Standard_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Just returns b-a as this is the derivative of the interpolate
  --   function of the parameters.  Is the 2nd default to use in the
  --   instantiation of Parameter_Continuation.

  generic
    with function Evaluate_Parameters 
                    ( t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector;
    -- returns value of parameters at t
    with function Differentiate_Parameters
                    ( t : Standard_Complex_Numbers.Complex_Number )
                    return Standard_Complex_Vectors.Vector;
    -- returns derivatives of parameters with respect to t
  procedure Standard_Parameter_Continuation
              ( file : in file_type;
                n : in integer32;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                pars : in Standard_Integer_Vectors.Vector;
                vars : in Standard_Integer_Vectors.Vector;
		sols : in out Standard_Complex_Solutions.Solution_List;
                output : in boolean );

  -- DESCRIPTION :
  --   Executes the parameter continuation for the homotopy defined by
  --   the system p and for parameters defined by the generics functions.

  -- ON ENTRY :
  --   file     for intermediate output and diagnostics;
  --   n        number of unknowns and parameters in p.
  --   p        polynomial system with parameters;
  --   pars     indices to the parameter variables in p;
  --   vars     indices to the unknown variables in p;
  --   sols     start solutions for start values of parameters;
  --   output   true if intermediate output during continuation.

  -- ON RETURN :
  --   sols     solutions for target values of parameters.

  procedure Coefficient_Parameter_Homotopy_Continuation
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 );

  -- DESCRIPTION :
  --   Interactive set up a coefficient-parameter homotopy 
  --   and calls then the path trackers.

  -- REQUIRED :
  --   All solutions have the same values for the parameters.

  -- ON ENTRY :
  --   file     file created for output;
  --   p        a polynomial system;
  --   sols     solutions for particular values of the parameters;
  --   nb_equ   number of equations in the system;
  --   nb_unk   number of unknowns in the system;
  --   nb_par   number of parameters: nb_unk - nb_equ.

  function Complex_Sweep_Line
              ( n,k : integer32;
                start,target : Standard_Complex_Numbers.Complex_Number )
              return Standard_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a polynomial in n+1 variables where the k-th
  --   parameter moves from start to target as the (n+1)-th
  --   continuation parameter moves from 0 to 1.

  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Link_to_Solution );
  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in Standard_Complex_Solutions.Link_to_Solution );
  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Floating_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Link_to_Solution );
  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Floating_Poly_Systems.Poly_Sys;
                f : in Standard_Floating_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Floating_Jaco_Matrices.Eval_Jaco_Mat;
                s : in Standard_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Does a sweep in complex or real arithmetic till the target
  --   at t = 1 is reached or till a singularity is encountered,
  --   starting at the given solution.

  -- REQUIRED : s.v is real in a real sweep.

  -- ON ENTRY :
  --   file     output file for intermediate diagnostics and results;
  --   output   if intermediate output along a path is needed;
  --   k        index of the solution that is swept (for output);
  --   nq       number of equations in the homotopy h;
  --   nv       number of variables in the homotopy h;
  --   h        a sweeping homotopy with as many sweeping lines
  --            as the total number of parameters;
  --   f        homotopy in evaluable form;
  --   jf       evaluable Jacobi map for the homotopy;
  --   s        values for variables and parameters at t = 0.

  -- ON RETURN :
  --   s        updated values for the solution, either at t = 1,
  --            or in case t < 1, close to a singular solution.
 
  procedure Sweep
              ( file : in file_type; isreal : in out boolean;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 );

  -- DESCRIPTION :
  --   Determines the start and target values for the parameters
  --   and then calls the path trackers to do a sweep.

  -- ON ENTRY :
  --   file     for writing values of parameters to output file
  --   isreal   true if the input coefficients of the system
  --            all have nonzero imaginary parts;
  --   p        system with more variables than equations;
  --   sols     solutions, all with the same value(s) for parameter(s);
  --   nb_equ   number of equations in the system;
  --   nb_unk   total number of unknowns in the solutions;
  --   nb_par   number of prameters equals nb_unk - nb_equ.

  -- ON RETURN :
  --   isreal   true if isreal on input and all starting values
  --            for the input parameters are real, false otherwise.

  function Show_Menu return character;

  -- DESCRIPTION :
  --   Shows the user the menu for parameter continuation and prompts
  --   for '1' (complex parameter continuation) or '2' (real sweep).
  --   The choice of the user is returned as either '1' or '2'.

end Parameter_Homotopy_Continuation;
