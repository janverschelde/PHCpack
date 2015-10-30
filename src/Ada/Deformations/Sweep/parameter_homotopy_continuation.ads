with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Floating_Poly_Systems;
with Standard_Floating_Poly_SysFun;
with Standard_Floating_Jaco_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with Double_Double_Poly_Systems;
with Double_Double_Poly_SysFun;
with Double_Double_Jaco_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with Quad_Double_Poly_Systems;
with Quad_Double_Poly_SysFun;
with Quad_Double_Jaco_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Parameter_Homotopy_Continuation is

-- DESCRIPTION :
--   Parameter continuation assumes a system with strictly fewer equations
--   than unknowns.  The definition of start and target values for the
--   parameters lead to naturally to interpolation between the two.
--   This interpolation may be seen as the default function to instantiate
--   the general parameter continuation routine with.

  function Define_Start
             ( v : Standard_Complex_Vectors.Vector;
               ip : Standard_Integer_Vectors.Vector )
             return Standard_Complex_Vectors.Vector;
  function Define_Start
             ( v : DoblDobl_Complex_Vectors.Vector;
               ip : Standard_Integer_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector;
  function Define_Start
             ( v : QuadDobl_Complex_Vectors.Vector;
               ip : Standard_Integer_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the elements in v indexed by the parameters in ip, for
  --   vectors v in standard double, double double, or quad double precision.

  function Define_Complex_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;
  function Define_Complex_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return DoblDobl_Complex_Vectors.Vector;
  function Define_Complex_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Given the indices to the parameters, this function prompts the
  --   user for complex target values of the parameters,
  --   in standard double, double double, or quad double precision.

  function Define_Real_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return Standard_Complex_Vectors.Vector;
  function Define_Real_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return DoblDobl_Complex_Vectors.Vector;
  function Define_Real_Target
              ( ip : Standard_Integer_Vectors.Vector )
              return QuadDobl_Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   Given the indices to the parameters, this function prompts the
  --   user for real target values of the parameters,
  --   in standard double, double double, or quad double precision.

  procedure Write_Complex_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in Standard_Complex_Vectors.Vector );
  procedure Write_Complex_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in DoblDobl_Complex_Vectors.Vector );
  procedure Write_Complex_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the values for the parameters in a nice format to file,
  --   for standard double, double double and quad double precision.
  --
  procedure Write_Real_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in Standard_Complex_Vectors.Vector );
  procedure Write_Real_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in DoblDobl_Complex_Vectors.Vector );
  procedure Write_Real_Parameter_Values
              ( file : in file_type;
                labels : in Standard_Integer_Vectors.Vector;
                values : in QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the values for the parameters in a nice format to file,
  --   omitting the imaginary parts of the complex numbers,
  --   in standard double, double double, and quad double precision.

  procedure Determine_Parameter_Values
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List;
                par : in Standard_Integer_Vectors.Vector;
                isreal : in out boolean;
                start,target : out Standard_Complex_Vectors.Vector );
  procedure Determine_Parameter_Values
              ( file : in file_type;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                par : in Standard_Integer_Vectors.Vector;
                isreal : in out boolean;
                start,target : out DoblDobl_Complex_Vectors.Vector );
  procedure Determine_Parameter_Values
              ( file : in file_type;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                par : in Standard_Integer_Vectors.Vector;
                isreal : in out boolean;
                start,target : out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Determines the values for the parameters,
  --   in standard double, double double, or quad double precision.

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

  procedure Coefficient_Parameter_Homotopy_Continuation
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 );
  procedure Coefficient_Parameter_Homotopy_Continuation
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 );
  procedure Coefficient_Parameter_Homotopy_Continuation
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 );

  -- DESCRIPTION :
  --   Interactive set up a coefficient-parameter homotopy,
  --   in standard double, double double, or quad double precision,
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
  function Complex_Sweep_Line
              ( n,k : integer32;
                start,target : DoblDobl_Complex_Numbers.Complex_Number )
              return DoblDobl_Complex_Polynomials.Poly;
  function Complex_Sweep_Line
              ( n,k : integer32;
                start,target : QuadDobl_Complex_Numbers.Complex_Number )
              return QuadDobl_Complex_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns a polynomial in n+1 variables, in standard double,
  --   double double, or quad double precision,  where the k-th
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
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Does a sweep in complex arithmetic, in standard double,
  --   double double, or quad double precision,  till the target
  --   at t = 1 is reached or till a singularity is encountered,
  --   starting at the given solution.

  -- ON ENTRY :
  --   file     output file for intermediate diagnostics and results;
  --   output   if intermediate output along a path is needed;
  --   k        index of the solution that is swept (for output);
  --   nq       number of equations in the homotopy h;
  --   nv       number of variables in the homotopy h;
  --   h        a sweeping homotopy with as many sweeping lines
  --            as the total number of parameters;
  --   s        values for variables and parameters at t = 0.

  -- ON RETURN :
  --   s        updated values for the solution, either at t = 1,
  --            or in case t < 1, close to a singular solution.

  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Complex_Poly_Systems.Poly_Sys;
                f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in Standard_Complex_Solutions.Link_to_Solution );
  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure Run_Complex_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Does a sweep in complex arithmetic, in standard double,
  --   double double, or quad double precision,  till the target
  --   at t = 1 is reached or till a singularity is encountered,
  --   starting at the given solution.

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

  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Floating_Poly_Systems.Poly_Sys;
                s : in Standard_Complex_Solutions.Link_to_Solution );
  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Double_Double_Poly_Systems.Poly_Sys;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Quad_Double_Poly_Systems.Poly_Sys;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Does a sweep in real arithmetic, in standard double,
  --   double double or quad double precision, till the target
  --   at t = 1 is reached or till a singularity is encountered,
  --   starting at the given solution.

  -- REQUIRED : s.v is real because this is a real sweep.

  -- ON ENTRY :
  --   file     output file for intermediate diagnostics and results;
  --   output   if intermediate output along a path is needed;
  --   k        index of the solution that is swept (for output);
  --   nq       number of equations in the homotopy h;
  --   nv       number of variables in the homotopy h;
  --   h        a sweeping homotopy with as many sweeping lines
  --            as the total number of parameters;
  --   s        values for variables and parameters at t = 0.

  -- ON RETURN :
  --   s        updated values for the solution, either at t = 1,
  --            or in case t < 1, close to a singular solution.

  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Standard_Floating_Poly_Systems.Poly_Sys;
                f : in Standard_Floating_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Floating_Jaco_Matrices.Eval_Jaco_Mat;
                s : in Standard_Complex_Solutions.Link_to_Solution );
  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Double_Double_Poly_Systems.Poly_Sys;
                f : in Double_Double_Poly_SysFun.Eval_Poly_Sys;
                jf : in Double_Double_Jaco_Matrices.Eval_Jaco_Mat;
                s : in DoblDobl_Complex_Solutions.Link_to_Solution );
  procedure Run_Real_Sweep
              ( file : in file_type; output : in boolean;
                k,nq,nv : in natural32;
                h : in Quad_Double_Poly_Systems.Poly_Sys;
                f : in Quad_Double_Poly_SysFun.Eval_Poly_Sys;
                jf : in Quad_Double_Jaco_Matrices.Eval_Jaco_Mat;
                s : in QuadDobl_Complex_Solutions.Link_to_Solution );

  -- DESCRIPTION :
  --   Does a sweep in real arithmetic, in standard double,
  --   double double or quad double precision, till the target
  --   at t = 1 is reached or till a singularity is encountered,
  --   starting at the given solution.

  -- REQUIRED : s.v is real because this is a real sweep.

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
  procedure Sweep
              ( file : in file_type; isreal : in out boolean;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 );
  procedure Sweep
              ( file : in file_type; isreal : in out boolean;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                nb_equ,nb_unk,nb_par : in integer32 );

  -- DESCRIPTION :
  --   Determines the start and target values for the parameters
  --   and then calls the path trackers to do a sweep,
  --   in standard double, double double, or quad double precision.

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
