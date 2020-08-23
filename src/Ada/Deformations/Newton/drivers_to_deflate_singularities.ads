with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Multprec_Complex_Vectors;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Drivers_to_Deflate_Singularities is

-- DESCRIPTION :
--  This packages provides driver routines to apply Newton's method
--  with deflation to isolated singularities.

  procedure Read_Tolerance ( tol : in out double_float );

  -- DESCRIPTION :
  --   Presenting the user with a given tolerance,
  --   the user has the opportunity to give new values for tol.

  procedure Multiple_Standard_Deflations
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys );
  procedure Multiple_Multprec_Deflations
               ( file : in file_type;
                 p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 size : in natural32 );

  -- DESCRIPTION :
  --   Interactive generation of a sequence of deflated systems,
  --   where the user is prompted for the number of multipliers.

  procedure Set_Default_Parameters
              ( symbolic,output : out boolean;
                maxitr,maxdef,nbdgts : out natural32;
	        tolerr,tolres,tolrnk : out double_float );

  -- DESCRIPTION :
  --   Sets the default values for the deflator parameters.

  -- ON RETURN :
  --   symbolic is flag for method of evaluation;
  --   output   flag for intermediate output;
  --   maxitr   upper bound on #iterations/root;
  --   maxdef   upper bound on #deflations/root;
  --   nbdgts   number of decimal places in working precision.
  --   tolerr   tolerance for error on the root;
  --   tolres   tolerance for the residual;
  --   tolrnk   tolerance to decide numerical rank.

  procedure Display_Parameters
              ( file : in file_type;
                symbolic,output : in boolean;
                maxitr,maxdef,nbdgts : in natural32;
                tolerr,tolres,tolrnk : in double_float );

  -- DESCRIPTION :
  --   Displays the menu with the current settings of the deflator.

  procedure Deflate_Singularities
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 vrblvl : in integer32 := 0 );
  procedure Deflate_Singularities
               ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 vrblvl : in integer32 := 0 );
  procedure Deflate_Singularities
               ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies deflation to the system p with solution in sols,
  --   using standard double precision and default values for
  --   the numerical tolerances.
  --   The verbose level is in vrblvl.

  procedure Deflate_Singularities
               ( file : in file_type; outfilename : in string;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 vrblvl : in integer32 := 0 );
  procedure Deflate_Singularities
               ( file : in file_type; outfilename : in string;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 vrblvl : in integer32 := 0 );
  procedure Deflate_Singularities
               ( file : in file_type; outfilename : in string;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 vrblvl : in integer32 := 0 );

  -- DESCRIPTION :
  --   Applies Newton's method with deflation to a list of solutions
  --   with standard precision coefficients.  After interactively
  --   determining the parameters, the user is no longer bothered.
  --   Systems and solutions are written to separate files, all
  --   starting with outfilename and ending with suffix _dk, where
  --   k is the number of the deflation step.
  --   The verbose level is in vrblvl.

end Drivers_to_Deflate_Singularities;
