with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Black_Box_Root_Refiners is

-- DESCRIPTION :
--   Wraps the root refiners with basic settins for the tolerances.

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the result of the root refining on the system p,
  --   starting at the roots in sols to file.
  --   Computations occur in double, double double, or quad double precision.

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Laur_Systems.Laur_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List );
  procedure Refine_Roots
               ( file : in file_type;
                 p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Applies the root refinement for Laurent systems,
  --   which is more basic as no deflation is available.

end Black_Box_Root_Refiners;
