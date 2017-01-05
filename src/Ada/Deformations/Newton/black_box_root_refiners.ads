with text_io;                            use text_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;

package Black_Box_Root_Refiners is

-- DESCRIPTION :
--   Wraps the root refiners with basic settins for the tolerances.

  procedure Refine_Roots
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in out Standard_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the result of the root refining on the system p,
  --   starting at the roots in sols to file.

end Black_Box_Root_Refiners;
