with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Poly_Systems;      use QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;      use QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Solutions;         use QuadDobl_Complex_Solutions;

package QuadDobl_BlackBox_Refiners is

-- DESCRIPTION :
--   A blackbox root refiner applies Newton's method with default settings
--   of the tolerances in quad double precision.
--   There are three choices to make:
--   1) silent or reporting to file;
--   2) single or multitasked;
--   3) on polynomial or Laurent systems.

  procedure Silent_Black_Box_Refine
              ( p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions given in sols of the polynomial system p,
  --   without output written to file.

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions given in sols of the polynomial system p,
  --   without output written to file, with nt tasks.

  procedure Silent_Black_Box_Refine
              ( p : in Laur_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the Laurent system p,
  --   without output written to file.

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the Laurent system p,
  --   without output written to file, with nt tasks.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the polynomial system p.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the polynomial system p
  --   with nt tasks.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Laur_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the Laurent system p.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the Laurent system p with nt tasks.

end QuadDobl_BlackBox_Refiners;
