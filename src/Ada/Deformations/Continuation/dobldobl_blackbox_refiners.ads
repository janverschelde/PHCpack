with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;      use DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;

package DoblDobl_BlackBox_Refiners is

-- DESCRIPTION :
--   A blackbox root refiner applies Newton's method with default settings
--   of the tolerances in double double precision.
--   There are three choices to make:
--   1) silent or reporting to file;
--   2) single or multitasked;
--   3) on polynomial or Laurent systems.

  procedure Silent_Black_Box_Refine
              ( p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solution given in sols of the polynomial system p,
  --   without output written to file.

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solution given in sols of the polynomial system p,
  --   without output written to file, with nt tasks.

  procedure Silent_Black_Box_Refine
              ( p : in Laur_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solution given in sols of the Laurent system p,
  --   without output written to file.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the system polynomial p.

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

end DoblDobl_BlackBox_Refiners;
