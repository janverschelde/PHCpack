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
  --   Refines the roots in sols w.r.t. the system p,
  --   without output written to file.
  --   By default, deflation is applied.

  procedure Silent_Black_Box_Refine
              ( p : in Laur_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the roots in sols w.r.t. the system p,
  --   without output written to file.
  --   For Laurent systems, deflation is not yet available.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; p : in Poly_Sys;
                sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the roots in sols w.r.t. the system p.
  --   By default, deflation is applied.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the roots in sols w.r.t. the system p, with nt tasks.
  --   With multitasking, deflation is not yet available...

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; p : in Laur_Sys;
                sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the roots in sols w.r.t. the system p.
  --   Deflation is not yet available for Laurent systems.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List );

  -- DESCRIPTION :
  --   Refines the roots in sols w.r.t. the system p with nt tasks.
  --   Deflation is not yet available for Laurent systems.

end DoblDobl_BlackBox_Refiners;
