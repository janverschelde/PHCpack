with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Standard_BlackBox_Refiners is

-- DESCRIPTION :
--   A blackbox root refiner applies Newton's method with default settings
--   of the tolerances in standard double precision.
--   There are three choices to make:
--   1) silent or reporting to file;
--   2) single or multitasked;
--   3) on polynomial or Laurent systems.

-- ON POLYNOMIAL SYSTEMS :

  procedure Silent_Black_Box_Refine
              ( p : in Poly_Sys; sols : in out Solution_List;
                deflate : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions given in sols of the polynomial system p,
  --   without output written to file.
  --   If deflate is false, then no deflation is applied.
  --   The verbose level is given by verbose.

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List;
                deflate : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions given in sols of the polynomial system p,
  --   without output written to file, with nt tasks.
  --   The verbose level is given by verbose.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Poly_Sys; sols : in out Solution_List;
                deflate : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions given in sols of the polynomial system p.
  --   The verbose level is given by verbose.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Poly_Sys; sols : in out Solution_List;
                deflate : in boolean; verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the polynomial system p
  --   with nt tasks.
  --   If deflate is false, then no deflation will be applied.
  --   The verbose level is given by verbose.

-- ON LAURENT SYSTEMS :

  procedure Silent_Black_Box_Refine
              ( p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions given in sols of the Laurent system p,
  --   without output written to file.
  --   The verbose level is given by verbose.

  procedure Silent_Black_Box_Refine
              ( nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions given in sols of the Laurent system p,
  --   without output written to file, with nt tasks.
  --   The verbose level is given by verbose.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type;
                p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the Laurent system p.
  --   The verbose level is given by verbose.

  procedure Reporting_Black_Box_Refine
              ( file : in file_type; nt : in integer32;
                p : in Laur_Sys; sols : in out Solution_List;
                verbose : in integer32 := 0 );

  -- DESCRIPTION :
  --   Refines the solutions in sols of the Laurent system p with nt tasks.
  --   The verbose level is given by verbose.

end Standard_BlackBox_Refiners;
