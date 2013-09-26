with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Drivers_for_Root_Counts is

-- DESCRIPTION :
--   This package offers interactive drivers to count the roots
--   and to construct a start system for general polynomial systems
--   and for Laurent systems.

  procedure Driver_for_Root_Counts
               ( file : in file_type; p,q : in out Poly_Sys;
                 own : in boolean;
                 qsols : in out Solution_List; roco : out natural32 );
  procedure Driver_for_Root_Counts
               ( file : in file_type; p,q : in out Laur_Sys;
                 qsols : in out Solution_List; roco : out natural32 );

  -- DESCRIPTION :
  --   This procedure implements an interactive driver for several
  --   root counting methods.  If the system is a Laurent system,
  --   then only the polyhedral methods are available.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   p         a polynomial system;
  --   own       if true, then the user has the possibility to give
  --             an own start system.

  -- ON RETURN :
  --   p         has eventually been made homogeneous;
  --   q         a start system based on the chosen root count;
  --   qsols     the solutions of q;
  --   roco      the root count.

end Drivers_for_Root_Counts;
