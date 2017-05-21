with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Homotopy_Series_Readers is

-- DESCRIPTION :
--   Provides interactive procedures to setup of homotopies of series,
--   in double, double double, and quad double precision.

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List );
  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List );
  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system with
  --   start solutions, returned in sols.
  --   The number of equations is returned in nbequ.
  --   The target and start system are stored in the Homotopy package.

end Homotopy_Series_Readers;
