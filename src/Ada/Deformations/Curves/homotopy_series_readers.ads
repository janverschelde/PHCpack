with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Homotopy_Series_Readers is

-- DESCRIPTION :
--   Provides interactive procedures to setup of homotopies of series,
--   in double, double double, and quad double precision.
--   The homotopy is an artificial parameter homotopy.

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 );
  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 );
  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 );

  -- DESCRIPTION :
  --   Prompts for a target system, a start system with start solutions.
  --   The target and start system are stored in the Homotopy package.

  -- ON ENTRY :
  --   tpow     power of the continuation parameter
  --            in the artificial parameter homotopy.

  -- ON RETURN :
  --   nbequ    number of equations in the systems in the homotopy;
  --   sols     start solutions in the homotopy.

end Homotopy_Series_Readers;
