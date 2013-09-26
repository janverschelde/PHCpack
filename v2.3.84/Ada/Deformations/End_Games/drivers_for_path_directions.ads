with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_VecVecs;          use Standard_Floating_VecVecs;  
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;

package Drivers_for_Path_Directions is

-- DESCRIPTION :
--   This package provides driver routines for the computation of the
--   direction of the solution paths.

  procedure Init_Path_Directions
               ( n,nv : in natural32; v : in out Link_to_VecVec;
                 errv : in out Link_to_Vector );

  -- DESCRIPTION :
  --   Initializes the data for path directions.

  -- ON ENTRY :
  --   n         dimension of the solution vectors;
  --   nv        number of solution paths.
 
  -- ON RETURN :
  --   v         nv zero vectors of range 1..n;
  --   errv      nv entries equal to 1.0.

  procedure Toric_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj,report : in boolean; v : in out VecVec;
                 errv : in out Vector; target : in Complex_Number );

  -- DESCRIPTION :
  --   Performs the continuation with computation of path directions.

  procedure Write_Directions 
               ( file : in file_type; v : in VecVec; errv : in Vector );

  -- DESCRIPTION :
  --   Writes the directions of the paths with their errors to the file.

end Drivers_for_Path_Directions;
