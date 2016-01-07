with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Double_Vectors;
with Double_Double_VecVecs;
with Quad_Double_Vectors;
with Quad_Double_VecVecs;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Drivers_for_Path_Directions is

-- DESCRIPTION :
--   This package provides driver routines for the computation of the
--   direction of the solution paths.

  procedure Init_Path_Directions
               ( n,nv : in natural32;
                 v : in out Standard_Floating_VecVecs.Link_to_VecVec;
                 errv : in out Standard_Floating_Vectors.Link_to_Vector );
  procedure Init_Path_Directions
               ( n,nv : in natural32;
                 v : in out Double_Double_VecVecs.Link_to_VecVec;
                 errv : in out Double_Double_Vectors.Link_to_Vector );
  procedure Init_Path_Directions
               ( n,nv : in natural32;
                 v : in out Quad_Double_VecVecs.Link_to_VecVec;
                 errv : in out Quad_Double_Vectors.Link_to_Vector );

  -- DESCRIPTION :
  --   Initializes the data for path directions,
  --   in standard double, double double, or quad double precision.

  -- ON ENTRY :
  --   n         dimension of the solution vectors;
  --   nv        number of solution paths.
 
  -- ON RETURN :
  --   v         nv zero vectors of range 1..n;
  --   errv      nv entries equal to 1.0.

  procedure Toric_Continue
               ( file : in file_type;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 proj,report : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Standard_Floating_VecVecs.VecVec;
                 errv : in out Standard_Floating_Vectors.Vector;
                 target : in Standard_Complex_Numbers.Complex_Number );
  procedure Toric_Continue
               ( file : in file_type;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 proj,report : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Double_Double_VecVecs.VecVec;
                 errv : in out Double_Double_Vectors.Vector;
                 target : in DoblDobl_Complex_Numbers.Complex_Number );
  procedure Toric_Continue
               ( file : in file_type;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 proj,report : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Quad_Double_VecVecs.VecVec;
                 errv : in out Quad_Double_Vectors.Vector;
                 target : in QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Performs the continuation with computation of path directions,
  --   in standard double, double double, or quad double precision.

  procedure Write_Direction
               ( file : in file_type; w : in integer32;
                 v : in Standard_Floating_Vectors.Vector;
                 error : in double_float; i : in integer32 );
  procedure Write_Direction
               ( file : in file_type; w : in integer32;
                 v : in Double_Double_Vectors.Vector;
                 error : in double_double; i : in integer32 );
  procedure Write_Direction
               ( file : in file_type; w : in integer32;
                 v : in Quad_Double_Vectors.Vector;
                 error : in quad_double; i : in integer32 );

  -- DESCRIPTION :
  --   Defines the output format to write a direction v,
  --   with the estimated winding number w
  --   and the estimated error for the i-th path,
  --   in standard double, double double, or quad double precision.

  procedure Write_Directions 
               ( file : in file_type;
                 w : in Standard_Integer_Vectors.Vector;
                 v : in Standard_Floating_VecVecs.VecVec;
                 errv : in Standard_Floating_Vectors.Vector );
  procedure Write_Directions 
               ( file : in file_type;
                 w : in Standard_Integer_Vectors.Vector;
                 v : in Double_Double_VecVecs.VecVec;
                 errv : in Double_Double_Vectors.Vector );
  procedure Write_Directions 
               ( file : in file_type;
                 w : in Standard_Integer_Vectors.Vector;
                 v : in Quad_Double_VecVecs.VecVec;
                 errv : in Quad_Double_Vectors.Vector );

  -- DESCRIPTION :
  --   Writes the directions of the paths in v with corresponding
  --   estimated winding numbers in w and their errors to the file,
  --   in standard double, double double, or quad double precision.

end Drivers_for_Path_Directions;
