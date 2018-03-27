with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;

package Random_Test_Points is

-- DESCRIPTION :
--   To test the homotopy membership test for a point to belong to
--   a solution set represented by a witness set, a random point is
--   generated, computed by tracking one solution path in a homotopy
--   going from the given slice to another slice.
--   The procedures in this package compute one random point
--   for witness sets for ordinary and Laurent polynomial systems,
--   in double, double double, and quad double precision.

  function Standard_Random_Point
             ( file : file_type;
               p : Standard_Complex_Poly_Systems.Poly_Sys;
               s : Standard_Complex_Solutions.Solution_List;
               d : natural32 )
             return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

  function Standard_Random_Point
             ( file : file_type;
               p : Standard_Complex_Laur_Systems.Laur_Sys;
               s : Standard_Complex_Solutions.Solution_List;
               d : natural32 )
             return Standard_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

  function DoblDobl_Random_Point
             ( file : file_type;
               p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : DoblDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

  function DoblDobl_Random_Point
             ( file : file_type;
               p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               s : DoblDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return DoblDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

  function QuadDobl_Random_Point
             ( file : file_type;
               p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : QuadDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

  function QuadDobl_Random_Point
             ( file : file_type;
               p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               s : QuadDobl_Complex_Solutions.Solution_List;
               d : natural32 )
             return QuadDobl_Complex_Solutions.Solution_List;

  -- DESCRIPTION :
  --   For a given witness set, takes the first point and applies
  --   path tracking to another set of slices to compute another point.

  -- ON ENTRY :
  --   file    for output during sampling;
  --   p       embedded polynomial system;
  --   s       generic points on the solution set;
  --   d       dimension of the solution set;

end Random_Test_Points;
