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

package Cascade_Homotopies_io is

-- DESCRIPTION :
--   Cascade homotopies produce candidate witness points.
--   All solutions at the end of a cascade homotopy which have a zero slack
--   variable form a witness super set, as all witness points are contained
--   in the list of all solutions at the end of a cascade homotopy.
--   The procedures in this package write witness supersets to file.

  procedure Write_Super_Witness_Points
              ( file : in file_type;
                sols : in Standard_Complex_Solutions.Solution_List );
  procedure Write_Super_Witness_Points
              ( file : in file_type;
                sols : in DoblDobl_Complex_Solutions.Solution_List );
  procedure Write_Super_Witness_Points
              ( file : in file_type;
                sols : in QuadDobl_Complex_Solutions.Solution_List );

  -- DESCRIPTION :
  --   Writes the solutions to file.  The assumption is that all
  --   solutions in sols have a zero value for the last slack variable.

  function Append_ck ( name : string; k : natural32 ) return string;

  -- DESCRIPTION :
  --   Appends "_swk" to the name.

  procedure Write_Witness_Superset
              ( name : in string;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32 );
  procedure Write_Witness_Superset
              ( name : in string;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                k : in natural32 );
  procedure Write_Witness_Superset
              ( name : in string;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32 );
  procedure Write_Witness_Superset
              ( name : in string;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                k : in natural32 );
  procedure Write_Witness_Superset
              ( name : in string;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32 );
  procedure Write_Witness_Superset
              ( name : in string;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                k : in natural32 );

  -- DESCRIPTION :
  --   Writes the embedded polynomial system along with its
  --   candidate witness points on the file name_ck.
  --   The variable k is the dimension of the witness super set.

end Cascade_Homotopies_io;
