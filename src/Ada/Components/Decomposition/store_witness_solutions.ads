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

package Store_Witness_Solutions is

-- DESCRIPTION :
--   This packages provides default callback procedures for use as
--   the Report_Witness_Set procedures in the Witness_Solutions packages.

  procedure Store ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 );
  procedure Store ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 );
  procedure Store ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 );
  procedure Store ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 );
  procedure Store ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 );
  procedure Store ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 );

  -- DESCRIPTION :
  --   Stores the witness set in the corresponding Witness_Solutions,
  --   corressponding to the precision and type of polynomial system.

  -- ON ENTRY :
  --   ep       embedded (Laurent) polynomial system;
  --   ws       generic points as solutions of ep;
  --   dim      dimension of the represented solution set.

end Store_Witness_Solutions;
