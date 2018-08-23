with Standard_Witness_Solutions;
with DoblDobl_Witness_Solutions;
with QuadDobl_Witness_Solutions;

package body Store_Witness_Solutions is

  procedure Store ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) is
  begin
    Standard_Witness_Solutions.Save_Embedded_System(ep,dim);
    Standard_Witness_Solutions.Save_Witness_Points(ws,dim);
  end Store;

  procedure Store ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) is
  begin
    Standard_Witness_Solutions.Save_Embedded_System(ep,dim);
    Standard_Witness_Solutions.Save_Witness_Points(ws,dim);
  end Store;

  procedure Store ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) is
  begin
    DoblDobl_Witness_Solutions.Save_Embedded_System(ep,dim);
    DoblDobl_Witness_Solutions.Save_Witness_Points(ws,dim);
  end Store;

  procedure Store ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) is
  begin
    DoblDobl_Witness_Solutions.Save_Embedded_System(ep,dim);
    DoblDobl_Witness_Solutions.Save_Witness_Points(ws,dim);
  end Store;

  procedure Store ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) is
  begin
    QuadDobl_Witness_Solutions.Save_Embedded_System(ep,dim);
    QuadDobl_Witness_Solutions.Save_Witness_Points(ws,dim);
  end Store;

  procedure Store ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) is
  begin
    QuadDobl_Witness_Solutions.Save_Embedded_System(ep,dim);
    QuadDobl_Witness_Solutions.Save_Witness_Points(ws,dim);
  end Store;

end Store_Witness_Solutions;
