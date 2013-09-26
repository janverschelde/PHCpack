with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Witness_Sets_io;                    use Witness_Sets_io;
with Driver_to_Factor_Components;

procedure ts_irco is

-- DESCRIPTION :
--   Calls the driver to factor pure dimensional solution sets.

  file : file_type;
  lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  sols : Standard_Complex_Solutions.Solution_List;
  dim : natural32;

begin
  new_line;
  put_line
    ("Factor a pure dimensional solution set into irreducible components");
  Standard_Read_Embedding(lp,sols,dim);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  Driver_to_Factor_Components(file,lp.all,sols,dim);
end ts_irco;
