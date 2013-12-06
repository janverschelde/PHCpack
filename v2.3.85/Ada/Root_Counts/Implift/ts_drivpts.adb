with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Drivers_for_Vertex_Points;          use Drivers_for_Vertex_Points;

procedure ts_drivpts is

-- DESCRIPTION :
--   This procedure calls the driver for extracting vertex points.

  lp : Link_to_Poly_Sys;
  file : file_type;

begin
  new_line;
  put_line("Test on extracting the vertices.");
  new_line;
  get(lp);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  Vertex_Points(file,lp.all);
end ts_drivpts;
