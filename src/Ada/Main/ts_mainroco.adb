with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Main_Root_Counters;

procedure ts_mainroco is

-- DESCRIPTION :
--   This procedure calls the driver to root counting.

  lp,lq : Link_to_Poly_Sys;
  qsols : Solution_List;
  b,nt : natural32 := 0;
  file : file_type;

begin
  new_line;
  put_line("Test on driver for root counting.");
  new_line;
  get(lp);
  new_line;
  put_line("Reading the name of the output file.");
  Read_Name_and_Create_File(file);
  new_line;
  put("Give the number of tasks (0 for no multitasking) : "); get(nt);
  skip_line; -- next reading will be character
  lq := new Poly_Sys(lp'range);
  Main_Root_Counters.Polynomial_Main(file,nt,lp.all,lq.all,false,qsols,b);
end ts_mainroco;
