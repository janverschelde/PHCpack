with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;

procedure ts_qdsols is

-- DESCRIPTION :
--   Test on solutions with quad double numbers.

  procedure Main is

    st_sols : Standard_Complex_Solutions.Solution_List;
    qd_sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    Read(st_sols);
    new_line;
    put("read ");
    put(Standard_Complex_Solutions.Length_Of(st_sols),1);
    put_line(" solutions ...");
    qd_sols := QuadDobl_Complex_Solutions.Create(st_sols);
    put_line("The solution list of double doubles :");
    put(qd_sols);
  end Main;

begin
  new_line;
  put_line("Testing solutions with quad double numbers.");
  new_line;
  Main;
end ts_qdsols;
