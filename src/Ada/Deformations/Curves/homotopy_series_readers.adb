with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_System_and_Solutions_io;
with DoblDobl_Homotopy;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Homotopy;
with QuadDobl_System_and_Solutions_io;

package body Homotopy_Series_Readers is

  procedure Standard_Reader
              ( nbequ : out integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 ) is

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := tpow;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    nbequ := target'last;
    new_line;
    put_line("Reading the start system and its solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    Standard_Homotopy.Create(target.all,start.all,k,gamma);
  end Standard_Reader;

  procedure DoblDobl_Reader
              ( nbequ : out integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 ) is

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := tpow;
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    nbequ := target'last;
    new_line;
    put_line("Reading the start system and its solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    DoblDobl_Homotopy.Create(target.all,start.all,k,gamma);
  end DoblDobl_Reader;

  procedure QuadDobl_Reader
              ( nbequ : out integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                tpow : in natural32 := 2 ) is

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := tpow;
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    nbequ := target'last;
    new_line;
    put_line("Reading the start system ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    QuadDobl_Homotopy.Create(target.all,start.all,k,gamma);
  end QuadDobl_Reader;

end Homotopy_Series_Readers;
