with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with TripDobl_Complex_Solutions;
with TripDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with PentDobl_Complex_Solutions;
with PentDobl_Complex_Solutions_io;
with OctoDobl_Complex_Solutions;
with OctoDobl_Complex_Solutions_io;
with DecaDobl_Complex_Solutions;
with DecaDobl_Complex_Solutions_io;
with HexaDobl_Complex_Solutions;
with HexaDobl_Complex_Solutions_io;

package body Test_Solutions_io is

  procedure Double_Read_Write is

    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Standard_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(Standard_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in double precision :");
    Standard_Complex_Solutions_io.put(sols);
  end Double_Read_Write;

  procedure DoblDobl_Read_Write is

    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    DoblDobl_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(DoblDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in double double precision :");
    DoblDobl_Complex_Solutions_io.put(sols);
  end DoblDobl_Read_Write;

  procedure TripDobl_Read_Write is

    sols : TripDobl_Complex_Solutions.Solution_List;

  begin
    TripDobl_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(TripDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in triple double precision :");
    TripDobl_Complex_Solutions_io.put(sols);
  end TripDobl_Read_Write;

  procedure QuadDobl_Read_Write is

    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    QuadDobl_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(QuadDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in quad double precision :");
    QuadDobl_Complex_Solutions_io.put(sols);
  end QuadDobl_Read_Write;

  procedure PentDobl_Read_Write is

    sols : PentDobl_Complex_Solutions.Solution_List;

  begin
    PentDobl_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(PentDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in penta double precision :");
    PentDobl_Complex_Solutions_io.put(sols);
  end PentDobl_Read_Write;

  procedure OctoDobl_Read_Write is

    sols : OctoDobl_Complex_Solutions.Solution_List;

  begin
    OctoDobl_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(OctoDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in octo double precision :");
    OctoDobl_Complex_Solutions_io.put(sols);
  end OctoDobl_Read_Write;

  procedure DecaDobl_Read_Write is

    sols : DecaDobl_Complex_Solutions.Solution_List;

  begin
    DecaDobl_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(DecaDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in deca double precision :");
    DecaDobl_Complex_Solutions_io.put(sols);
  end DecaDobl_Read_Write;

  procedure HexaDobl_Read_Write is

    sols : HexaDobl_Complex_Solutions.Solution_List;

  begin
    HexaDobl_Complex_Solutions_io.Read(sols);
    new_line;
    put("-> read ");
    put(HexaDobl_Complex_Solutions.Length_Of(sols),1);
    put_line(" solutions.");
    new_line;
    put_line("The solution list in hexa double precision :");
    HexaDobl_Complex_Solutions_io.put(sols);
  end HexaDobl_Read_Write;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test input/output of solution lists:");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. triple double precision");
    put_line("  3. quad double precision");
    put_line("  4. penta double precision");
    put_line("  5. octo double precision");
    put_line("  6. deca double precision");
    put_line("  7. hexa double precision");
    put("Type 0, 1, 2, 3, 4, 5, 6, or 7 to select the test : ");
    Ask_Alternative(ans,"01234567");
    new_line;
    case ans is
      when '0' => Double_Read_Write;
      when '1' => DoblDobl_Read_Write;
      when '2' => TripDobl_Read_Write;
      when '3' => QuadDobl_Read_Write;
      when '4' => PentDobl_Read_Write;
      when '5' => OctoDobl_Read_Write;
      when '6' => DecaDobl_Read_Write;
      when '7' => HexaDobl_Read_Write;
      when others => null;
    end case;
  end Main;

end Test_Solutions_io;
