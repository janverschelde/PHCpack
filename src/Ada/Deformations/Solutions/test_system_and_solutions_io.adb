with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with PentDobl_Complex_Poly_Systems;
with OctoDobl_Complex_Poly_Systems;
with DecaDobl_Complex_Poly_Systems;
with HexaDobl_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with TripDobl_Complex_Solutions;
with TripDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with PentDobl_Complex_Solutions;
with PentDobl_System_and_Solutions_io;
with OctoDobl_Complex_Solutions;
with OctoDobl_System_and_Solutions_io;
with DecaDobl_Complex_Solutions;
with DecaDobl_System_and_Solutions_io;
with HexaDobl_Complex_Solutions;
with HexaDobl_System_and_Solutions_io;

package body Test_System_and_Solutions_io is

  procedure Double_Test_Get is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    Standard_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(Standard_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end Double_Test_Get;

  procedure DoblDobl_Test_Get is

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(DoblDobl_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end DoblDobl_Test_Get;

  procedure TripDobl_Test_Get is

    lp : TripDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : TripDobl_Complex_Solutions.Solution_List;

  begin
    TripDobl_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(TripDobl_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end TripDobl_Test_Get;

  procedure QuadDobl_Test_Get is

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(QuadDobl_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end QuadDobl_Test_Get;

  procedure PentDobl_Test_Get is

    lp : PentDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : PentDobl_Complex_Solutions.Solution_List;

  begin
    PentDobl_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(PentDobl_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end PentDobl_Test_Get;

  procedure OctoDobl_Test_Get is

    lp : OctoDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : OctoDobl_Complex_Solutions.Solution_List;

  begin
    OctoDobl_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(OctoDobl_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end OctoDobl_Test_Get;

  procedure DecaDobl_Test_Get is

    lp : DecaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DecaDobl_Complex_Solutions.Solution_List;

  begin
    DecaDobl_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(DecaDobl_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end DecaDobl_Test_Get;

  procedure HexaDobl_Test_Get is

    lp : HexaDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : HexaDobl_Complex_Solutions.Solution_List;

  begin
    HexaDobl_System_and_Solutions_io.get(lp,sols);
    put("-> read ");
    put(natural32(lp'last),1); put(" polynomials and ");
    put(HexaDobl_Complex_Solutions.Length_Of(sols),1);
    put(" solutions.");
  end HexaDobl_Test_Get;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test reading system with solutions :");
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
    case ans is
      when '0' => Double_Test_Get;
      when '1' => DoblDobl_Test_Get;
      when '2' => TripDobl_Test_Get;
      when '3' => QuadDobl_Test_Get;
      when '4' => PentDobl_Test_Get;
      when '5' => OctoDobl_Test_Get;
      when '6' => DecaDobl_Test_Get;
      when '7' => HexaDobl_Test_Get;
      when others => null;
    end case;
  end Main;

end Test_System_and_Solutions_io;
