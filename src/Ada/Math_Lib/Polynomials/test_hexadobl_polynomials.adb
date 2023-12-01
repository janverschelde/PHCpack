with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Hexa_Double_Polynomials;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Laurentials;
with HexaDobl_Polynomial_Convertors;     use HexaDobl_Polynomial_Convertors;

package body Test_HexaDobl_Polynomials is

  procedure Test_Polynomial_Convertor is

    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    dd_rp : Hexa_Double_Polynomials.Poly;
    dd_cp : HexaDobl_Complex_Polynomials.Poly;

  begin
    new_line;
    put_line("Testing conversion with standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    new_line;
    put_line("Your polynomial : "); put(p); new_line;
    put_line("converting to a hexa double polynomial ...");
    dd_rp := Standard_Polynomial_to_Hexa_Double(p);
    put_line("converting to a double complex polynomial ...");
    dd_cp := Standard_Polynomial_to_HexaDobl_Complex(p);
    Standard_Complex_Polynomials.Clear(p);
    put_line("Writing original back after conversions ...");
    p := Hexa_Double_to_Standard_Polynomial(dd_rp);
    put("hexa double : "); put(p); new_line;
    Standard_Complex_Polynomials.Clear(p);
    p := HexaDobl_Complex_to_Standard_Polynomial(dd_cp);
    put("HexaDobl complex : "); put(p); new_line;
    Symbol_Table.Clear;
  end Test_Polynomial_Convertor;

  procedure Test_Laurent_Polynomial_Convertor is

    n : natural32 := 0;
    p : Standard_Complex_Laurentials.Poly;
    dd_p : HexaDobl_Complex_Laurentials.Poly;

  begin
    new_line;
    put_line("Testing conversion with standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    new_line;
    put_line("Your polynomial : "); put(p); new_line;
    put_line("converting to a double complex polynomial ...");
    dd_p := Standard_Laurential_to_HexaDobl_Complex(p);
    Standard_Complex_Laurentials.Clear(p);
    put_line("Writing original back after conversions ...");
    p := HexaDobl_Complex_to_Standard_Laurential(dd_p);
    put("HexaDobl complex : "); put(p); new_line;
    Symbol_Table.Clear;
  end Test_Laurent_Polynomial_Convertor;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing polynomials with hexa double coefficients ...");
    new_line;
    put_line("MENU with testing operations :");
    put_line("  1. test convertor for polynomials");
    put_line("  2. test convertor for Laurent polynomials");
    put("Type 1 or 2 to choose a test : ");
    Ask_Alternative(ans,"12");
    case ans is
      when '1' => Test_Polynomial_Convertor;
      when '2' => Test_Laurent_Polynomial_Convertor;
      when others => null;
    end case;
  end Main;

end Test_HexaDobl_Polynomials;
