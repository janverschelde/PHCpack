with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Penta_Double_Numbers_io;            use Penta_Double_Numbers_io;
with Multprec_PentDobl_Convertors;       use Multprec_PentDobl_Convertors;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with PentDobl_Complex_Numbers;
with PentDobl_Complex_Numbers_io;        use PentDobl_Complex_Numbers_io;
with PentDobl_Complex_Numbers_cv;        use PentDobl_Complex_Numbers_cv;

package body Test_Multprec_PentDobl_Casts is

  procedure Multprec_to_Penta_Double is

    f,df : Floating_Number;
    d : penta_double;

  begin
    new_line;
    put("Give a multiprecision float : "); get(f);
    d := to_penta_double(f);
    put_line("-> your multiprecision float : "); put(f); new_line;
    put_line(" the converted penta double : "); put(d); new_line;
    df := to_floating_number(d);
    put_line("-> converted back to multiprecision float : ");
    put(df); new_line;
  end Multprec_to_Penta_Double;

  procedure Penta_Double_to_Multprec is

    d,fd : penta_double;
    f : Floating_Number;

  begin
    new_line;
    put("Give a penta double : "); get(d);
    f := to_floating_number(d);
    put_line("-> your penta double : "); put(d); new_line;
    put_line(" the converted multiprecision float : "); put(f); new_line;
    fd := to_penta_double(f);
    put_line("-> converted back to penta double : "); put(fd); new_line;
  end Penta_Double_to_Multprec;

  procedure Complex_to_Penta_Double is
    
    c : Standard_Complex_Numbers.Complex_Number;
    dd_c : PentDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put("Give a complex number : "); get(c);
    dd_c := Standard_to_PentDobl_Complex(c);
    put_line("-> penta double : "); put(dd_c); new_line;
  end Complex_to_Penta_Double;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Conversions between penta double and multiprecision float ...");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line(" 0. exit the program;");
      put_line(" 1. multprecision -> penta double -> multiprecision;");
      put_line(" 2. penta double -> multprecision -> penta double;");
      put_line(" 3. convert complex numbers.");
      put("Type 0, 1, 2, or 3 : "); Ask_Alternative(ans,"0123");
      exit when ans = '0';
      case ans is
        when '1' => Multprec_to_Penta_Double;
        when '2' => Penta_Double_to_Multprec;
        when '3' => Complex_to_Penta_Double;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Multprec_PentDobl_Casts;
