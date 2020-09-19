with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Numbers_io;           use Triple_Double_Numbers_io;
with Multprec_TripDobl_Convertors;       use Multprec_TripDobl_Convertors;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with TripDobl_Complex_Numbers;
with TripDobl_Complex_Numbers_io;        use TripDobl_Complex_Numbers_io;
with TripDobl_Complex_Numbers_cv;        use TripDobl_Complex_Numbers_cv;

package body Test_Multprec_TripDobl_Casts is

  procedure Multprec_to_Triple_Double is

    f,df : Floating_Number;
    d : triple_double;

  begin
    new_line;
    put("Give a multiprecision float : "); get(f);
    d := to_triple_double(f);
    put_line("-> your multiprecision float : "); put(f); new_line;
    put_line(" the converted triple double : "); put(d); new_line;
    df := to_floating_number(d);
    put_line("-> converted back to multiprecision float : ");
    put(df); new_line;
  end Multprec_to_Triple_Double;

  procedure Triple_Double_to_Multprec is

    d,fd : triple_double;
    f : Floating_Number;

  begin
    new_line;
    put("Give a triple double : "); get(d);
    f := to_floating_number(d);
    put_line("-> your triple double : "); put(d); new_line;
    put_line(" the converted multiprecision float : "); put(f); new_line;
    fd := to_triple_double(f);
    put_line("-> converted back to triple double : "); put(fd); new_line;
  end Triple_Double_to_Multprec;

  procedure Complex_to_Triple_Double is
    
    c : Standard_Complex_Numbers.Complex_Number;
    dd_c : TripDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put("Give a complex number : "); get(c);
    dd_c := Standard_to_TripDobl_Complex(c);
    put_line("-> triple double : "); put(dd_c); new_line;
  end Complex_to_Triple_Double;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Conversions between triple double and multiprecision float ...");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line(" 0. exit the program;");
      put_line(" 1. multprecision -> triple double -> multiprecision;");
      put_line(" 2. triple double -> multprecision -> triple double;");
      put_line(" 3. convert complex numbers.");
      put("Type 0, 1, 2, or 3 : "); Ask_Alternative(ans,"0123");
      exit when ans = '0';
      case ans is
        when '1' => Multprec_to_Triple_Double;
        when '2' => Triple_Double_to_Multprec;
        when '3' => Complex_to_Triple_Double;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Multprec_TripDobl_Casts;
