with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Deca_Double_Numbers;               use Deca_Double_Numbers;
with Deca_Double_Numbers_io;            use Deca_Double_Numbers_io;
with Multprec_DecaDobl_Convertors;       use Multprec_DecaDobl_Convertors;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers_io;        use DecaDobl_Complex_Numbers_io;
with DecaDobl_Complex_Numbers_cv;        use DecaDobl_Complex_Numbers_cv;

package body Test_Multprec_DecaDobl_Casts is

  procedure Multprec_to_Deca_Double is

    f,df : Floating_Number;
    d : deca_double;

  begin
    new_line;
    put("Give a multiprecision float : "); get(f);
    d := to_deca_double(f);
    put_line("-> your multiprecision float : "); put(f); new_line;
    put_line(" the converted deca double : "); put(d); new_line;
    df := to_floating_number(d);
    put_line("-> converted back to multiprecision float : ");
    put(df); new_line;
  end Multprec_to_Deca_Double;

  procedure Deca_Double_to_Multprec is

    d,fd : deca_double;
    f : Floating_Number;

  begin
    new_line;
    put("Give a deca double : "); get(d);
    f := to_floating_number(d);
    put_line("-> your deca double : "); put(d); new_line;
    put_line(" the converted multiprecision float : "); put(f); new_line;
    fd := to_deca_double(f);
    put_line("-> converted back to deca double : "); put(fd); new_line;
  end Deca_Double_to_Multprec;

  procedure Complex_to_Deca_Double is
    
    c : Standard_Complex_Numbers.Complex_Number;
    dd_c : DecaDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put("Give a complex number : "); get(c);
    dd_c := Standard_to_DecaDobl_Complex(c);
    put_line("-> deca double : "); put(dd_c); new_line;
  end Complex_to_deca_double;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Conversions between deca double and multiprecision float ...");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line(" 0. exit the program;");
      put_line(" 1. multprecision -> deca double -> multiprecision;");
      put_line(" 2. deca double -> multprecision -> deca double;");
      put_line(" 3. convert complex numbers.");
      put("Type 0, 1, 2, or 3 : "); Ask_Alternative(ans,"0123");
      exit when ans = '0';
      case ans is
        when '1' => Multprec_to_Deca_double;
        when '2' => Deca_Double_to_Multprec;
        when '3' => Complex_to_Deca_double;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Multprec_DecaDobl_Casts;
