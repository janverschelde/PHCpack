with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Octo_Double_Numbers_io;            use Octo_Double_Numbers_io;
with Multprec_OctoDobl_Convertors;       use Multprec_OctoDobl_Convertors;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;        use OctoDobl_Complex_Numbers_io;
with OctoDobl_Complex_Numbers_cv;        use OctoDobl_Complex_Numbers_cv;

package body Test_Multprec_OctoDobl_Casts is

  procedure Multprec_to_Octo_Double is

    f,df : Floating_Number;
    d : octo_double;

  begin
    new_line;
    put("Give a multiprecision float : "); get(f);
    d := to_octo_double(f);
    put_line("-> your multiprecision float : "); put(f); new_line;
    put_line(" the converted octo double : "); put(d); new_line;
    df := to_floating_number(d);
    put_line("-> converted back to multiprecision float : ");
    put(df); new_line;
  end Multprec_to_Octo_Double;

  procedure Octo_Double_to_Multprec is

    d,fd : octo_double;
    f : Floating_Number;

  begin
    new_line;
    put("Give a octo double : "); get(d);
    f := to_floating_number(d);
    put_line("-> your octo double : "); put(d); new_line;
    put_line(" the converted multiprecision float : "); put(f); new_line;
    fd := to_octo_double(f);
    put_line("-> converted back to octo double : "); put(fd); new_line;
  end Octo_Double_to_Multprec;

  procedure Complex_to_Octo_Double is
    
    c : Standard_Complex_Numbers.Complex_Number;
    dd_c : OctoDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put("Give a complex number : "); get(c);
    dd_c := Standard_to_OctoDobl_Complex(c);
    put_line("-> octo double : "); put(dd_c); new_line;
  end Complex_to_Octo_Double;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Conversions between octo double and multiprecision float ...");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line(" 0. exit the program;");
      put_line(" 1. multprecision -> octo double -> multiprecision;");
      put_line(" 2. octo double -> multprecision -> octo double;");
      put_line(" 3. convert complex numbers.");
      put("Type 0, 1, 2, or 3 : "); Ask_Alternative(ans,"0123");
      exit when ans = '0';
      case ans is
        when '1' => Multprec_to_Octo_Double;
        when '2' => Octo_Double_to_Multprec;
        when '3' => Complex_to_Octo_Double;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Multprec_OctoDobl_Casts;
