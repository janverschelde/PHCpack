with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_QuadDobl_Convertors;       use Multprec_QuadDobl_Convertors;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_cv;        use QuadDobl_Complex_Numbers_cv;

package body Test_Multprec_QuadDobl_Casts is

  procedure Multprec_to_Quad_Double is

    f,df : Floating_Number;
    d : quad_double;

  begin
    new_line;
    put("Give a multiprecision float : "); get(f);
    d := to_quad_double(f);
    put_line("-> your multiprecision float : "); put(f); new_line;
    put_line(" the converted quad double : "); put(d); new_line;
    df := to_floating_number(d);
    put_line("-> converted back to multiprecision float : ");
    put(df); new_line;
  end Multprec_to_Quad_Double;

  procedure Quad_Double_to_Multprec is

    d,fd : quad_double;
    f : Floating_Number;

  begin
    new_line;
    put("Give a quad double : "); get(d);
    f := to_floating_number(d);
    put_line("-> your quad double : "); put(d); new_line;
    put_line(" the converted multiprecision float : "); put(f); new_line;
    fd := to_quad_double(f);
    put_line("-> converted back to quad double : "); put(fd); new_line;
  end Quad_Double_to_Multprec;

  procedure Complex_to_Quad_Double is

    c : Standard_Complex_Numbers.Complex_Number;
    qd_c : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put("Give a complex number : "); get(c);
    qd_c := Standard_to_QuadDobl_Complex(c);
    put_line("-> quad double : "); put(qd_c); new_line;
  end Complex_to_Quad_Double;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Conversions between quad double and multiprecision float ...");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line(" 0. exit the program;");
      put_line(" 1. multprecision -> quad double -> multiprecision;");
      put_line(" 2. quad double -> multprecision -> quad double;");
      put_line(" 3. convert complex numbers.");
      put("Type 0, 1, 2, or 3 : "); Ask_Alternative(ans,"0123");
      exit when ans = '0';
      case ans is
        when '1' => Multprec_to_Quad_Double;
        when '2' => Quad_Double_to_Multprec;
        when '3' => Complex_to_Quad_Double;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Multprec_QuadDobl_Casts;
