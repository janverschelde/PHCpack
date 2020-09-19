with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;

package body Test_Standard_Laurentials is

  procedure Test_Standard_Laurent_Eval
              ( p : in Standard_Complex_Laurentials.Poly;
                e : in Standard_Complex_Laur_Functions.Eval_Poly;
                x : in Standard_Complex_Vectors.Vector;
                output_of_results : in boolean; bug : out boolean ) is

    y1 : constant Standard_Complex_Numbers.Complex_Number := Eval(p,x);
    y2 : constant Standard_Complex_Numbers.Complex_Number := Eval(e,x);
    tol : constant double_float := 1.0E-12;

  begin
    if AbsVal(y1-y2) < tol
     then put_line("Test on evaluation is successful."); bug := false;
     else put_line("Different results!  Bug detected."); bug := true;
    end if;
    if output_of_results or bug then
      put("p(x) : "); put(y1); new_line;
      put("e(x) : "); put(y2); new_line;
    end if;
  end Test_Standard_Laurent_Eval;

  procedure Interactive_Standard_Laurent_Eval is

    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    lp : Standard_Complex_Laurentials.Poly;
    elp : Standard_Complex_Laur_Functions.Eval_Poly;
    bug : boolean;
    ans : character;

  begin
    new_line;
    put_line("Interactive evaluation of standard complex Laurent polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    put_line("Your polynomial p : "); put(p); new_line;
    lp := Polynomial_to_Laurent_Polynomial(p);
    elp := Create(lp);
    loop
      declare
        x : Standard_Complex_Vectors.Vector(1..integer32(n));
      begin
        put("Give "); put(n,1); put_line(" complex numbers : "); get(x);
        Test_Standard_Laurent_Eval(lp,elp,x,true,bug);
      end;
      put("Do you wish to evaluate at other points (y/n) ? "); get(ans);
      exit when (ans /= 'y');
    end loop;
    Symbol_Table.Clear;
    Clear(p); Clear(lp); Clear(elp);
  end Interactive_Standard_Laurent_Eval;

  procedure Random_Standard_Laurent_Eval is

    n,nb : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    lp : Standard_Complex_Laurentials.Poly;
    elp : Standard_Complex_Laur_Functions.Eval_Poly;
    bug : boolean;

  begin
    new_line;
    put_line("Random evaluation of standard complex Laurent polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    put_line("Your polynomial p : "); put(p); new_line;
    lp := Polynomial_to_Laurent_Polynomial(p);
    elp := Create(lp);
    put("Give the number of samples : "); get(nb);
    for i in 1..nb loop
      declare
        x : constant Standard_Complex_Vectors.Vector 
          := Random_Vector(1,integer32(n));
      begin
        Test_Standard_Laurent_Eval(lp,elp,x,false,bug);
      end;
      exit when bug;
    end loop;
    Symbol_Table.Clear;
    Clear(p); Clear(lp); Clear(elp);
  end Random_Standard_Laurent_Eval;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test complex Laurent polynomials.");
    loop
      new_line;
      put_line("Choose one of the following :                            ");
      put_line("  0. Exit this program.                                  ");
      put_line("  1. evaluation of standard Laurent polynomials.         ");
      put_line("  2. Random evaluation of standard Laurent polynomials.  ");
      put("Type 0, 1, or 2 to select : ");
      Ask_Alternative(ans,"012");
      exit when (ans = '0');
      case ans is
        when '1' => Interactive_Standard_Laurent_Eval;
        when '2' => Random_Standard_Laurent_Eval;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Standard_Laurentials;
