with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Chebychev_Polynomials;              use Chebychev_Polynomials;

procedure ts_cheby is

  procedure Test_Create is

    k : natural32 := 0;

  begin
    new_line;
    put_line("Testing the creation of Chebychev polynomials.");
    new_line;
    put("Give the degree k : "); get(k);
    for i in 0..k loop
      declare
        pi : constant Vector := Create(i);
      begin
        put("The coefficients of the "); put(i,1);
        put_line("-degree Chebychev polynomial : ");
        put_line(pi);
      end;
    end loop;
  end Test_Create;

  procedure Test_Diff is

    k : natural32 := 0;
    ans : character;

  begin
    new_line;   
    put_line("Testing the differentiation of Chebychev polynomials.");
    loop
      new_line;
      put("Give the degree k : "); get(k);
      declare
        pk : constant Vector := Create(k);
        dp : constant Vector := Diff(pk);
      begin
        put("The coefficients of the "); put(k,1);
        put_line("-degree Chebychev polynomial : ");
        put_line(pk);
        put_line("The coefficients of its derivative : ");
        put_line(dp);
        for i in 1..k loop
          declare
            dpi : constant Vector := Diff(pk,i);
          begin
            put("The coefficients after deriving "); put(i,1);
            put_line(" times : "); put_line(dpi);
          end;
        end loop;
      end;
      put("Do you want more tests (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Diff;

  procedure Test_Int is

    k : natural32 := 0;
    ans : character;

  begin
    new_line;   
    put_line("Testing the Antidifferentiation of Chebychev polynomials.");
    loop
      new_line;
      put("Give the degree k : "); get(k);
      declare
        pk : constant Vector := Create(k);
        dp : constant Vector := Int(pk);
      begin
        put("The coefficients of the "); put(k,1);
        put_line("-degree Chebychev polynomial : ");
        put_line(pk);
        put_line("The coefficients of its antiderivative : ");
        put_line(dp);
        for i in 1..k loop
          declare
            dpi : constant Vector := Int(pk,i);
          begin
            put("The coefficients after antideriving "); put(i,1);
            put_line(" times : "); put_line(dpi);
          end;
        end loop;
      end;
      put("Do you want more tests (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Int;

  procedure Test_Eval is

    k : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing the evaluation of Chebychev polynomials.");
    loop
      new_line;
      put("Give the degree k : "); get(k);
      declare
        p : constant Vector := Create(k);
        x : double_float := 0.0;
      begin
        loop
          put("Give x : "); get(x);
          put("p(x) : "); put(Eval(p,x)); new_line;
          put("COS-ARCCOS-eval : "); put(Eval(k,x)); new_line;
          put("Do you want more evaluations ? (y/n) "); get(ans);
          exit when (ans /= 'y');
        end loop;
      end;
      put("Do you want other polynomials to evaluate ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Eval;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing the manipulation of Chebychev polynomials.");
    loop
      new_line;
      put_line("Choose one of the following :");
      put_line("  0. Exit this program.");
      put_line("  1. Creation of Chebychev polynomials.");
      put_line("  2. Differentiation of Chebychev polynomials.");
      put_line("  3. Antidifferentiation of Chebychev polynomials.");
      put_line("  4. Evaluation of Chebychev polynomials.");
      put("Type 0,1,2,3 or 4 to make your choice : "); get(ans);
      exit when (ans = '0');
      case ans is 
        when '1' => Test_Create;
        when '2' => Test_Diff;
        when '3' => Test_Int;
        when '4' => Test_Eval;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_cheby;
