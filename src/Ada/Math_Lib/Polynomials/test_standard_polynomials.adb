with text_io;                            use text_io;
with Symbol_Table;
with Symbol_Table_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Vectors;      use Standard_Complex_Poly_Vectors;
with Standard_Complex_Poly_Vectors_io;   use Standard_Complex_Poly_Vectors_io;
with Standard_Complex_Poly_Matrices;     use Standard_Complex_Poly_Matrices;
with Standard_Complex_Poly_Matrices_io;  use Standard_Complex_Poly_Matrices_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;

package body Test_Standard_Polynomials is

  procedure Test_Standard_io is

    n,m : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;

  begin
    new_line;
    put_line
      ("Interactive testing of input/output of standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    m := Number_of_Unknowns(p);
    put("the number of unknowns : "); put(m,1); new_line;
    put("the number of terms : "); put(Number_of_Terms(p),1); new_line;
    put("the size of the support : "); put(Size_of_Support(p),1); new_line;
    put("the degree of p : "); put(Degree(p),1);
    put("  max degrees : ");
    for i in 1..integer32(m) loop
      put(Degree(p,i),1); put(" ");
    end loop;
    new_line;
    put_line("Your polynomial : "); put(p); new_line;
    Symbol_Table.Clear;
    Clear(p);
  end Test_Standard_io;

  procedure Test_Vector_io is

    n : natural32 := 0;

  begin
    new_line;
    put_line("Interactive testing of i/o of vectors of complex polynomials.");
    new_line;
    put("Give the dimension of the vector : "); get(n);
    Symbol_Table.Init(n);
    declare
      p : Standard_Complex_Poly_Vectors.Vector(1..integer32(n));
    begin
      put("Give "); put(n,1); put(" polynomials in "); put(n,1);
      put_line(" variables : "); get(p);
      put_line("Your polynomials : "); put_line(p);
    end;
  end Test_Vector_io;

  procedure Test_Matrix_io is

    n : integer32 := 0;

  begin
    new_line;
    put_line("Interactive testing of i/o of matrices of complex polynomials.");
    new_line;
    put("Give the dimension of the matrix : "); get(n);
    Symbol_Table.Init(natural32(n));
    declare
      p : Matrix(1..n,1..n);
    begin
      put("Give "); put(n,1); put("x"); put(n,1);
      put(" polynomial matrix in "); put(n,1);
      put_line(" variables : "); get(p);
      put_line("Your polynomial matrix : "); put(p);
    end;
  end Test_Matrix_io;

  procedure Test_Standard_Eval
              ( p : in Standard_Complex_Polynomials.Poly;
                e : in Standard_Complex_Poly_Functions.Eval_Poly;
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
  end Test_Standard_Eval;

  procedure Interactive_Standard_Eval is

    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    e : Standard_Complex_Poly_Functions.Eval_Poly;
    bug : boolean;
    ans : character;

  begin
    new_line;
    put_line("Interactive evaluation of standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    put_line("Your polynomial p : "); put(p); new_line;
    e := Create(p);
    loop
      declare
        x : Standard_Complex_Vectors.Vector(1..integer32(n));
      begin
        put("Give "); put(n,1); put_line(" complex numbers : "); get(x);
        Test_Standard_Eval(p,e,x,true,bug);
      end;
      put("Do you wish to evaluate at other points (y/n) ? "); get(ans);
      exit when (ans /= 'y');
    end loop;
    Symbol_Table.Clear;
    Clear(p); Clear(e);
  end Interactive_Standard_Eval;

  procedure Random_Standard_Eval is

    n,nb : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    e : Standard_Complex_Poly_Functions.Eval_Poly;
    bug : boolean;

  begin
    new_line;
    put_line("Random evaluation of standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    put_line("Your polynomial p : "); put(p); new_line;
    e := Create(p);
    put("Give the number of samples : "); get(nb);
    for i in 1..nb loop
      declare
        x : constant Standard_Complex_Vectors.Vector
          := Random_Vector(1,integer32(n));
      begin
        Test_Standard_Eval(p,e,x,false,bug);
      end;
      exit when bug;
    end loop;
    Symbol_Table.Clear;
    Clear(p); Clear(e);
  end Random_Standard_Eval;

  procedure Test_Standard_Diff is

    n,m : natural32 := 0;
    p,dp : Standard_Complex_Polynomials.Poly;

  begin
    new_line;
    put_line("Test on differentiation of standard complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    m := Number_of_Unknowns(p);
    put("Your polynomial p : "); put(p); new_line;
    put("The number of unknowns : "); put(m,1); new_line;
    for i in 1..integer32(m) loop
      dp := Diff(p,i);
      put("Diff(p,"); put(i,1); put(") : ");
      put(dp); new_line;
      Clear(dp);
    end loop;
    Symbol_Table.Clear;
    Clear(p);
  end Test_Standard_Diff;

  procedure Test_System_io is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Interactive testing on input/output of polynomial systems.");
    new_line;
    get(lp);
    put_line("Your system : "); put(lp.all);
    new_line;
    put_line("The symbols read :");
    for i in 1..Symbol_Table.Number loop
      put(" "); Symbol_Table_io.put(Symbol_Table.get(i));
    end loop;
    new_line;
  end Test_System_io;

  procedure Test_Eval_Standard_System is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the evaluation of standard polynomial systems.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    declare
      n : constant integer32 := lp'last;
      x,y1,y2,y3 : Standard_Complex_Vectors.Vector(1..n);
      ep : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(1..n) := Create(lp.all);

      function Evaluate ( x : Standard_Complex_Vectors.Vector )
                        return Standard_Complex_Vectors.Vector is
      begin
        return Eval(ep,x);
      end Evaluate;
    begin
      put("Give "); put(n,1); put_line(" complex numbers :"); get(x);
      y1 := Eval(lp.all,x);
      y2 := Eval(ep,x);
      y3 := Evaluate(x);
      put("p(x) : "); put(y1); new_line;
      put("e(x) : "); put(y2); new_line;
      put("E(x) : "); put(y2); new_line;
      if Equal(y1,y2) and Equal(y2,y3)
       then put_line("Test on evaluation of system is successful.");
       else put_line("Different evaluations!  Bug detected.");
      end if;
      Clear(ep);
    end;
  end Test_Eval_Standard_System;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing the operations on complex polynomials ...");
    loop
      new_line;
      put_line("Choose one of the following :                               ");
      put_line("  0. Exit this program.                                     ");
      put_line("  1. i/o of standard complex polynomials.                   ");
      put_line("  2. i/o of vectors of standard complex polynomials.        ");
      put_line("  3. i/o of matrices of standard complex polynomials.       ");
      put_line("  4. Interactive evaluation of standard complex polynomials.");
      put_line("  5. Random evaluation of standard complex polynomials.     ");
      put_line("  6. Differentiation of standard complex polynomials.       ");
      put_line("  7. i/o of systems of standard complex polynomials.        ");
      put_line("  8. Evaluation of systems of standard complex polynomials. ");
      put("Type 0, 1, 2, 3, 4, 5, 6, 7, or 8 : ");
      get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Test_Standard_io;
        when '2' => Test_Vector_io;
        when '3' => Test_Matrix_io;
        when '4' => Interactive_Standard_Eval;
        when '5' => Random_Standard_Eval;
        when '6' => Test_Standard_Diff;
        when '7' => Test_System_io;
        when '8' => Test_Eval_Standard_System;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Standard_Polynomials;
