with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_to_Multprec_Convertors;    use Standard_to_Multprec_Convertors;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Multprec_Complex_Polynomials;       use Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;    use Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Systems;      use Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;       use Multprec_Complex_Poly_SysFun;

package body Test_Multprec_Polynomials is

  procedure Test_Multprec_io is

    n,m,size : natural32 := 0;
    p : Multprec_Complex_Polynomials.Poly;

  begin
    new_line;
    put_line
      ("Interactive testing of input/output of multprec complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    put("Give the size of the numbers : "); get(size);
    Symbol_Table.Init(n);
    Multprec_Complex_Polynomials_io.Set_Working_Precision(size);
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
  end Test_Multprec_io;

  procedure Test_Multprec_Eval is

    ans : character;
    n : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    mp : Multprec_Complex_Polynomials.Poly;
    ep : Multprec_Complex_Poly_Functions.Eval_Poly;
    ecp : Multprec_Complex_Poly_Functions.Eval_Coeff_Poly;

  begin
    new_line;
    put_line("Testing the evaluation of multiprecision complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    put_line("Your polynomial p : "); put(p); new_line;
    mp := Convert(p);
    ep := Create(mp);
    ecp := Create(mp);
    declare
      x : Multprec_Complex_Vectors.Vector(1..integer32(n));
      cp : constant Multprec_Complex_Vectors.Vector
         := Multprec_Complex_Poly_Functions.Coeff(mp); 
      sz : natural32 := 0;
      y1,y2,y3 : Multprec_Complex_Numbers.Complex_Number;
    begin
      put_line("The coefficient vector : "); put_line(cp);
      put("Give "); put(n,1); put_line(" complex numbers : "); get(x);
      loop
        put("Give the size of the numbers : "); get(sz);
        Set_Size(x,sz);
        y1 := Eval(mp,x);     put("p(x)   : "); put(y1); new_line;
        y2 := Eval(ep,x);     put("e(x)   : "); put(y2); new_line;
        y3 := Eval(ecp,cp,x); put("f(c,x) : "); put(y3); new_line;
        if Equal(y1,y2) and Equal(y1,y3)
         then put_line("Test on evaluation is successful.");
         else put_line("Different results!  Bug detected.");
        end if;
        put("Do you want to evaluate for other precisions ? (y/n) ");
        get(ans);
        exit when ans /= 'y';
      end loop;
    end;
    Symbol_Table.Clear;
    Clear(p); Clear(mp); Clear(ep);
  end Test_Multprec_Eval;

  procedure Test_Multprec_Diff is

    n,m : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    mp,dp : Multprec_Complex_Polynomials.Poly;

  begin
    new_line;
    put_line("Test on differentiation of multiprecision complex polynomials.");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Give a polynomial (terminate with ;) : "); get(p);
    m := Number_of_Unknowns(p);
    put("Your polynomial p : "); put(p); new_line;
    mp := Convert(p);
    put("As multiprecision poly : "); put(mp); new_line;
    put("The number of unknowns : "); put(m,1); new_line;
    for i in 1..integer32(m) loop
      dp := Diff(mp,i);
      put("Diff(p,"); put(i,1); put(") : ");
      put(dp); new_line;
      Clear(dp);
    end loop;
    Symbol_Table.Clear;
    Clear(p); Clear(mp);
  end Test_Multprec_Diff;

  procedure Test_Eval_Multprec_System is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the evaluation of multiprecision polynomial systems.");
    new_line;
    get(lp);
    put_line("The system : "); put(lp.all);
    declare
      n : constant integer32 := lp'last;
      x,y1,y2 : Multprec_Complex_Vectors.Vector(1..n);
      mp : Multprec_Complex_Poly_Systems.Poly_Sys(1..n) := Convert(lp.all);
      ep : Multprec_Complex_Poly_SysFun.Eval_Poly_Sys(1..n) := Create(mp);
    begin
      put("Give "); put(n,1); put_line(" complex numbers :"); get(x);
      y1 := Eval(mp,x);
      y2 := Eval(ep,x);
      put("p(x) : "); put(y1); new_line;
      put("e(x) : "); put(y2); new_line;
      if Equal(y1,y2)
       then put_line("Test on evaluation of system is successful.");
       else put_line("Different evaluations!  Bug detected.");
      end if;
      Clear(mp); Clear(ep);
    end;
  end Test_Eval_Multprec_System;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing multiprecision complex polynomials ...");
    loop
      new_line;
      put_line("Choose one of the following :                            ");
      put_line("  0. Exit this program.                                  ");
      put_line("  1. i/o of multiprecision complex polynomials.         ");
      put_line("  2. Evaluation of multiprecision complex polynomials.  ");
      put_line("  3. Differentiation of multiprecision polynomials.      ");
      put_line("  4. Evaluation of systems of multiprecision polynomials.");
      put("Type 0, 1, 2, 3, 4, 5, or 6 to select : ");
      Ask_Alternative(ans,"0123456");
      exit when (ans = '0');
      case ans is
        when '1' => Test_Multprec_io;
        when '2' => Test_Multprec_Eval;
        when '3' => Test_Multprec_Diff;
        when '4' => Test_Eval_Multprec_System;
        when others => null;
      end case;
    end loop;
  end Main;

end Test_Multprec_Polynomials;
