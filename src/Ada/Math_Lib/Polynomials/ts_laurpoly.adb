with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Tools;
with Multprec_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laurentials_io;    use Multprec_Complex_Laurentials_io;
with Multprec_Complex_Laur_Functions;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_Systems_io;   use Multprec_Complex_Laur_Systems_io;

procedure ts_laurpoly is

-- DESCRIPTION :
--   Interactive testing routines on working with Laurent polynomials.

  procedure Read ( n : out natural32;
                   p : out Standard_Complex_Laurentials.Poly ) is

  begin
    n := 0;
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Reading a polynomial, terminate with semicolon ; ");
    get(p);
  end Read;

  procedure Read ( n : out natural32;
                   p : out Multprec_Complex_Laurentials.Poly ) is

    size : natural32 := 0;

  begin
    n := 0;
    new_line;
    put("Give the number of variables : "); get(n);
    put("Give the size of the numbers : "); get(size);
    Symbol_Table.Init(n);
    Set_Working_Precision(size);
    put_line("Reading a polynomial, terminate with semicolon ; ");
    get(p);
  end Read;

  procedure Standard_Test_Read_Write is

  -- DESCRIPTION :
  --   Does a test on reading and writing of a Laurent polynomial
  --   with standard complex coefficients.

    n : natural32 := 0;
    p : Standard_Complex_Laurentials.Poly;

  begin
    Read(n,p);
    put_line("Your polynomial is "); put(p); new_line;
    put_line("Your polynomial is ... "); put_line(p);
  end Standard_Test_Read_Write;

  procedure Standard_Test_Read_Write_System is

  -- DESCRIPTION :
  --   Does a test on reading and writing of a Laurent polynomial system
  --   with standard complex coefficients.

    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    new_line;
    get(lp);
    put_line("Your Laurent polynomial system is ");
    put(lp.all);
    put_line("Your Laurent polynomial system is ");
    put_line(lp.all);
  end Standard_Test_Read_Write_System;

  procedure Multprec_Test_Read_Write is

  -- DESCRIPTION :
  --   Does a test on reading and writing of a Laurent polynomial
  --   with multiprecision complex coefficients.

    n : natural32 := 0;
    p : Multprec_Complex_Laurentials.Poly;

  begin
    Read(n,p);
    put_line("Your polynomial is "); put(p); new_line;
    put_line("Your polynomial is ... "); put_line(p);
  end Multprec_Test_Read_Write;

  procedure Multprec_Test_Read_Write_System is

  -- DESCRIPTION :
  --   Does a test on reading and writing of a Laurent polynomial system
  --   with multiprecision complex coefficients.

    lp : Multprec_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    new_line;
    get(lp);
    put_line("Your Laurent polynomial system is ");
    put(lp.all);
    put_line("Your Laurent polynomial system is ");
    put_line(lp.all);
  end Multprec_Test_Read_Write_System;

  procedure Test_Evaluation
               ( n : in natural32;
                 p : in Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Functions;

    x : Standard_Complex_Vectors.Vector(1..integer32(n));
    y1,y2,y3,y4 : Complex_Number;
    ep : constant Eval_Poly := Create(p);
    ecp : constant Eval_Coeff_Poly := Create(p);
    cv : constant Standard_Complex_Vectors.Vector := Coeff(p);
    ans : character;

  begin
    put_line("The coefficient vector : "); put_line(cv);
    loop
      put("Evaluate at random or user given point ? (r/u) ");
      Ask_Alternative(ans,"ru");
      if ans = 'u' then
        put("Give "); put(n,1); put_line(" complex numbers : "); get(x);
      else
        x := Standard_Random_Vectors.Random_Vector(1,integer32(n));
      end if;
      y1 := Eval(p,x);
      put_line("At the vector : "); put_line(x);
      put("the value is "); put(y1); new_line;
      y2 := Eval(ep,x);
      put("   and again "); put(y2); new_line;
      y3 := Eval(p,cv,x);
      put("   and again "); put(y3); new_line;
      y4 := Eval(ecp,cv,x);
      put("   and again "); put(y4); new_line;
      put("Evaluate at more points ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Evaluation;

  procedure Test_Evaluation
               ( n : in natural32;
                 p : in Multprec_Complex_Laurentials.Poly ) is

    use Multprec_Complex_Numbers;
    use Multprec_Complex_Laurentials;
    use Multprec_Complex_Laur_Functions;

    x : Multprec_Complex_Vectors.Vector(1..integer32(n));
    y1,y2,y3,y4 : Complex_Number;
    ep : constant Eval_Poly := Create(p);
    ecp : constant Eval_Coeff_Poly := Create(p);
    cv : constant Multprec_Complex_Vectors.Vector := Coeff(p);
    ans : character;
    size : natural32 := 0;

  begin
    put_line("The coefficient vector : "); put_line(cv);
    put("Give the size of the arguments : "); get(size);
    loop
      put("Evaluate at random or user given point ? (r/u) ");
      Ask_Alternative(ans,"ru");
      if ans = 'u' then
        put("Give "); put(n,1); put_line(" complex numbers : "); get(x);
        Multprec_Complex_Vector_Tools.Set_Size(x,size);
      else
        x := Multprec_Random_Vectors.Random_Vector(1,integer32(n),size);
      end if;
      y1 := Eval(p,x);
      put_line("At the vector : "); put_line(x);
      put("the value is "); put(y1); new_line;
      y2 := Eval(ep,x);
      put("   and again "); put(y2); new_line;
      y3 := Eval(p,cv,x);
      put("   and again "); put(y3); new_line;
      y4 := Eval(ecp,cv,x);
      put("   and again "); put(y4); new_line;
      put("Evaluate at more points ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_Evaluation;

  procedure Standard_Test_Evaluation is

  -- DESCRIPTION :
  --   Test on evaluation of Laurent polynomial
  --   with standard complex coefficients.

    n : natural32 := 0;
    p : Standard_Complex_Laurentials.Poly;

  begin
    Read(n,p);
    put("Your polynomial : "); put(p); new_line;
    Test_Evaluation(n,p);
  end Standard_Test_Evaluation;

  procedure Multprec_Test_Evaluation is

  -- DESCRIPTION :
  --   Test on evaluation of Laurent polynomial
  --   with standard complex coefficients.

    n : natural32 := 0;
    p : Multprec_Complex_Laurentials.Poly;

  begin
    Read(n,p);
    put("Your polynomial : "); put(p); new_line;
    Test_Evaluation(n,p);
  end Multprec_Test_Evaluation;

  procedure Standard_Test_Differentiation is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial and computes
  --   all partial derivatives.

    n : natural32 := 0;
    p : Standard_Complex_Laurentials.Poly;

  begin
    Read(n,p);
    put("Your polynomial : "); put(p); new_line;
    for i in 1..integer32(n) loop
      put("Derivative w.r.t. variable "); put(i,1); put_line(" :");
      declare
        dp : constant Standard_Complex_Laurentials.Poly
           := Standard_Complex_Laurentials.Diff(p,i); 
      begin
        put(dp); new_line;
      end;
    end loop;
  end Standard_Test_Differentiation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU for reading and writing Laurent polynomials... ");
    put_line("  1. read/write of one single Laurent polynomial;");
    put_line("  2. read/write of a Laurent polynomial system;");
    put_line("  3. test evaluation of a standard Laurent polynomial;");
    put_line("  4. test differentiation of a Laurent polynomial;");
    put_line("  5. test i/o for multiprecision Laurent polynomial;");
    put_line("  6. test i/o for multiprecision Laurent system;");
    put_line("  7. test evaluation of a multprecision Laurent polynomial.");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to make your choice : ");
    Ask_Alternative(ans,"1234567");
    case ans is
      when '1' => Standard_Test_Read_Write;
      when '2' => Standard_Test_Read_Write_System;
      when '3' => Standard_Test_Evaluation;
      when '4' => Standard_Test_Differentiation;
      when '5' => Multprec_Test_Read_Write;
      when '6' => Multprec_Test_Read_Write_System;
      when '7' => Multprec_Test_Evaluation;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_laurpoly;
