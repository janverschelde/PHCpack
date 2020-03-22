with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Complex_Series_Vectors;
with Standard_Random_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_CSeries_Vector_Functions;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Random_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_CSeries_Vector_Functions;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Random_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with QuadDobl_CSeries_Vector_Functions;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_SysFun;
with Complex_Series_and_Polynomials_io;  use Complex_Series_and_Polynomials_io;
with Series_and_Homotopies;

procedure ts_serhom is

-- DESCRIPTION :
--   Test on the operations in the package series_and_homotopies.

  procedure Standard_Test_Creation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Test on converting Standard_Homotopy.Homotopy_System of
  --   range 1..nq to a series system.
  --   The output of put_line(h) is the same as put(s,1),
  --   except for that the last variable in h has become
  --   the first variable in s.

    h : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq);

  begin
    put_line("The homotopy system :"); put_line(h);
    s := Series_and_Homotopies.Create(h,nq+1,true);
    put_line("The series system :"); put(s,1);
  end Standard_Test_Creation;

  procedure DoblDobl_Test_Creation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Test on converting DoblDobl_Homotopy.Homotopy_System of
  --   range 1..nq to a series system.
  --   The output of put_line(h) is the same as put(s,1),
  --   except for that the last variable in h has become
  --   the first variable in s.

    h : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);

  begin
    put_line("The homotopy system :"); put_line(h);
    s := Series_and_Homotopies.Create(h,nq+1,true);
    put_line("The series system :"); put(s,1);
  end DoblDobl_Test_Creation;

  procedure QuadDobl_Test_Creation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Test on converting QuadDobl_Homotopy.Homotopy_System of
  --   range 1..nq to a series system.
  --   The output of put_line(h) is the same as put(s,1),
  --   except for that the last variable in h has become
  --   the first variable in s.

    h : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);

  begin
    put_line("The homotopy system :"); put_line(h);
    s := Series_and_Homotopies.Create(h,nq+1,true);
    put_line("The series system :"); put(s,1);
  end QuadDobl_Test_Creation;

  procedure Standard_Test_Evaluation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in Standard_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the evaluated
  --   polynomial system is shown.

    p : Standard_Complex_Poly_Systems.Poly_Sys(1..nq);
    h : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : constant Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
    t : double_float := 0.0;
    ans : character;

  begin
    loop
      put("Give a value for t : "); get(t);
      p := Series_and_Homotopies.Eval(s,t);
      put_line("The evaluated system :"); put_line(p);
      put("Continue for other t values ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Test_Evaluation;

  procedure DoblDobl_Test_Evaluation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in DoblDobl_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the evaluated
  --   polynomial system is shown.

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    h : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
    t : double_double := create(0.0);
    ans : character;

  begin
    loop
      put("Give a value for t : "); get(t);
      p := Series_and_Homotopies.Eval(s,t);
      put_line("The evaluated system :"); put_line(p);
      put("Continue for other t values ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end DoblDobl_Test_Evaluation;

  procedure QuadDobl_Test_Evaluation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in QuadDobl_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the evaluated
  --   polynomial system is shown.

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    h : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
    t : quad_double := create(0.0);
    ans : character;

  begin
    loop
      put("Give a value for t : "); get(t);
      p := Series_and_Homotopies.Eval(s,t);
      put_line("The evaluated system :"); put_line(p);
      put("Continue for other t values ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end QuadDobl_Test_Evaluation;

  procedure Standard_Shift_Eval
              ( hom : in Standard_CSeries_Poly_Systems.Poly_Sys;
                sfh : in Standard_CSeries_Poly_Systems.Poly_Sys;
                c : in double_float; d : in integer32 ) is

  -- DESCRIPTION :
  --   Evaluates the polynomial system at a random series
  --   with respect to the shift value in c.

  -- ON ENTRY :
  --   hom      square homotopy;
  --   sfh      shifted homotopy with respected to c;
  --   c        value of the shift;
  --   d        degree of the series used for evaluation.

  -- REQUIRED :
  --   The homotopy in hom has as many equations as unknowns.

    n : constant integer32 := hom'last;
    x : constant Standard_Complex_Series_Vectors.Vector
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    cffhom : constant Standard_Complex_Series_VecVecs.VecVec(hom'range)
           := Standard_CSeries_Poly_SysFun.Coeff(hom);
    cshift : constant Standard_Complex_Series_VecVecs.VecVec(hom'range)
           := Standard_CSeries_Vector_Functions.Shift(cffhom,c);
    f : constant Standard_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(hom'range)
      := Standard_CSeries_Poly_SysFun.Create(hom);
    y : constant Standard_Complex_Series_Vectors.Vector
      := Standard_CSeries_Poly_SysFun.Eval(hom,x);
    sy : constant Standard_Complex_Series_Vectors.Vector
       := Standard_CSeries_Poly_SysFun.Eval(sfh,x);
    z : constant Standard_Complex_Series_Vectors.Vector
      := Standard_CSeries_Poly_SysFun.Eval(f,cffhom,x);
    sz : constant Standard_Complex_Series_Vectors.Vector
       := Standard_CSeries_Poly_SysFun.Eval(f,cshift,x);
    ans : character;

  begin
    Symbol_Table.Init(1);
    Symbol_Table.Add_String("t");
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Evaluations after the shift ...");
      put_line("The value at the polynomial :"); put(sy);
      put_line("The value at the coefficient polynomial :"); put(sz);
    end if;
  end Standard_Shift_Eval;

  procedure Standard_Test_Shift ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in Standard_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the system is shifted.

    h : constant Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : constant Standard_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
    c : double_float := 0.0;
    shifteds : Standard_CSeries_Poly_Systems.Poly_Sys(1..nq);
    y,z : Standard_Complex_Poly_Systems.Poly_Sys(1..nq);
    ans : character;
    deg : integer32 := 0;

  begin
    put("Give a value for the shift : "); get(c);
    shifteds := Series_and_Homotopies.Shift(s,c);
    y := Series_and_Homotopies.Eval(s,c);
    z := Series_and_Homotopies.Eval(shifteds,0.0);
    put_line("system(shift value) :"); put_line(y);
    put_line("shifted system(0.0) :"); put_line(z);
    put("Continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the degree of the series in x : "); get(deg);
      Standard_Shift_Eval(s,shifteds,c,deg);
    end if;
  end Standard_Test_Shift;

  procedure DoblDobl_Shift_Eval
              ( hom : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                sfh : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
                c : in double_double; d : in integer32 ) is

  -- DESCRIPTION :
  --   Evaluates the polynomial system at a random series
  --   with respect to the shift value in c.

  -- ON ENTRY :
  --   hom      square homotopy;
  --   sfh      shifted homotopy with respected to c;
  --   c        value of the shift;
  --   d        degree of the series used for evaluation.

  -- REQUIRED :
  --   The homotopy in hom has as many equations as unknowns.

    n : constant integer32 := hom'last;
    x : constant DoblDobl_Complex_Series_Vectors.Vector
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    cffhom : constant DoblDobl_Complex_Series_VecVecs.VecVec(hom'range)
           := DoblDobl_CSeries_Poly_SysFun.Coeff(hom);
    cshift : constant DoblDobl_Complex_Series_VecVecs.VecVec(hom'range)
           := DoblDobl_CSeries_Vector_Functions.Shift(cffhom,c);
    f : constant DoblDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(hom'range)
      := DoblDobl_CSeries_Poly_SysFun.Create(hom);
    y : constant DoblDobl_Complex_Series_Vectors.Vector
      := DoblDobl_CSeries_Poly_SysFun.Eval(hom,x);
    sy : constant DoblDobl_Complex_Series_Vectors.Vector
       := DoblDobl_CSeries_Poly_SysFun.Eval(sfh,x);
    z : constant DoblDobl_Complex_Series_Vectors.Vector
      := DoblDobl_CSeries_Poly_SysFun.Eval(f,cffhom,x);
    sz : constant DoblDobl_Complex_Series_Vectors.Vector
       := DoblDobl_CSeries_Poly_SysFun.Eval(f,cshift,x);
    ans : character;

  begin
    Symbol_Table.Init(1);
    Symbol_Table.Add_String("t");
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Evaluations after the shift ...");
      put_line("The value at the polynomial :"); put(sy);
      put_line("The value at the coefficient polynomial :"); put(sz);
    end if;
  end DoblDobl_Shift_Eval;

  procedure DoblDobl_Test_Shift ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in DoblDobl_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the system is shifted.

    h : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := DoblDobl_Homotopy.Homotopy_System;
    s : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
    zero : constant double_double := create(0.0);
    c : double_double := zero;
    shifteds : DoblDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);
    y,z : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    ans : character;
    deg : integer32 := 0;

  begin
    put("Give a value for the shift : "); get(c);
    shifteds := Series_and_Homotopies.Shift(s,c);
    y := Series_and_Homotopies.Eval(s,c);
    z := Series_and_Homotopies.Eval(shifteds,zero);
    put_line("system(shift value) :"); put_line(y);
    put_line("shifted system(0.0) :"); put_line(z);
    put("Continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the degree of the series in x : "); get(deg);
      DoblDobl_Shift_Eval(s,shifteds,c,deg);
    end if;
  end DoblDobl_Test_Shift;

  procedure QuadDobl_Shift_Eval
              ( hom : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                sfh : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
                c : in quad_double; d : in integer32 ) is

  -- DESCRIPTION :
  --   Evaluates the polynomial system at a random series
  --   with respect to the shift value in c.

  -- ON ENTRY :
  --   hom      square homotopy;
  --   sfh      shifted homotopy with respected to c;
  --   c        value of the shift;
  --   d        degree of the series used for evaluation.

  -- REQUIRED :
  --   The homotopy in hom has as many equations as unknowns.

    n : constant integer32 := hom'last;
    x : constant QuadDobl_Complex_Series_Vectors.Vector
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,n,d);
    cffhom : constant QuadDobl_Complex_Series_VecVecs.VecVec(hom'range)
           := QuadDobl_CSeries_Poly_SysFun.Coeff(hom);
    cshift : constant QuadDobl_Complex_Series_VecVecs.VecVec(hom'range)
           := QuadDobl_CSeries_Vector_Functions.Shift(cffhom,c);
    f : constant QuadDobl_CSeries_Poly_SysFun.Eval_Coeff_Poly_Sys(hom'range)
      := QuadDobl_CSeries_Poly_SysFun.Create(hom);
    y : constant QuadDobl_Complex_Series_Vectors.Vector
      := QuadDobl_CSeries_Poly_SysFun.Eval(hom,x);
    sy : constant QuadDobl_Complex_Series_Vectors.Vector
       := QuadDobl_CSeries_Poly_SysFun.Eval(sfh,x);
    z : constant QuadDobl_Complex_Series_Vectors.Vector
      := QuadDobl_CSeries_Poly_SysFun.Eval(f,cffhom,x);
    sz : constant QuadDobl_Complex_Series_Vectors.Vector
       := QuadDobl_CSeries_Poly_SysFun.Eval(f,cshift,x);
    ans : character;

  begin
    Symbol_Table.Init(1);
    Symbol_Table.Add_String("t");
    put_line("Evaluation at a random vector ..."); put(x);
    put_line("The value at the polynomial :"); put(y);
    put_line("The value at the coefficient polynomial :"); put(z);
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Evaluations after the shift ...");
      put_line("The value at the polynomial :"); put(sy);
      put_line("The value at the coefficient polynomial :"); put(sz);
    end if;
  end QuadDobl_Shift_Eval;

  procedure QuadDobl_Test_Shift ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in QuadDobl_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the system is shifted.

    h : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq)
      := QuadDobl_Homotopy.Homotopy_System;
    s : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
    zero : constant quad_double := create(0.0);
    c : quad_double := zero;
    shifteds : QuadDobl_CSeries_Poly_Systems.Poly_Sys(1..nq);
    y,z : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..nq);
    ans : character;
    deg : integer32 := 0;

  begin
    put("Give a value for the shift : "); get(c);
    shifteds := Series_and_Homotopies.Shift(s,c);
    y := Series_and_Homotopies.Eval(s,c);
    z := Series_and_Homotopies.Eval(shifteds,zero);
    put_line("system(shift value) :"); put_line(y);
    put_line("shifted system(0.0) :"); put_line(z);
    put("Continue ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("Give the degree of the series in x : "); get(deg);
      QuadDobl_Shift_Eval(s,shifteds,c,deg);
    end if;
  end QuadDobl_Test_Shift;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in standard double precision.

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := 2;
    gamma : constant Standard_Complex_Numbers.Complex_Number
          := Standard_Random_Numbers.Random1;
    ans : character;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system ..."); get(start);
    Standard_Homotopy.Create(target.all,start.all,k,gamma);
    new_line;
    put_line("MENU of testing operations : ");
    put_line("  0. test the creation of the series homotopy;");
    put_line("  1. test the evaluation of the series homotopy;");
    put_line("  2. test the shift of the series homotopy.");
    put("Type 0, 1, or 2 to select a test : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test_Creation(target'last);
      when '1' => Standard_Test_Evaluation(target'last);
      when '2' => Standard_Test_Shift(target'last);
      when others => null;
    end case;
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in double double precision.

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := 2;
    gamma : constant DoblDobl_Complex_Numbers.Complex_Number
          := DoblDobl_Random_Numbers.Random1;
    ans : character;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system ..."); get(start);
    DoblDobl_Homotopy.Create(target.all,start.all,k,gamma);
    new_line;
    put_line("MENU of testing operations : ");
    put_line("  0. test the creation of the series homotopy;");
    put_line("  1. test the evaluation of the series homotopy;");
    put_line("  2. test the shift of the series homotopy.");
    put("Type 0, 1, or 2 to select a test : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => DoblDobl_Test_Creation(target'last);
      when '1' => DoblDobl_Test_Evaluation(target'last);
      when '2' => DoblDobl_Test_Shift(target'last);
      when others => null;
    end case;
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on the operations of a homotopy with series coefficients,
  --   in quad double precision.

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    k : constant natural32 := 2;
    gamma : constant QuadDobl_Complex_Numbers.Complex_Number
          := QuadDobl_Random_Numbers.Random1;
    ans : character;

  begin
    new_line;
    put_line("Testing the creation of a homotopies as a series system ...");
    new_line;
    put_line("Reading the target system ..."); get(target);
    new_line;
    put_line("Reading the start system ..."); get(start);
    QuadDobl_Homotopy.Create(target.all,start.all,k,gamma);
    new_line;
    put_line("MENU of testing operations : ");
    put_line("  0. test the creation of the series homotopy;");
    put_line("  1. test the evaluation of the series homotopy;");
    put_line("  2. test the shift of the series homotopy.");
    put("Type 0, 1, or 2 to select a test : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => QuadDobl_Test_Creation(target'last);
      when '1' => QuadDobl_Test_Evaluation(target'last);
      when '2' => QuadDobl_Test_Shift(target'last);
      when others => null;
    end case;
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the working precision and launches the tests.

    ans : constant character := Prompt_for_Precision;

  begin
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_serhom;
