with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Homotopy;
with Standard_Series_Poly_Systems;
with Series_and_Polynomials;
with Series_and_Polynomials_io;          use Series_and_Polynomials_io;
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

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_Series_Poly_Systems.Poly_Sys(1..nq);

  begin
    put_line("The homotopy system :"); put_line(h);
    s := Series_and_Homotopies.Create(h,nq+1,true);
    put_line("The series system :"); put(s,1);
  end Standard_Test_Creation;

  procedure Standard_Test_Evaluation ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in Standard_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the evaluated
  --   polynomial system is shown.

    p : Standard_Complex_Poly_Systems.Poly_Sys(1..nq);
    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_Series_Poly_Systems.Poly_Sys(1..nq)
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

  procedure Standard_Test_Shift ( nq : in integer32 ) is

  -- DESCRIPTION :
  --   Given in Standard_Homotopy is a homotopy system of nq equations.
  --   The user is prompted for a value of t and the system is shifted.

    h : Standard_Complex_Poly_Systems.Poly_Sys(1..nq)
      := Standard_Homotopy.Homotopy_System;
    s : Standard_Series_Poly_Systems.Poly_Sys(1..nq)
      := Series_and_Homotopies.Create(h,nq+1);
    c : double_float := 0.0;
    shifteds : Standard_Series_Poly_Systems.Poly_Sys(1..nq);
    y,z : Standard_Complex_Poly_Systems.Poly_Sys(1..nq);

  begin
    put("Give a value for the shift : "); get(c);
    shifteds := Series_and_Homotopies.Shift(s,c);
    y := Series_and_Homotopies.Eval(s,c);
    z := Series_and_Homotopies.Eval(shifteds,0.0);
    put_line("system(shift value) :"); put_line(y);
    put_line("shifted system(0.0) :"); put_line(z);
  end Standard_Test_Shift;

  procedure Main is

  -- DESCRIPTION :
  --   Test on converting a homotopy system to a series system.

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
  end Main;

begin
  Main;
end ts_serhom;
