with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Homotopy;                 use Standard_Homotopy;
with Standard_Laurent_Homotopy;         use Standard_Laurent_Homotopy;

procedure ts_homotopy is

-- DESCRIPTION :
--   Test on homotopy creation and clearing.

  procedure Test_Homotopy ( p,q : in Poly_Sys ) is

    a : constant Complex_Number := Random1;
    n : natural32 := 0;

  begin
    new_line;
    put_line("creating and clearing homotopies...");
    new_line;
    put("Give the number of tests : "); get(n); new_line;
    for i in 1..n loop
      put("before create #"); put(i,1); put_line(" ...");
      Standard_Homotopy.Create(p,q,2,a);
      put("after create and before clear #"); put(i,1); put_line(" ...");
      Standard_Homotopy.Clear;
      put("after clear #"); put(i,1); put_line(" ...");
    end loop;
  end Test_Homotopy;

  procedure Random_Test_Homotopy_Eval ( p,q : in Poly_Sys ) is

    a : constant Complex_Number := Random1;
    n : natural32 := 0;
    nbv : constant integer32
        := integer32(Number_of_Unknowns(p(p'first)));

  begin
    new_line;
    put_line("testing evaluation in homotopy...");
    new_line;
    put("Give the number of samples : "); get(n); new_line;
    Standard_Homotopy.Create(p,q,2,a);
    for i in 1..n loop
      declare
        t : constant Complex_Number := Random1;
	x : Vector(1..nbv) := Random_Vector(1,nbv);
        y : Vector(p'range) := Standard_Homotopy.Eval(x,t);
      begin
        put_line("The evaluation at a random point :");
        put_line(y);
      end;
    end loop;
    Standard_Homotopy.Clear;
  end Random_Test_Homotopy_Eval;

  procedure Interactive_Test ( p,q : in Poly_Sys ) is

    a : Complex_Number := Create(1.0);
    ans : character;
    nbv : constant integer32
        := integer32(Number_of_Unknowns(p(p'first)));

  begin
    new_line;
    put_line("testing evaluation in homotopy...");
    Standard_Homotopy.Create(p,q,2,a);
    put("Dimension of the homotopy : ");
    put(Standard_Homotopy.Dimension,1); new_line;
    loop
      declare
        t : Complex_Number;
	x : Vector(1..nbv);
	y : Vector(p'range);
      begin
        new_line;
        put("Give a complex number for t : ");
	get(t);
        put("Give "); put(nbv,1); put_line(" complex numbers for x :");
        get(x);
        put("t = "); put(t); new_line;
        put_line("The point x : "); put_line(x);
        y := Standard_Homotopy.Eval(x,t);
        put_line("The evaluation at the point (x,t) :");
        put_line(y);
      end;
      new_line;
      put("Do you want to evaluate at another point ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Homotopy.Clear;
  end Interactive_Test;

  procedure Interactive_Test ( p,q : in Laur_Sys ) is

    a : Complex_Number := Create(1.0);
    ans : character;
    nbv : constant integer32
        := integer32(Number_of_Unknowns(p(p'first)));

  begin
    new_line;
    put_line("testing evaluation in homotopy...");
    new_line;
    put("Number of polynomials : "); put(p'last,1); new_line;
    put("Number of variables : "); put(nbv,1); new_line;
    Standard_Laurent_Homotopy.Create(p,q,2,a);
    loop
      declare
        t : Complex_Number;
	x : Vector(1..nbv);
        y : Vector(p'range);
      begin
        new_line;
        put("Give a complex number for t : ");
	get(t);
        put("Give "); put(nbv,1); put_line(" complex numbers for x :");
        get(x);
        put("t = "); put(t); new_line;
        put_line("The point x : "); put_line(x);
        y := Standard_Laurent_Homotopy.Eval(x,t);
        put_line("The evaluation at the point (x,t) :");
        put_line(y);
      end;
      new_line;
      put("Do you want to evaluate at another point ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Standard_Laurent_Homotopy.Clear;
  end Interactive_Test;

  procedure Main is

    start,target : Link_to_Poly_Sys;
    lstart,ltarget : Link_to_Laur_Sys;
    ans : character;
    laurent : boolean;

  begin
    new_line;
    put_line("Testing homotopy create and clear operations...");
    new_line;
    put("Laurent polynomial systems ? (y/n) ");
    Ask_Yes_or_No(ans);
    laurent := (ans = 'y');
    new_line;
    put_line("Reading the target system...");
    if laurent
     then get(ltarget);
     else get(target);
    end if;
    new_line;
    put_line("Reading the start system...");
    if laurent
     then get(lstart);
     else get(start);
    end if;
    if laurent then
      ans := '4';
    else
      new_line;
      put_line("MENU to test the homotopy package :");
      put_line("  1. test creation and clear of homotopy;");
      put_line("  2. test evaluation at random points;");
      put_line("  3. test evaluation at user given points.");
      put("Type 1, 2, or 3 to select : ");
      Ask_Alternative(ans,"123");
    end if;
    case ans is
      when '1' => Test_Homotopy(target.all,start.all);
      when '2' => Random_Test_Homotopy_Eval(target.all,start.all);
      when '3' => Interactive_Test(target.all,start.all);
      when '4' => Interactive_Test(ltarget.all,lstart.all);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_homotopy;
