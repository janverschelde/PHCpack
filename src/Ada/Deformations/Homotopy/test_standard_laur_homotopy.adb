with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Laurent_Homotopy;

package body Test_Standard_Laur_Homotopy is

  procedure Interactive_Test ( p,q : in Laur_Sys ) is

    a : constant Complex_Number := Create(1.0);
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

    start,target : Link_to_Laur_Sys;

  begin
    new_line;
    put_line("Reading the target system..."); get(target);
    new_line;
    put_line("Reading the start system..."); get(start);
    Interactive_Test(target.all,start.all);
  end Main;

end Test_Standard_Laur_Homotopy;
