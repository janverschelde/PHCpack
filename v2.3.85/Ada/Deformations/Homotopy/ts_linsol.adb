with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Linear_Poly_Solvers;       use Standard_Linear_Poly_Solvers;

procedure ts_linsol is

  procedure Solve_Linear_System ( p : in Poly_Sys ) is

    neq : constant integer32 := p'last;
    nvr : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    A : Matrix(1..neq,1..nvr);
    b : Vector(1..neq);
    s : Solution(nvr);

  begin
    put("Solving a linear system of ");
    put(neq,1); put(" equations in ");
    put(nvr,1); put_line(" variables...");
    Coefficients(p,A,b);
    put_line("The coefficients of the system :"); put(A);
    put_line("The right-hand side vector : "); put_line(b);
    if neq /= nvr
     then put_line("nonsquare solving not yet implemented...");
     else s := Square_Solve(A,b);
          put_line("The solution :"); put(s);
    end if;
  end Solve_Linear_System;

  procedure Main is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Testing the solution of a linear system ...");
    new_line;
    get(lp);
    new_line;
    if not Is_Linear(lp.all)
     then put_line("The system is NONlinear.");
     else put_line("The system is linear.");
          Solve_Linear_System(lp.all);
    end if;
  end Main;
 
begin
  Main;
end ts_linsol;
