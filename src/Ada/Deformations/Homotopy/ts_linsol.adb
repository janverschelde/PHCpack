with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Standard_Linear_Poly_Solvers;
with DoblDobl_Linear_Poly_Solvers;
with QuadDobl_Linear_Poly_Solvers;

procedure ts_linsol is

-- DESCRIPTION :
--   Interactive test on the solving of linear systems given
--   in symbolic form as polynomial systems.

  procedure Solve_Linear_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Given a linear system in symbolic polynomial format,
  --   extracts the coefficient matrix and right hand side vector,
  --   and then solves the linear system with standard double arithmetic.

    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Solutions;
    use Standard_Linear_Poly_Solvers;

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
    if neq /= nvr then
      put_line("nonsquare solving not yet implemented...");
    else
      s := Square_Solve(A,b);
      put_line("The solution :"); put(s);
    end if;
  end Solve_Linear_System;

  procedure Solve_Linear_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Given a linear system in symbolic polynomial format,
  --   extracts the coefficient matrix and right hand side vector,
  --   and then solves the linear system with double double arithmetic.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Solutions;
    use DoblDobl_Linear_Poly_Solvers;

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
    if neq /= nvr then
      put_line("nonsquare solving not yet implemented...");
    else
      s := Square_Solve(A,b);
      put_line("The solution :"); put(s);
    end if;
  end Solve_Linear_System;

  procedure Solve_Linear_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Given a linear system in symbolic polynomial format,
  --   extracts the coefficient matrix and right hand side vector,
  --   and then solves the linear system with quad double arithmetic.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Solutions;
    use QuadDobl_Linear_Poly_Solvers;

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
    if neq /= nvr then
      put_line("nonsquare solving not yet implemented...");
    else
      s := Square_Solve(A,b);
      put_line("The solution :"); put(s);
    end if;
  end Solve_Linear_System;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with coefficients and
  --   standard double precision.  If the system parses into a linear
  --   systems, then it is solved with standard double arithmetic.

    use Standard_Complex_Poly_Systems;
    use Standard_Linear_Poly_Solvers;

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    get(lp);
    new_line;
    if not Is_Linear(lp.all) then
      put_line("The system is NONlinear.");
    else
      put_line("The system is linear.");
      Solve_Linear_System(lp.all);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with coefficients and
  --   standard double precision.  If the system parses into a linear
  --   systems, then it is solved with double double arithmetic.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Linear_Poly_Solvers;

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    get(lp);
    new_line;
    if not Is_Linear(lp.all) then
      put_line("The system is NONlinear.");
    else
      put_line("The system is linear.");
      Solve_Linear_System(lp.all);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system with coefficients and
  --   standard double precision.  If the system parses into a linear
  --   systems, then it is solved with quad double arithmetic.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Linear_Poly_Solvers;

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    get(lp);
    new_line;
    if not Is_Linear(lp.all) then
      put_line("The system is NONlinear.");
    else
      put_line("The system is linear.");
      Solve_Linear_System(lp.all);
    end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then calls the corresponding tester.

    ans : character;

  begin
    new_line;
    put_line("Testing the solution of a linear system ...");
    put_line("  0. test in standard double precision;");
    put_line("  1. test in double double precision;");
    put_line("  2. test in quad double precision.");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;
 
begin
  Main;
end ts_linsol;
