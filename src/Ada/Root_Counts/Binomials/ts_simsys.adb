with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;    use Standard_Integer64_Matrices_io;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with DoblDobl_Random_Vectors;           use DoblDobl_Random_Vectors;
with DoblDobl_Random_Matrices;          use DoblDobl_Random_Matrices;
with QuadDobl_Random_Vectors;           use QuadDobl_Random_Vectors;
with QuadDobl_Random_Matrices;          use QuadDobl_Random_Matrices;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laurentials;      use DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laurentials;      use QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Standard_Simplex_Systems;
with Standard_Simplex_Solvers;          use Standard_Simplex_Solvers;
with DoblDobl_Simplex_Systems;
with DoblDobl_Simplex_Solvers;          use DoblDobl_Simplex_Solvers;
with QuadDobl_Simplex_Systems;
with QuadDobl_Simplex_Solvers;          use QuadDobl_Simplex_Solvers;

procedure ts_simsys is

-- DESCRIPTION :
--   Interactive test of solving a user given or a random simplex system
--   in standard double, double double, and quad double arithmetic.

  procedure Standard_Parse_and_Solve
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 ) is

    use Standard_Complex_Solutions,Standard_Simplex_Systems;

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b : Standard_Complex_Vectors.Vector(1..nq);
    tol : constant double_float := 1.0E-8;
    sum,rcond : double_float;
    fail,zero_y : boolean;
    r : integer32;
    sols : Solution_List;
    file : file_type;
    ans : character;

  begin
    put_line("The system on input : "); put(p);
    Parse(p,nv,A,C,b,fail);
    if fail then
      put_line("The given system is not a simplex system...");
    else
      put_line("The given system is a simplex system,");
      put_line("its exponent matrix is "); put(A);
      put_line("Solving the simplex system ...");
      Solve(standard_output,A,C,b,tol,zero_y,r,rcond,sols);
      if zero_y
       then put_line("very small intermediate auxiliary solution components");
       else put_line("all intermediate solution components are nonzero");
      end if;
      put("estimate for inverse condition number : ");
      put(rcond,3); new_line;
      sum := Sum_Residuals(A,C,b,sols);
      put("The sum of the residuals : "); put(sum,3);
      put_line(".");
      put("Do you wish to see all residuals ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_Residuals(standard_output,A,C,b,sols);
      end if;
      new_line;
      put("Do you want to write the solutions to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(file);
        new_line;
        put_line("Writing the solutions to file...");
        new_line;
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end Standard_Parse_and_Solve;

  procedure DoblDobl_Parse_and_Solve
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 ) is

    use DoblDobl_Complex_Solutions,DoblDobl_Simplex_Systems;

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : DoblDobl_Complex_Vectors.Vector(1..nq);
    tol : constant double_double := create(1.0E-8);
    sum,rcond : double_double;
    fail,zero_y : boolean;
    r : integer32;
    sols : Solution_List;
    file : file_type;
    ans : character;

  begin
    put_line("The system on input : "); put(p);
    Parse(p,nv,A,C,b,fail);
    if fail then
      put_line("The given system is not a simplex system...");
    else
      put_line("The given system is a simplex system,");
      put_line("its exponent matrix is "); put(A);
      put_line("Solving the simplex system ...");
      Solve(standard_output,A,C,b,tol,zero_y,r,rcond,sols);
      if zero_y
       then put_line("very small intermediate auxiliary solution components");
       else put_line("all intermediate solution components are nonzero");
      end if;
      put("estimate for inverse condition number : ");
      put(rcond,3); new_line;
      sum := Sum_Residuals(A,C,b,sols);
      put("The sum of the residuals : "); put(sum,3);
      put_line(".");
      put("Do you wish to see all residuals ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_Residuals(standard_output,A,C,b,sols);
      end if;
      new_line;
      put("Do you want to write the solutions to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(file);
        new_line;
        put_line("Writing the solutions to file...");
        new_line;
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end DoblDobl_Parse_and_Solve;

  procedure QuadDobl_Parse_and_Solve
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv : in integer32 ) is

    use QuadDobl_Complex_Solutions,QuadDobl_Simplex_Systems;

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C : QuadDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b : QuadDobl_Complex_Vectors.Vector(1..nq);
    tol : constant quad_double := create(1.0E-8);
    sum,rcond : quad_double;
    fail,zero_y : boolean;
    r : integer32;
    sols : Solution_List;
    file : file_type;
    ans : character;

  begin
    put_line("The system on input : "); put(p);
    Parse(p,nv,A,C,b,fail);
    if fail then
      put_line("The given system is not a simplex system...");
    else
      put_line("The given system is a simplex system,");
      put_line("its exponent matrix is "); put(A);
      put_line("Solving the simplex system ...");
      Solve(standard_output,A,C,b,tol,zero_y,r,rcond,sols);
      if zero_y
       then put_line("very small intermediate auxiliary solution components");
       else put_line("all intermediate solution components are nonzero");
      end if;
      put("estimate for inverse condition number : ");
      put(rcond,3); new_line;
      sum := Sum_Residuals(A,C,b,sols);
      put("The sum of the residuals : "); put(sum,3);
      put_line(".");
      put("Do you wish to see all residuals ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Write_Residuals(standard_output,A,C,b,sols);
      end if;
      new_line;
      put("Do you want to write the solutions to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(file);
        new_line;
        put_line("Writing the solutions to file...");
        new_line;
        put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
    end if;
  end QuadDobl_Parse_and_Solve;

  procedure Standard_Parse_and_Solve_Simplex_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system 
  --   and calls the parse and solve procedure.

    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("#equations : "); put(nq,1); put("  ");
    put("#variables : "); put(nv,1); new_line;
    Standard_Parse_and_Solve(lp.all,nq,nv);
  end Standard_Parse_and_Solve_Simplex_System;

  procedure DoblDobl_Parse_and_Solve_Simplex_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system 
  --   and calls the parse and solve procedure.

    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("#equations : "); put(nq,1); put("  ");
    put("#variables : "); put(nv,1); new_line;
    DoblDobl_Parse_and_Solve(lp.all,nq,nv);
  end DoblDobl_Parse_and_Solve_Simplex_System;

  procedure QuadDobl_Parse_and_Solve_Simplex_System is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system 
  --   and calls the parse and solve procedure.

    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    nq,nv : integer32;

  begin
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    new_line;
    put("#equations : "); put(nq,1); put("  ");
    put("#variables : "); put(nv,1); new_line;
    QuadDobl_Parse_and_Solve(lp.all,nq,nv);
  end QuadDobl_Parse_and_Solve_Simplex_System;

  procedure Standard_Solve_Random_Simplex_System ( nq,nv : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for lower and upper bounds on the exponents,
  --   and generates random complex coefficients for a simplex system
  --   with nq equations and nv variables, which is then solved.

    res,parsed_res : Standard_Complex_Laur_Systems.Laur_Sys(1..nq);
    lower,upper,info : integer32 := 0;
    A,pA : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    C,pC : Standard_Complex_Matrices.Matrix(1..nq,1..nq);
    b,pb : Standard_Complex_Vectors.Vector(1..nq);
    tol : constant double_float := 1.0E-8;
    fail,zero_y : boolean;
    r : integer32;
    sols : Standard_Complex_Solutions.Solution_List;
    ans : character;

    use Standard_Complex_Laur_Systems,Standard_Simplex_Systems;
 
  begin
    put("  give lower bound for exponents : "); get(lower);
    put("  give upper bound for exponents : "); get(upper);
    A := Random_Matrix(natural32(nv),natural32(nq),lower,upper);
    put_line("The exponent matrix : "); put(A);
    C := Random_Matrix(natural32(nq),natural32(nq));
    b := Random_Vector(1,nq);
    res := Create(A,C,b);
    Parse(res,nv,pA,pC,pb,fail);
    if fail then
      put_line("Failed to parse as a simplex system! Bug!");
    else
      put_line("The parsed exponent matrix : "); put(pA);
      parsed_res := Create(pA,pC,pb);
      put_line("The parsed polynomial system :"); put_line(parsed_res);
      put_line("The created polynomial system :"); put_line(res);
      put_line("difference between the two systems : ");
      put_line(res-parsed_res);
      Standard_Complex_Laur_Systems.Clear(parsed_res);
    end if;
    new_line;
    put("Continue with solving ? "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Solve(standard_output,A,C,b,tol,zero_y,r,info,sols);
      Write_Residuals(standard_output,A,C,b,sols);
    end if;
  end Standard_Solve_Random_Simplex_System;

  procedure DoblDobl_Solve_Random_Simplex_System ( nq,nv : in integer32) is

    res,parsed_res : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..nq);
    info : integer32 := 0;
    lower,upper : integer64 := 0;
    A,pA : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C,pC : DoblDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b,pb : DoblDobl_Complex_Vectors.Vector(1..nq);
    tol : constant double_double := create(1.0E-8);
    fail,zero_y : boolean;
    r : integer32;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    ans : character;

    use DoblDobl_Complex_Laur_Systems,DoblDobl_Simplex_Systems;
 
  begin
    put("  give lower bound for exponents : "); get(lower);
    put("  give upper bound for exponents : "); get(upper);
    A := Random_Matrix(natural32(nv),natural32(nq),lower,upper);
    put_line("The exponent matrix : "); put(A);
    C := Random_Matrix(natural32(nq),natural32(nq));
    b := Random_Vector(1,nq);
    res := Create(A,C,b);
    Parse(res,nv,pA,pC,pb,fail);
    if fail then
      put_line("Failed to parse as a simplex system! Bug!");
    else
      put_line("The parsed exponent matrix : "); put(pA);
      parsed_res := Create(pA,pC,pb);
      put_line("The parsed polynomial system :"); put_line(parsed_res);
      put_line("The created polynomial system :"); put_line(res);
      put_line("difference between the two systems : ");
      put_line(res-parsed_res);
      DoblDobl_Complex_Laur_Systems.Clear(parsed_res);
    end if;
    new_line;
    put("Continue with solving ? "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Solve(standard_output,A,C,b,tol,zero_y,r,info,sols);
      Write_Residuals(standard_output,A,C,b,sols);
    end if;
  end DoblDobl_Solve_Random_Simplex_System;

  procedure QuadDobl_Solve_Random_Simplex_System ( nq,nv : in integer32) is

    res,parsed_res : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..nq);
    info : integer32 := 0;
    lower,upper : integer64 := 0;
    A,pA : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    C,pC : QuadDobl_Complex_Matrices.Matrix(1..nq,1..nq);
    b,pb : QuadDobl_Complex_Vectors.Vector(1..nq);
    tol : constant quad_double := create(1.0E-8);
    fail,zero_y : boolean;
    r : integer32;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    ans : character;

    use QuadDobl_Complex_Laur_Systems,QuadDobl_Simplex_Systems;
 
  begin
    put("  give lower bound for exponents : "); get(lower);
    put("  give upper bound for exponents : "); get(upper);
    A := Random_Matrix(natural32(nv),natural32(nq),lower,upper);
    put_line("The exponent matrix : "); put(A);
    C := Random_Matrix(natural32(nq),natural32(nq));
    b := Random_Vector(1,nq);
    res := Create(A,C,b);
    Parse(res,nv,pA,pC,pb,fail);
    if fail then
      put_line("Failed to parse as a simplex system! Bug!");
    else
      put_line("The parsed exponent matrix : "); put(pA);
      parsed_res := Create(pA,pC,pb);
      put_line("The parsed polynomial system :"); put_line(parsed_res);
      put_line("The created polynomial system :"); put_line(res);
      put_line("difference between the two systems : ");
      put_line(res-parsed_res);
      QuadDobl_Complex_Laur_Systems.Clear(parsed_res);
    end if;
    new_line;
    put("Continue with solving ? "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      Solve(standard_output,A,C,b,tol,zero_y,r,info,sols);
      Write_Residuals(standard_output,A,C,b,sols);
    end if;
  end QuadDobl_Solve_Random_Simplex_System;

  procedure Main is

    ans : character;
    nq,nv : integer32 := 0;

  begin
    new_line;
    put_line("Testing the operations in simplex systems ...");
    new_line;
    put_line("MENU to test the simplex system solvers :");
    put_line("  1. solve a random simplex system in standard arithmetic;");
    put_line("  2. solve a random simplex system with double doubles;");
    put_line("  3. solve a random simplex system with quad doubles;");
    put_line("  4. parse and solve a simplex system in standard arithmetic;");
    put_line("  5. parse and solve a simplex system with double doubles;");
    put_line("  6. parse and solve a simplex system with quad doubles.");
    put("Type 1, 2, 3, 4, 5, or 6 to make your choice : ");
    Ask_Alternative(ans,"123456");
    new_line;
    if ans = '1' or ans = '2' or ans = '3' then
      put("Give the number of equations : "); get(nq);
      nv := nq; -- put("Give the number of variables : "); get(nv);
    end if;
    case ans is
      when '1' => Standard_Solve_Random_Simplex_System(nq,nv);
      when '2' => DoblDobl_Solve_Random_Simplex_System(nq,nv);
      when '3' => QuadDobl_Solve_Random_Simplex_System(nq,nv);
      when '4' => Standard_Parse_and_Solve_Simplex_System;
      when '5' => DoblDobl_Parse_and_Solve_Simplex_System;
      when '6' => QuadDobl_Parse_and_Solve_Simplex_System;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_simsys;
