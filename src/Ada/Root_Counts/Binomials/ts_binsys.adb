with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with QuadDobl_Random_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Multprec_Integer_Matrices;
with Standard_Random_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with Standard_Binomial_Systems;          use Standard_Binomial_Systems;
with Standard_Binomial_Solvers;          use Standard_Binomial_Solvers;
with DoblDobl_Binomial_Systems;          use DoblDobl_Binomial_Systems;
with DoblDobl_Binomial_Solvers;          use DoblDobl_Binomial_Solvers;
with QuadDobl_Binomial_Systems;          use QuadDobl_Binomial_Systems;
with QuadDobl_Binomial_Solvers;          use QuadDobl_Binomial_Solvers;

procedure ts_binsys is

-- DESCRIPTION :
--   Testing the binomial system solvers.

  function Convert ( A : Standard_Integer_Matrices.Matrix )
                   return Multprec_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the matrix A converted to multiprecision type.

    res : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Create(A(i,j));
      end loop;
    end loop;
    return res;
  end Convert;

  function Convert ( A : Standard_Integer64_Matrices.Matrix )
                   return Multprec_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the matrix A converted to multiprecision type.

    res : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Create(integer32(A(i,j)));
      end loop;
    end loop;
    return res;
  end Convert;

  procedure Standard_Solver
              ( A : in Standard_Integer_Matrices.Matrix;
                c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the solver for x^A = c and computes residuals.

    use Standard_Complex_Solutions;

    M,U : Standard_Integer_Matrices.Matrix(A'range(1),A'range(2));
    d : natural32;
    r : integer32;
    deg : integer32;
    sols : Solution_List;

  begin
    Solve(A,c,r,M,U,sols);
    deg := Standard_Binomial_Solvers.Degree(U);
    put("The determinant of the upper triangular matrix : ");
    put(deg,1); new_line;
    put("The rank of the exponent matrix is "); put(r,1); new_line;
    if r /= 0 then
      d := Length_Of(sols);
      put("Found "); put(d,1); put_line(" solutions.");
      if d > 0 then
        put("Sum of residuals : ");
        put(Sum_Residuals(A,c,sols),3); new_line;
      end if;
    end if;
  end Standard_Solver;

  procedure DoblDobl_Solver
              ( A : in Standard_Integer64_Matrices.Matrix;
                c : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the solver for x^A = c and computes residuals.

    use DoblDobl_Complex_Solutions;

    M,U : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2));
    d : natural32;
    r : integer32;
    deg : integer64;
    sols : Solution_List;

  begin
    Solve(A,c,r,M,U,sols);
    deg := DoblDobl_Binomial_Solvers.Degree(U);
    put("The determinant of the upper triangular matrix : ");
    put(deg,1); new_line;
    put("The rank of the exponent matrix is "); put(r,1); new_line;
    if r /= 0 then
      d := Length_Of(sols);
      put("Found "); put(d,1); put_line(" solutions.");
      if d > 0 then
        put("Sum of residuals : ");
        put(Sum_Residuals(A,c,sols),3); new_line;
      end if;
    end if;
  end DoblDobl_Solver;

  procedure QuadDobl_Solver
              ( A : in Standard_Integer64_Matrices.Matrix;
                c : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the solver for x^A = c and computes residuals.

    use QuadDobl_Complex_Solutions;

    M,U : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2));
    d : natural32;
    r : integer32;
    deg : integer64;
    sols : Solution_List;

  begin
    Solve(A,c,r,M,U,sols);
    deg := QuadDobl_Binomial_Solvers.Degree(U);
    put("The determinant of the upper triangular matrix : ");
    put(deg,1); new_line;
    put("The rank of the exponent matrix is "); put(r,1); new_line;
    if r /= 0 then
      d := Length_Of(sols);
      put("Found "); put(d,1); put_line(" solutions.");
      if d > 0 then
        put("Sum of residuals : ");
        put(Sum_Residuals(A,c,sols),3); new_line;
      end if;
    end if;
  end QuadDobl_Solver;

  procedure Standard_Solver_with_Multiprecision_Hermite
              ( A : in Standard_Integer_Matrices.Matrix;
                c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the Hermite normal form with multiprecision arithmetic.

    use Standard_Complex_Solutions;

    mpA,M,U : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
    d : natural32;
    r : integer32;
    sols : Solution_List;

  begin
    mpA := Convert(A);
    Solve(mpA,c,r,M,U,sols);
    put("The rank of the exponent matrix is "); put(r,1); new_line;
    if r /= 0 then
      d := Length_Of(sols);
      put("Found "); put(d,1); put_line(" solutions.");
      put(standard_Output,d,natural32(c'last),sols);
      if d > 0 then
        put("Sum of residuals : ");
        put(Sum_Residuals(A,c,sols),3); new_line;
      end if;
    end if;
  end Standard_Solver_with_Multiprecision_Hermite;

  procedure DoblDobl_Solver_with_Multiprecision_Hermite
              ( A : in Standard_Integer64_Matrices.Matrix;
                c : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the Hermite normal form with multiprecision arithmetic.

    use DoblDobl_Complex_Solutions;

    mpA,M,U : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
    d : natural32;
    r : integer32;
    sols : Solution_List;

  begin
    mpA := Convert(A);
    Solve(mpA,c,r,M,U,sols);
    put("The rank of the exponent matrix is "); put(r,1); new_line;
    if r /= 0 then
      d := Length_Of(sols);
      put("Found "); put(d,1); put_line(" solutions.");
      put(standard_Output,d,natural32(c'last),sols);
      if d > 0 then
        put("Sum of residuals : ");
        put(Sum_Residuals(A,c,sols),3); new_line;
      end if;
    end if;
  end DoblDobl_Solver_with_Multiprecision_Hermite;

  procedure Standard_Parse_and_Solve
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                nq,nv : in integer32; mp : in boolean ) is

  -- DESCRIPTION :
  --   Parses the system p with nq equations in nv variables
  --   into a binomial system and if parsing succeeds,
  --   then calls the solver.
  --   The flag mp indicates whether to compute the Hermite normal
  --   form with multiprecision arithmetic or not.

    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    c : Standard_Complex_Vectors.Vector(1..nq);
    fail : boolean;

  begin
    Parse(p,nq,A,c,fail);
    if fail then
      put_line("The system is not a binomial system.");
    else
      put_line("The system is a binomial system.");
      if mp 
       then Standard_Solver_with_Multiprecision_Hermite(A,c);
       else Standard_Solver(A,c);
      end if;
    end if;
  end Standard_Parse_and_Solve;

  procedure DoblDobl_Parse_and_Solve
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                nq,nv : in integer32 ) is

  -- DESCRIPTION :
  --   Parses the system into a binomial system and if parsing succeeds,
  --   then calls the solver.

    use DoblDobl_Complex_Solutions;

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    c : DoblDobl_Complex_Vectors.Vector(1..nq);
    fail : boolean;
    sols : Solution_List;

  begin
    Parse(p,nq,A,c,fail);
    if fail then
      put_line("The system is not a binomial system.");
    else
      put_line("The system is a binomial system.");
      DoblDobl_Solver(A,c);
    end if;
  end DoblDobl_Parse_and_Solve;

  procedure QuadDobl_Parse_and_Solve
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                nq,nv : in integer32 ) is

  -- DESCRIPTION :
  --   Parses the system into a binomial system and if parsing succeeds,
  --   then calls the solver.

    use QuadDobl_Complex_Solutions;

    A : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    c : QuadDobl_Complex_Vectors.Vector(1..nq);
    fail : boolean;
    sols : Solution_List;

  begin
    Parse(p,nq,A,c,fail);
    if fail then
      put_line("The system is not a binomial system.");
    else
      put_line("The system is a binomial system.");
      QuadDobl_Solver(A,c);
    end if;
  end QuadDobl_Parse_and_Solve;

  procedure Standard_Test_Solver ( multprec_hermite : in boolean ) is

  -- DESCRIPTION :
  --   Prompts the user for a system and calls parse & solve.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    nq,nv : integer32;

  begin
    new_line;
    get(lp);
    new_line;
    put_line("The polynomial system :"); put(lp.all);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'last)));
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    Standard_Parse_and_Solve(lp.all,nq,nv,multprec_hermite);
  end Standard_Test_Solver;

  procedure DoblDobl_Test_Solver is

  -- DESCRIPTION :
  --   Prompts the user for a system and calls parse & solve.

    lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nq,nv : integer32;

  begin
    new_line;
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'last)));
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    DoblDobl_Parse_and_Solve(lp.all,nq,nv);
  end DoblDobl_Test_Solver;

  procedure QuadDobl_Test_Solver is

  -- DESCRIPTION :
  --   Prompts the user for a system and calls parse & solve.

    lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    nq,nv : integer32;

  begin
    new_line;
    get(lp);
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'last)));
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    QuadDobl_Parse_and_Solve(lp.all,nq,nv);
  end QuadDobl_Test_Solver;

  procedure Standard_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for dimension and degrees, generates a random
  --   integer exponent matrix A, corresponding coefficient vector c,
  --   and then solves the system x^A = c.

    dim : natural32 := 0;
    low,upp : integer32 := 0;

  begin
    new_line;
    put("give the dimension : "); get(dim);
    put("give a lower bound on the exponents : "); get(low);
    put("give an upper bound on the exponents : "); get(upp);
    declare
      n : constant integer32 := integer32(dim);
      A : Standard_Integer_Matrices.Matrix(1..n,1..n)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,low,upp);
      c : Standard_Complex_Vectors.Vector(1..n)
        := Standard_Random_Vectors.Random_Vector(1,n);
    begin
      put_line("The exponents : "); put(A);
      put_line("The coefficients : "); put_line(c);
      Standard_Solver(A,c);
    end;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for dimension and degrees, generates a random
  --   integer exponent matrix A, corresponding coefficient vector c,
  --   and then solves the system x^A = c.

    dim : natural32 := 0;
    low,upp : integer64 := 0;

  begin
    new_line;
    put("give the dimension : "); get(dim);
    put("give a lower bound on the exponents : "); get(low);
    put("give an upper bound on the exponents : "); get(upp);
    declare
      n : constant integer32 := integer32(dim);
      A : Standard_Integer64_Matrices.Matrix(1..n,1..n)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,low,upp);
      c : DoblDobl_Complex_Vectors.Vector(1..n)
        := DoblDobl_Random_Vectors.Random_Vector(1,n);
    begin
      put_line("The exponents : "); put(A);
      put_line("The coefficients : "); put_line(c);
      DoblDobl_Solver(A,c);
    end;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test is

  -- DESCRIPTION :
  --   Prompts the user for dimension and degrees, generates a random
  --   integer exponent matrix A, corresponding coefficient vector c,
  --   and then solves the system x^A = c.

    dim : natural32 := 0;
    low,upp : integer64 := 0;

  begin
    new_line;
    put("give the dimension : "); get(dim);
    put("give a lower bound on the exponents : "); get(low);
    put("give an upper bound on the exponents : "); get(upp);
    declare
      n : constant integer32 := integer32(dim);
      A : Standard_Integer64_Matrices.Matrix(1..n,1..n)
        := Standard_Random_Matrices.Random_Matrix(dim,dim,low,upp);
      c : QuadDobl_Complex_Vectors.Vector(1..n)
        := QuadDobl_Random_Vectors.Random_Vector(1,n);
    begin
      put_line("The exponents : "); put(A);
      put_line("The coefficients : "); put_line(c);
      QuadDobl_Solver(A,c);
    end;
  end QuadDobl_Random_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test binomial system solvers :");
    put_line("  1. parse and solve in standard arithmetic;");
    put_line("  2. parse and solve in double double arithmetic;");
    put_line("  3. parse and solve in quad double arithmetic;");
    put_line("  4. parse and solve with multiprecision Hermite normal form;");
    put_line("  5. solve a random system in standard arithmetic.");
    put_line("  6. solve a random system in double double arithmetic.");
    put_line("  7. solve a random system in quad double arithmetic.");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to choose : ");
    Ask_Alternative(ans,"1234567");
    case ans is
      when '1' => Standard_Test_Solver(false);
      when '2' => DoblDobl_Test_Solver;
      when '3' => QuadDobl_Test_Solver;
      when '4' => Standard_Test_Solver(true);
      when '5' => Standard_Random_Test;
      when '6' => DoblDobl_Random_Test;
      when '7' => QuadDobl_Random_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_binsys;
