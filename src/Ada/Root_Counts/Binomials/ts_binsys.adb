with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with Standard_Integer_Matrices;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with Standard_Binomial_Systems;          use Standard_Binomial_Systems;
with Standard_Binomial_Solvers;          use Standard_Binomial_Solvers;
with DoblDobl_Binomial_Systems;          use DoblDobl_Binomial_Systems;
with DoblDobl_Binomial_Solvers;          use DoblDobl_Binomial_Solvers;

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

  procedure Standard_Solver
              ( A : in Standard_Integer_Matrices.Matrix;
                c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Calls the solver for x^A = c and computes residuals.

    use Standard_Complex_Solutions;

    M,U : Standard_Integer_Matrices.Matrix(A'range(1),A'range(2));
    d : natural32;
    r : integer32;
    sols : Solution_List;

  begin
    Solve(A,c,r,M,U,sols);
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

  procedure Solver_with_Multiprecision_Hermite
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
  end Solver_with_Multiprecision_Hermite;

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
       then Solver_with_Multiprecision_Hermite(A,c);
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

    A,M,U : Standard_Integer64_Matrices.Matrix(1..nv,1..nq);
    c : DoblDobl_Complex_Vectors.Vector(1..nq);
    r : integer32;
    d : natural32;
    fail : boolean;
    sols : Solution_List;

  begin
    Parse(p,nq,A,c,fail);
    if fail then
      put_line("The system is not a binomial system.");
    else
      put_line("The system is a binomial system.");
      Solve(A,c,r,M,U,sols);
      put("The rank of the exponent matrix is "); put(r,1); new_line;
      if r /= 0 then
        d := Length_Of(sols);
        put("Found "); put(d,1); put_line(" solutions.");
        if d > 0 then
          put("Sum of residuals : ");
          put(Sum_Residuals(A,c,sols),3); new_line;
        end if;
      end if;
    end if;
  end DoblDobl_Parse_and_Solve;

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

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test binomial system solvers :");
    put_line("  1. parse and solve in standard arithmetic;");
    put_line("  2. parse and solve with double doubles;");
    put_line("  3. parse and solve with multiprecision Hermite normal form.");
    put("Type 1, 2, or 3 to choose : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Standard_Test_Solver(false);
      when '2' => DoblDobl_Test_Solver;
      when '3' => Standard_Test_Solver(true);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_binsys;
