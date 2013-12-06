with text_io,integer_io;                use text_io,integer_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Norms;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Integer_Matrices;         use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Standard_Integer_Linear_Solvers;   use Standard_Integer_Linear_Solvers;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_Binomial_Systems;         use Standard_Binomial_Systems;
with Standard_Binomial_Solvers;         use Standard_Binomial_Solvers;
with Standard_Radial_Solvers;           use Standard_Radial_Solvers;

procedure ts_binwit is

-- DESCRIPTION :
--   Interactive development of binomial system solver.
--   The key operation is the calculation of the Hermite normal form.

  function Rank_Upper ( U : Matrix ) return natural is

  -- DESCRIPTION :
  --   Given an upper triangular matrix, returns its rank.

    res : natural := 0;
    pivot : integer := U'first(2);

  begin
    for i in U'range(1) loop
      while pivot <= U'last(2) and then U(i,pivot) = 0 loop
        pivot := pivot + 1;
      end loop;
      exit when pivot > U'last(2);
      res := i;
    end loop;
    return res;
  end Rank_Upper;

  procedure Normalize_Sign ( v : in out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   On return the orientation of the vector v might be changed 
  --   so the sign of its first nonzero element is positive.

    ispos : boolean := true; 

  begin
    for i in v'range loop
      if v(i) > 0 then
        return;
      elsif v(i) < 0 then
        exit;
      end if;
    end loop;
    Standard_Integer_Vectors.Min(v);
  end Normalize_Sign;

  procedure Kernel ( n,k : in natural; A : in Matrix;
                     r : out natural; V : out Matrix ) is

  -- DESCRIPTION :
  --   Computes the kernel of a k-by-n exponent matrix A.

  -- REQUIRED : k < n.

  -- ON ENTRY :
  --   n          the number of variables or the ambient dimension;
  --   k          the number of equations or the expected codimension;
  --   A          a matrix of k rows and n columns.

  -- ON RETURN :
  --   r          the rank of the matrix is expected to be k;
  --   V          a matrix of n rows and n-k columns, 
  --              the columns of V span the kernel of A.

    U : Matrix(A'range(1),A'range(2)) := A;
    B : Matrix(A'range(1),A'first(1)..A'last(1)+1);
    x : Standard_Integer_Vectors.Vector(B'range(2));

  begin
    Upper_Triangulate(U);
    put_line("After triangulation of A :"); put(U);
    r := Rank_Upper(U);
    put("The rank of the matrix : "); put(r,1); new_line;
    for i in B'range(1) loop
      for j in B'range(2) loop
        B(i,j) := U(i,j);
      end loop;
    end loop;
    for k in V'range(2) loop
      for i in B'range(1) loop
        B(i,B'last(2)) := U(i,U'last(1)+k);
      end loop;
      put_line("The matrix B :"); put(B);
      Solve0(B,x);
      put("The solution x : "); put(x); new_line;
      Normalize_Sign(x);
      for i in B'range(1) loop
        V(i,k) := x(i);
      end loop;
      for i in B'last(1)+1..V'last(1) loop
        V(i,k) := 0;
      end loop;
      V(B'last(1)+k,k) := x(x'last);
    end loop;
  end Kernel;

  procedure Compute_Kernel ( n,k : in natural ) is

  -- DESCRIPTION :
  --   Prompts for or generates a k-by-n matrix of exponents
  --   and computes its kernel.

  -- REQUIRED : k < n.

  -- ON ENTRY :
  --   n          the number of variables or the ambient dimension;
  --   k          the number of equations or the expected codimension.

    A : Matrix(1..k,1..n); 
    V : Matrix(1..n,1..n-k);
    r : natural;
    lower,upper : integer;
    ans : character;

  begin
    put("Generate a random matrix ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put("  give lower bound for exponents : "); get(lower);
      put("  give upper bound for exponents : "); get(upper);
      A := Random_Matrix(k,n,lower,upper);
      put("-> a random "); put(k,1); put("-by-"); put(n,1);
      put_line(" matrix :"); put(A);
    else
      put("-> give a "); put(k,1); put("-by-"); put(n,1);
      put_line(" matrix "); get(A);
    end if;
    Kernel(n,k,A,r,V);
    put_line("The kernel : "); put(V);
    put_line("The residuals : "); put(A*V);
  end Compute_Kernel;

  function Max_Delta ( x,y : Standard_Floating_Vectors.Vector )
                     return double_float is

  -- DESCRIPTION :
  --   Returns the maximal component of the difference x - y.

    res : double_float := abs(x(x'first) - y(y'first));
    dif : double_float;

  begin
    for i in x'first+1..x'last loop
      dif := abs(x(i) - y(i));
      if dif > res
       then res := dif;
      end if;
    end loop;
    return res;
  end Max_Delta;

  procedure Radial_Solve
              ( A : in Standard_Integer_Matrices.Matrix;
                c : in Standard_Complex_Vectors.Vector ) is

    r : constant Standard_Floating_Vectors.Vector(c'range) := Radii(c);
    logr : constant Standard_Floating_Vectors.Vector(c'range) := Log10(r);
    logx,x,y : Standard_Floating_Vectors.Vector(c'range);

  begin
    put_line("The radii of the complex numbers :"); put_line(r);
    put_line("The logarithms of the radii :"); put_line(logr);
    logx := Radial_Upper_Solve(A,logr);
    put_line("The logarithm of the solution : "); put_line(logx);
    y := Multiply(A,logx);
    put_line("The value A*x :"); put_line(y);
    put_line("the logarithms of the radii : "); put_line(logr);
    put("The residual : "); put(Max_Delta(y,logr),3); new_line;
    x := Exp10(logx);
    put_line("The solution x : "); put_line(x);
    y := Eval(A,x);
    put_line("The value x^A : "); put_line(y);
    put_line("The righthand side : "); put_line(r);
    put("The residual : "); put(Max_Delta(r,y),3); new_line;
  end Radial_Solve;

  procedure Compute_Tropism
               ( A : in Standard_Integer_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes a tropism to the exponent matrix of a binomial system
  --   in the standard format x^A - c = 0.
  --   The matrix has as many rows as there are unknowns in x
  --   and as many columns as there are equations, therefore
  --   we upper triangulate the transpose of the matrix A to
  --   compute a binomial tropism.

    B : Matrix(A'range(2),A'range(1)) := Transpose(A);
    U : Matrix(A'range(2),A'range(1)) := B;
    v : Standard_Integer_Vectors.Vector(A'range(1));
    r : Standard_Integer_Vectors.Vector(A'range(2));

  begin
    Upper_Triangulate(U);
    Solve0(U,v);
    Standard_Integer_Norms.Normalize(v);
    put("A tropism : "); put(v); new_line;
    r := B*v;
    put("the residual : "); put(r); new_line;
  end Compute_Tropism;

  procedure Test_Binomial_Solver
               ( A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

    U : Matrix(A'range(1),A'range(2)) := A;
    M : Matrix(A'range(1),A'range(1));
    r : natural;
    Usols,Asols : Solution_List;
    ans : character;

  begin
    Compute_Tropism(A);
    put("Do you want intermediate output of the solver ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Solve(standard_output,A,c,r,M,U,Usols,Asols);
     else Solve(A,c,r,M,U,Usols,Asols);
    end if;
    put("The rank of the matrix is "); put(r,1); put_line(".");
    if r < A'last(2) then
      put_line("-> system is rank deficient ...");
    else
      put("-> the solution set has dimension ");
      put(A'last(1)-A'last(2),1); put_line(".");
    end if;
    if not Is_Null(Usols) then
      put_line("The residuals at y^U - c, y = x^M : ");
      Write_Residuals(standard_output,U,c,Usols);
      put_line("The residuals at x^A - c : ");
      Write_Residuals(standard_output,A,c,Asols);
      put("Do you want to see the solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put(standard_output,Length_Of(Asols),Head_Of(Asols).n,Asols);
      end if;
    end if;
  end Test_Binomial_Solver;

  procedure Test_Solver ( p : in Poly_Sys ) is

    nq : constant natural := p'last;
    nv : constant natural := Number_of_Unknowns(p(p'first));
    fail : boolean;
    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    c : Standard_Complex_Vectors.Vector(1..nq);
 
  begin
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    Parse(p,nq,A,c,fail);
    if fail
     then put_line("The system is not a binomial system.");
     else put_line("The system is a binomial system.");
          put_line("The exponent matrix :"); put(A);
          put_line("The coefficient vector :"); put_line(c);
          Test_Binomial_Solver(A,c);
    end if;
  end Test_Solver;

  function Create_Random_Binomial_System
             ( nq,nv : natural ) return Poly_Sys is

    res : Poly_Sys(1..nq);
    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    c : constant Standard_Complex_Vectors.Vector(1..nq)
      := Random_Vector(1,nq);
    lower,upper : integer;

  begin
    put("  give lower bound for exponents : "); get(lower);
    put("  give upper bound for exponents : "); get(upper);
    A := Random_Matrix(nv,nq,lower,upper);
    put_line("The exponent matrix : "); put(A);
    put_line("The coefficient vector : "); put_line(c);
    res := Create(A,c);
    put_line("The binomial system : "); put(res);
    return res;
  end Create_Random_Binomial_System;

  procedure Read_and_Test_Binomial_Solver ( n : in natural ) is

    A : Standard_Integer_Matrices.Matrix(1..n,1..n);
    c : Standard_Complex_Vectors.Vector(1..n);

  begin
    put("Give an "); put(n,1); put("-by-"); put(n,1);
    put_line(" exponent matrix : "); get(A);
    put("Give an "); put(n,1);
    put_line("-vector for righthand side :"); get(c);
    put_line("The exponent matrix is "); put(A);
    put_line("The coefficient vector is "); put_line(c);
    Radial_Solve(A,c);
   -- Test_Binomial_Solver(A,c);
  end Read_and_Test_Binomial_Solver;
 
  procedure Main is

    lp : Link_to_Poly_Sys;
    ans : character;
    nq,nv : natural;

  begin
    new_line;
    put_line("Creating witness sets for binomial systems.");
    new_line;
    put_line("MENU for testing the operations on binomial systems :");
    put_line("  1. compute kernel of an exponent matrix;");
    put_line("  2. create and solve a random binomial system;");
    put_line("  3. give a polynomial system to test the solver;");
    put_line("  4. give a system x^A = b to test the solver.");
    put("Type 1, 2, 3, or 4 to select : "); 
    Ask_Alternative(ans,"1234");
    new_line;
    if ans = '1' then
      put_line("Computing the kernel of an exponent matrix...");
      put("  give the number of equations : "); get(nq);
      put("  give the number of variables : "); get(nv);
      Compute_Kernel(nv,nq);
    elsif ans = '4' then
      put("Give the dimension of the system : "); get(nq);
      Read_and_Test_Binomial_Solver(nq);
    else
      if ans = '2' then
        put_line("Creating a random binomial system...");
        put("  give the number of equations : "); get(nq);
        put("  give the number of variables : "); get(nv);
        lp := new Poly_Sys'(Create_Random_Binomial_System(nq,nv));
      else
        put_line("Reading a binomial system...");
        get(lp);
      end if;
      Test_Solver(lp.all);
    end if;
  end Main;

begin
  Main;
end ts_binwit;
