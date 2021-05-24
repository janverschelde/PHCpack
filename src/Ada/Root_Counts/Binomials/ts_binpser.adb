with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Common_Divisors;          use Standard_Common_Divisors;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Norms;
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
with Standard_Binomial_Systems;         use Standard_Binomial_Systems;
with Standard_Power_Transformations;    use Standard_Power_Transformations;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Poly_Laur_Convertors;     use Standard_Poly_Laur_Convertors;
with Standard_Laur_Poly_Convertors;     use Standard_Laur_Poly_Convertors;
with Standard_Integer32_Transformations;
 use Standard_Integer32_Transformations;
with Transforming_Laurent_Systems;      use Transforming_Laurent_Systems;

procedure ts_binpser is

-- DESCRIPTION :
--   Interactive development of Puiseux series for binomial systems
--   in the standard form x^A = c.

  function Tropism ( A : in Standard_Integer_Matrices.Matrix )
                   return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a tropism to the exponent matrix of a binomial system
  --   in the standard format x^A - c = 0.
  --   The matrix has as many rows as there are unknowns in x
  --   and as many columns as there are equations, therefore
  --   v'range equals A'range(1).

    res : Standard_Integer_Vectors.Vector(A'range(1));
    B : constant Matrix(A'range(2),A'range(1)) := Transpose(A);
    U : Matrix(A'range(2),A'range(1)) := B;

  begin
    Upper_Triangulate(U);
    Solve0(U,res);
    Standard_Integer_Norms.Normalize(res);
    if res(res'first) < 0
     then Standard_Integer_Vectors.Min(res);
    end if; 
    return res;
  end Tropism;

  function Degree ( v : Standard_Integer_Vectors.Vector ) return integer32 is

  -- DESCRIPTION :
  --   The degree of a space curve defined by a binomial system
  --   in the standard form x^A = c with tropism v: A'*v = 0 is
  --   defined as the maximal entry in t minus the minimal entry in v.

    max : integer32 := v(v'first);
    min : integer32 := v(v'first);
    d : integer32;

  begin
    for k in v'first+1..v'last loop
      if v(k) > max then
        max := v(k);
      elsif v(k) < min then
        min := v(k);
      end if;
    end loop;
    d := max - min;
    if d < 0
     then return -d;
     else return d;
    end if;
  end Degree;

  procedure Test_Binomial_Solver
               ( A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

    v : constant Standard_Integer_Vectors.Vector(A'range(1)) := Tropism(A);
    p : constant integer32 := Pivot(v);
    M : Matrix(A'range(1),A'range(1));
    B : Matrix(A'range(1),A'range(2));
    D : Matrix(A'first(1)..A'last(1)-1,A'range(2));

  begin
    put("A tropism : "); put(v); new_line;
    put("The residual : "); put(Transpose(A)*v); new_line;
    put("Degree of the space curve : "); put(Degree(v),1); new_line;
    if p >= v'first then
      M := Eliminate(v,p);
      put_line("The transformation matrix : "); put(M);
      put("-> its determinant : "); put(det(M),1); new_line;
      B := M*A;
      put_line("After applying the transformation : "); put(B);
      for i in D'range(1) loop
        for j in D'range(2) loop
          D(i,j) := B(i+1,j);
        end loop;
      end loop;
      put_line("After truncating the first row : "); put(D);
      put("Determinant : "); put(det(D),1); new_line;
    end if;
  end Test_Binomial_Solver;

  procedure Transform ( M : in Matrix; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Performs a unimodular coordinate transformation M to p.

    t : constant Transfo := Create(M);
    q : constant Laur_Sys(p'range) := Polynomial_to_Laurent_System(p);
    tq : constant Laur_Sys(q'range) := Transform(t,q);
    tp : constant Poly_Sys(p'range) := Laurent_to_Polynomial_System(tq);
    file : file_type;

  begin
    put_line("The transformed system : "); put(tq);
    new_line;
    skip_line;
    Read_Name_and_Create_File(file);
    put(file,tp);
    close(file);
  end Transform;

  procedure Test_Unimodular_Transformation
              ( nq,nv : in integer32; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Prompts the user for a direction to define a unimodular transformation.

    v : Standard_Integer_Vectors.Vector(1..nv);
    piv : integer32;
    M : Matrix(1..nv,1..nv);

  begin
    new_line;
    put("Give "); put(nv,1); put(" integers for a tropism : ");
    get(v);
    put("Your tropism : "); put(v);
    piv := Pivot(v);
    put(" its pivot is "); put(piv,1); new_line;
    M := Eliminate(v,piv);
    put_line("The corresponding unimodular transformation : "); put(M);
    Transform(M,p);
  end Test_Unimodular_Transformation;

  procedure Test_Solver ( p : in Poly_Sys ) is

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    fail : boolean;
    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    c : Standard_Complex_Vectors.Vector(1..nq);
 
  begin
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    Parse(p,nq,A,c,fail);
    if fail then
      put_line("The system is not a binomial system.");
      Test_Unimodular_Transformation(nq,nv,p);
    else
      put_line("The system is a binomial system.");
      put_line("The exponent matrix :"); put(A);
      put_line("The coefficient vector :"); put_line(c);
      Test_Binomial_Solver(A,c);
    end if;
  end Test_Solver;

  function Create_Random_Binomial_System
             ( nq,nv : integer32 ) return Poly_Sys is

    res : Poly_Sys(1..nq);
    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    c : constant Standard_Complex_Vectors.Vector(1..nq)
      := Random_Vector(1,nq);
    lower,upper : integer32 := 0;

  begin
    put("  give lower bound for exponents : "); get(lower);
    put("  give upper bound for exponents : "); get(upper);
    A := Random_Matrix(natural32(nv),natural32(nq),lower,upper);
    put_line("The exponent matrix : "); put(A);
    put_line("The coefficient vector : "); put_line(c);
    res := Create(A,c);
    put_line("The binomial system : "); put(res);
    return res;
  end Create_Random_Binomial_System;

  procedure Read_and_Test_Binomial_Solver ( n : in integer32 ) is

    A : Standard_Integer_Matrices.Matrix(1..n,1..n-1);
    c : Standard_Complex_Vectors.Vector(1..n-1);

  begin
    put("Give an "); put(n,1); put("-by-"); put(n-1,1);
    put_line(" exponent matrix : "); get(A);
    put("Give an "); put(n-1,1);
    put_line("-vector for righthand side :"); get(c);
    put_line("The exponent matrix is "); put(A);
    put_line("The coefficient vector is "); put_line(c);
    Test_Binomial_Solver(A,c);
  end Read_and_Test_Binomial_Solver;
 
  procedure Main is

    lp : Link_to_Poly_Sys;
    ans : character;
    nq,nv : integer32 := 0;

  begin
    new_line;
    put_line("Computing Puiseux series for binomial systems...");
    new_line;
    put_line("MENU for testing the operations on Puiseux series :");
    put_line("  1. create and solve a random binomial system;");
    put_line("  2. give a polynomial system to test the solver;");
    put_line("  3. give a system x^A = c to test the solver.");
    put("Type 1, 2, or 3 to select : "); 
    Ask_Alternative(ans,"123");
    new_line;
    if ans = '3' then
      put("Give the number of variables : "); get(nv);
      Read_and_Test_Binomial_Solver(nv);
    else
      if ans = '1' then
        put_line("Creating a random binomial system...");
        put("  give the number of equations : "); get(nq);
       -- put("  give the number of variables : "); get(nv);
        nv := nq + 1;
        put("-> the number of variables is "); put(nv,1); new_line;
        lp := new Poly_Sys'(Create_Random_Binomial_System(nq,nv));
      else
        put_line("Reading a polynomial system...");
        get(lp);
      end if;
      Test_Solver(lp.all);
    end if;
  end Main;

begin
  Main;
end ts_binpser;
