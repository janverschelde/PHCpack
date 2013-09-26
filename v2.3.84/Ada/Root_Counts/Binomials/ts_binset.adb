with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Multprec_Integer_Numbers;
with Multprec_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers; 
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io; 
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Integer_Linear_Solvers;
with Standard_Integer_Kernel;
with Standard_Integer64_Matrices;
with Multprec_Integer_Matrices;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Random_Matrices;          use Standard_Random_Matrices;
with Symbol_Table;
with Standard_Complex_Laurentials;      use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;     use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Exponent_Transformations;
with Standard_Binomial_Systems;         use Standard_Binomial_Systems;
with Standard_Binomial_Varieties;
with Standard_Binomial_Varieties_io;
with Multprec_Binomial_Varieties;

procedure ts_binset is

-- DESCRIPTION :
--   Interactive development of binomial system solver, to write an
--   explicit description of its positive dimensional solution set.

  procedure Show_Pivots ( V : in Standard_Integer_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   The pivots of V are the elements on the diagonal of the
  --   Hermite normal form of the matrix V.

    W : Standard_Integer_Matrices.Matrix(V'range(2),V'range(1))
      := Standard_Integer_Matrices.Transpose(V);
    rank : integer32;
    pivots : Standard_Integer_Vectors.Vector(W'range(1));

  begin
    Standard_Integer_Linear_Solvers.Upper_Triangulate(W);
    Standard_Integer_Kernel.Pivots_in_Upper(W,rank,pivots);
    put("Pivots of Hermite form : ");
    for i in pivots'range loop
      put(" "); put(W(i,pivots(i)),1);
    end loop;
    new_line;
  end Show_Pivots;

  procedure Write_Solution_to_File
               ( d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and then writes the algebraic
  --   set of dimension d, defined by M and c to that file.
 
    file : file_type;

  begin
    Read_Name_and_Create_File(file);
    Standard_Binomial_Varieties_io.Write_Header(file,natural32(M'last(1)),d);
    Standard_Binomial_Varieties_io.Write_Solution(file,d,M,c);
  end Write_Solution_to_File;

  procedure Write_Solution_to_File
               ( d : in natural32;
                 M : in Standard_Integer_Matrices.Matrix;
                 w : in Standard_Integer_Vectors.Vector;
                 c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and then writes the algebraic
  --   set of dimension d, defined by M and c to that file.
 
    file : file_type;

  begin
    Read_Name_and_Create_File(file);
    Standard_Binomial_Varieties_io.Write_Header(file,natural32(M'last(1)),d);
    Standard_Binomial_Varieties_io.Write_Solution(file,d,M,w,c);
  end Write_Solution_to_File;
 
  procedure Write_Solutions_to_File
               ( d : in natural32; M : in Standard_Integer_Matrices.Matrix;
                 c : in Solution_List ) is

    tmp : Solution_List := c;

  begin
    for i in 1..Length_Of(c) loop
      put("writing solution "); put(i,1); put_line(" to file ...");
      Write_Solution_to_File(d,M,Head_Of(tmp).v);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Solutions_to_File;

  procedure Write_Solutions_to_File
               ( d : in natural32; M : in Standard_Integer_Matrices.Matrix;
                 w : in Standard_Integer_Vectors.Vector;
                 c : in Solution_List ) is

    tmp : Solution_List := c;

  begin
    put("the denominators are "); put(w); new_line;
    put_line("the coordinate transformation is "); put(M);
    -- filling the symbol table is done in Random_Point_Filter...
    Standard_Binomial_Varieties_io.Fill_Symbol_Table(natural32(M'last(1)));
    for i in 1..Length_Of(c) loop
      put("writing solution "); put(i,1); put_line(" to file ...");
      Write_Solution_to_File(d,M,w,Head_Of(tmp).v);
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Solutions_to_File;

  procedure Test_Standard_Integer_Binomial_Solver
               ( A : in Standard_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes a basis for the kernel of the exponent matrix A
  --   using standard 32-bit integer arithmetic.

    r,d : integer32;
    V : Standard_Integer_Matrices.Link_to_Matrix;
    bug : boolean;
    ans : character;

  begin
    Standard_Binomial_Varieties.Cone_of_Tropisms(A,r,V);
    put("rank = "); put(r,1); new_line;
    if r = A'last(1) then
      put_line("The exponent matrix is of full rank, no tropisms.");
    else
      Standard_Binomial_Varieties.Check_Cone(standard_output,A,V.all,r,bug);
      if bug
       then put_line("Bug in computation of cone of tropisms.");
       else put_line("Computation of cone of tropisms is okay.");
      end if;
      Standard_Binomial_Varieties.Expected_Dimension
        (standard_output,A,V.all,r,d);
      if d > 0 then
       -- put("Extra output during computations ? (y/n) ");
       -- Ask_Yes_or_No(ans);
        Show_Pivots(V.all);
        declare
          M : Standard_Integer_Matrices.Matrix(V'range(1),V'range(1));
          w : Standard_Integer_Vectors.Vector(V'range(2));
          cff : Solution_List;
        begin
          Standard_Binomial_Varieties.Solve
            (standard_output,d,A,V.all,c,1.0E-8,M,w,cff);
          if not Is_Null(cff) then
            put("Write solutions to file ? (y/n) ");
            Ask_Yes_or_No(ans);
            if ans = 'y' then
              if w(w'first) = 0 
               then Write_Solutions_to_File(natural32(d),M,cff);
               else Write_Solutions_to_File(natural32(d),M,w,cff);
              end if;
            end if;
          end if;
        end;
      end if;
    end if;
  end Test_Standard_Integer_Binomial_Solver;

  procedure Test_Standard_Integer64_Binomial_Solver
               ( A : in Standard_Integer64_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes a basis for the kernel of the exponent matrix A
  --   using standard 64-bit integer arithmetic.

    r,d : integer32;
    V : Standard_Integer64_Matrices.Link_to_Matrix;
    bug : boolean;
    use Standard_Exponent_Transformations;

  begin
    Standard_Binomial_Varieties.Cone_of_Tropisms(A,r,V);
    Standard_Binomial_Varieties.Check_Cone(standard_output,A,V.all,r,bug);
    if bug
     then put_line("Bug in computation of cone of tropisms.");
     else put_line("Computation of cone of tropisms is okay.");
    end if;
    Standard_Binomial_Varieties.Expected_Dimension
      (standard_output,A,V.all,r,d);
    if d > 0 then
      declare
        M : Standard_Integer64_Matrices.Matrix(V'range(1),V'range(1));
        U : Standard_Integer64_Matrices.Matrix(1..M'last(1)-d,A'range(2));
        fail : boolean;
        nzr : natural32;
        rank : integer32;
        pivots : Standard_Integer_Vectors.Vector(c'range);
        pp : Standard_Integer_Numbers.integer64;
        use Standard_Binomial_Varieties;
      begin
        Unimodular_Coordinate_Transformation(standard_output,V.all,fail,M);
        if not fail then
          Test_Unimodular_Coordinate_Transformation
            (standard_output,A,M,d,nzr,fail);
          pivots := (pivots'range => 0);
          Upper_Transformed_Exponents(standard_output,A,M,d,U,rank,pivots);
          pp := Product_of_Pivots(U,pivots);
          put("The product of pivots : ");
          Standard_Integer_Numbers_io.put(pp,1); new_line;
        end if;
      end;
    end if;
  end Test_Standard_Integer64_Binomial_Solver;

  procedure Test_Multprec_Integer_Binomial_Solver
               ( A : in Multprec_Integer_Matrices.Matrix;
                 c : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes a basis for the kernel of the exponent matrix A
  --   using multiprecision integer arithmetic.

    r,d : integer32;
    V : Multprec_Integer_Matrices.Link_to_Matrix;
    bug : boolean;

  begin
    Multprec_Binomial_Varieties.Cone_of_Tropisms(A,r,V);
    Multprec_Binomial_Varieties.Check_Cone(A,V.all,r,true,bug);
    if bug
     then put_line("Bug in computation of cone of tropisms.");
     else put_line("Computation of cone of tropisms is okay.");
    end if;
    Multprec_Binomial_Varieties.Expected_Dimension(A,V.all,r,true,d);
    if d > 0 then
      declare
        M : Multprec_Integer_Matrices.Matrix(V'range(1),V'range(1));
        U : Multprec_Integer_Matrices.Matrix(1..M'last(1)-d,A'range(2));
        fail : boolean;
        nzr : natural32;
        rank : integer32;
        pivots : Standard_Integer_Vectors.Vector(c'range);
        pp : Multprec_Integer_Numbers.Integer_Number;
        use Multprec_Binomial_Varieties;
      begin
        Unimodular_Coordinate_Transformation(V.all,true,fail,M);
        if not fail then
          Test_Unimodular_Coordinate_Transformation(A,M,d,true,nzr,fail);
          pivots := (pivots'range => 0);
          Upper_Transformed_Exponents(A,M,d,true,U,rank,pivots);
          pp := Product_of_Pivots(U,pivots);
          put("The product of pivots : ");
          Multprec_Integer_Numbers_io.put(pp); new_line;
        end if;
      end;
    end if;
  end Test_Multprec_Integer_Binomial_Solver;

  function to_Standard_Integer64
               ( A : Standard_Integer_Matrices.Matrix )
               return Standard_Integer64_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the 64-bit version of the given
  --   standard 32-bit integer matrix.
  
    res : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2));
  
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Standard_Integer_Numbers.integer64(A(i,j));
      end loop;
    end loop;
    return res;
  end to_Standard_Integer64;

  function standard_to_multprec
               ( A : Standard_Integer_Matrices.Matrix )
               return Multprec_Integer_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the multiprecision version of the given
  --   standard 32-bit integer matrix.
  
    res : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
  
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := Multprec_Integer_Numbers.Create(A(i,j));
      end loop;
    end loop;
    return res;
  end standard_to_multprec;

  procedure Test_Solver ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Parses the system p into a binomial system, extracting
  --   exponent matrix and coefficient vector, and then calls
  --   the binomial system solver if p is a binomial system.

    nq : constant integer32 := p'last;
    nv : constant integer32 := integer32(Number_of_Unknowns(p(p'first)));
    fail : boolean;
    A : Standard_Integer_Matrices.Matrix(1..nv,1..nq);
    c : Standard_Complex_Vectors.Vector(1..nq);
    ans : character;
 
  begin
    put("number of equations : "); put(nq,1); new_line;
    put("number of variables : "); put(nv,1); new_line;
    Parse(p,nq,A,c,fail);
    if fail then
      put_line("The system is not a binomial system.");
    else
      put_line("The system is a binomial system.");
      put_line("The exponent matrix :"); put(A);
      put("Do you want to see the coefficient vector ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put_line("The coefficient vector :"); put_line(c);
      end if;
      put("Use 64-bit arithmetic on the exponent matrix ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        declare
          B : constant Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2))
            := to_Standard_Integer64(A);
        begin
          Test_Standard_Integer64_Binomial_Solver(B,c);
        end;
      else
        put("Use multiprecision on the exponent matrix ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          declare
            B : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2))
              := standard_to_multprec(A);
          begin
            Test_Multprec_Integer_Binomial_Solver(B,c);
            Multprec_Integer_Matrices.Clear(B);
          end;
        else
          Test_Standard_Integer_Binomial_Solver(A,c);
        end if;
      end if;
    end if;
  end Test_Solver;

  function Create_Random_Binomial_System
             ( nq,nv : integer32 ) return Laur_Sys is

  -- DESCRIPTION :
  --   Returns a binomial system of nq equations in nv variables,
  --   with randomly generated exponent matrix and coefficient vector.

    res : Laur_Sys(1..nq);
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

  -- DESCRIPTION :
  --   Reads an exponent matrix and coefficient vector
  --   of a binomial system and calls the binomial system solver.

    A : Standard_Integer_Matrices.Matrix(1..n,1..n);
    c : Standard_Complex_Vectors.Vector(1..n);

  begin
    put("Give an "); put(n,1); put("-by-"); put(n,1);
    put_line(" exponent matrix : "); get(A);
    put("Give an "); put(n,1);
    put_line("-vector for righthand side :"); get(c);
    put_line("The exponent matrix is "); put(A);
    put_line("The coefficient vector is "); put_line(c);
    Test_Standard_Integer_Binomial_Solver(A,c);
  end Read_and_Test_Binomial_Solver;

  function Evaluate_Algebraic_Set
              ( A,T : Standard_Integer_Matrices.Matrix;
                b,c,x : Standard_Complex_Vectors.Vector ) 
              return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns y^A - b for y = c*x^T, where T is the n-by-d matrix
  --   of tropisms and x contains values for the d parameters.

    y : Standard_Complex_Vectors.Vector(T'range(1))
      := Standard_Binomial_Varieties.Evaluate_Tropisms(T,x);

  begin
    for i in y'range loop
      y(i) := c(i)*y(i);
    end loop;
    return Standard_Binomial_Varieties.Evaluate_Binomial_System(A,b,y);
  end Evaluate_Algebraic_Set;

  procedure Random_Point_Test
              ( A : in Standard_Integer_Matrices.Matrix;
                b : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Reads a solution for the binomial system x^A = b,
  --   and evaluates a random point at the solution in the system.

    file : file_type; 
    sol : Link_to_Laur_Sys;
    n,d : integer32;
    T : Standard_Integer_Matrices.Link_to_Matrix;
    c : Standard_Complex_Vectors.Link_to_Vector;

  begin
    new_line;
    put_line("Reading a solution to a binomial system...");
    Read_Name_and_Open_File(file);
    Symbol_Table.Clear;
    get(file,sol);
    put_line("The solution : "); put(sol.all);
    n := sol'last;
    put("number of equations : "); put(n,1); new_line;
    d := integer32(Number_of_Unknowns(sol(sol'first))) - n;
    put("dimension of the set : "); put(d,1); new_line;
    Standard_Binomial_Varieties_io.Parse_Binomial_Variety(sol.all,T,c);
    put_line("The tropisms : "); put(T.all);
    put_line("The coefficients : "); put_line(c.all);
    declare
      x : constant Standard_Complex_Vectors.Vector(1..d) := Random_Vector(1,d);
      z : constant Standard_Complex_Vectors.Vector(b'range)
        := Evaluate_Algebraic_Set(A,T.all,b,c.all,x);
      res : constant double_float := Standard_Complex_Norms_Equals.Max_Norm(z);
    begin
      put_line("after evaluation : "); put_line(z);
      put("The max norm : "); put(res); new_line;
    end;
  end Random_Point_Test;

  procedure Random_Point_Test ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Parses the polynomial system p into a binomial system
  --   and calls the Random_Point_Test routine.
 
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
      put_line("The system is not a binomial system!");
    else
      put_line("The exponent matrix : "); put(A);
      put_line("The coefficient vector : "); put_line(c);
      Random_Point_Test(A,c);
    end if;
  end Random_Point_Test;

  procedure Compute_Degree ( s : in Laur_Sys ) is

  -- DESCRIPTION :
  --   On input is a solution to a binomial system,
  --   printed is the degree of this solution set.

    n : constant integer32 := s'last;
    d : constant integer32 := integer32(Number_of_Unknowns(s(s'first))) - n;
    T : Standard_Integer_Matrices.Link_to_Matrix;
    c : Standard_Complex_Vectors.Link_to_Vector;
    v : natural32;

  begin
    put_line("The solution : "); put(s);
    put("number of original variables : "); put(n,1); new_line;
    put("dimension of the set : "); put(d,1); new_line;
    Standard_Binomial_Varieties_io.Parse_Binomial_Variety(s,T,c);
    put_line("The tropisms : "); put(T.all);
    put_line("The coefficients : "); put_line(c.all);
    v := Standard_Binomial_Varieties.Degree(T.all);
    put("the degree is the volume : "); put(v,1); new_line;
  end Compute_Degree;

  procedure Call_Black_Box_Solver ( p : in Laur_Sys ) is

  -- DESCRIPTION :
  --   Calls the black box solver.

    fail : boolean;
    d : integer32;
    M : Standard_Integer_Matrices.Link_to_Matrix;
    cff,tmp : Solution_List;
    ans : character;
    use Standard_Binomial_Varieties;
 
  begin
    put("Do you want output ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Black_Box_Solver(standard_output,p,fail,d,M,cff);
     else Black_Box_Solver(p,fail,d,M,cff);
    end if;
    if not Is_Null(cff) then
      tmp := cff;
      for i in 1..Length_Of(cff) loop
        put("solution "); put(i,1); put_line(" : ");
        Standard_Binomial_Varieties_io.Write_Header
          (standard_output,natural32(M'last(1)),natural32(d));
        Standard_Binomial_Varieties_io.Write_Solution
          (standard_output,natural32(d),M.all,Head_Of(tmp).v);
        tmp := Tail_of(tmp);
      end loop;
      put("Write solutions to file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' 
       then Write_Solutions_to_File(natural32(d),M.all,cff);
      end if;
    end if;
  end Call_Black_Box_Solver;

  procedure Main is

    lp : Link_to_Laur_Sys;
    ans : character;
    nq,nv : integer32 := 0;

  begin
    new_line;
    put_line("Positive dimensional toric solution sets of binomial systems.");
    new_line;
    put_line("MENU for testing the operations on binomial systems :");
    put_line("  1. create and solve a random binomial system;");
    put_line("  2. give a polynomial system to test the solver;");
    put_line("  3. give a system x^A = b to test the solver;");
    put_line("  4. do random point test on system x^A = b and solution;");
    put_line("  5. compute the degree of a solution to x^A = b;");
    put_line("  6. black box solver on given polynomial system;");
    put("Type 1, 2, 3, 4, 5, or 6 to select : "); 
    Ask_Alternative(ans,"123456");
    new_line;
    case ans is
      when '1' =>
        put_line("Creating a random binomial system...");
        put("  give the number of equations : "); get(nq);
        put("  give the number of variables : "); get(nv);
        lp := new Laur_Sys'(Create_Random_Binomial_System(nq,nv));
      when '2' =>
        put_line("Reading a binomial system...");
        get(lp);
      when '3' =>
        put("Give the dimension of the system : "); get(nq);
        Read_and_Test_Binomial_Solver(nq);
      when '4' =>
        put_line("Reading a binomial system..."); get(lp);
        Random_Point_Test(lp.all);
      when '5' =>
        put_line("Reading a solution to a binomial system..."); get(lp);
        Compute_Degree(lp.all);
      when '6' =>
        put_line("Reading a binomial system..."); get(lp);
        Call_Black_Box_Solver(lp.all);
      when others => null;    
    end case;
    if ans = '1' or ans = '2' or ans = '3'
     then Test_Solver(lp.all);
    end if;
  end Main;

begin
  Main;
end ts_binset;
