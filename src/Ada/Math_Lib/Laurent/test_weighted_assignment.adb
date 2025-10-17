with Ada.Text_io;                       use Ada.Text_IO;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;        use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Double_Weighted_Assignment;

package body Test_Weighted_Assignment is

  procedure Check_Matching
              ( A : in Matrix; m : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Checks the value of the matching m, computing the value,
  --   which should equal the minimum cost.
  --
    val : double_float := 0.0;

  begin
    put_line("Checking the value of the matching ...");
    put_line("A : "); put(A,1,A'last(1),0);
    put("m : "); put(m,1); new_line;
    for i in A'range(2) loop
      put("("); put(m(i),1); put(","); put(i,1);
      put(") : "); put(A(i,m(i))); new_line;
      val := val + A(i,m(i));
    end loop;
    put("value : "); put(val); new_line;
  end Check_Matching;

  procedure Run_Brute_Force
              ( A : in Matrix; val : out double_float;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Runs the brute force method, enumerating all columns of A.
  --   The minimum cost is returned in val.

    columns,selection : Standard_Integer_Vectors.Vector(A'range(2));
    accu : constant double_float := 0.0;
    mincost : double_float := 1.0E+99;

  begin
    columns := (columns'range => 0);
    selection := (selection'range => 0);
    Double_Weighted_Assignment.enumerate
      (A,1,columns,accu,selection,mincost,vrblvl);
    put("Minimum : "); put(mincost,1,3,3);
    put(" at"); put(selection); new_line;
    Check_Matching(A,selection);
    val := Double_Weighted_Assignment.Value_Selection(A,selection);
    if val = mincost
     then put_line("The value matches the minimum cost, okay.");
     else put_line("The value does not match the minimum cost, bug!");
    end if;
  end Run_Brute_Force;

  procedure Run_Hungarian
              ( A : in Matrix; val : out double_float;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Runs the Hungarian algorithm on the matrix A.
  --   The minimum cost is returned in val.

    u : Standard_Floating_Vectors.Vector(0..A'last(1));
    v : Standard_Floating_Vectors.Vector(0..A'last(2));
    p : Standard_Integer_Vectors.Vector(0..A'last(2));
    w : Standard_Integer_Vectors.Vector(1..A'last(2));
    match : Standard_Integer_Vectors.Vector(A'range(1));
    mincost : double_float := 0.0;

  begin
    Double_Weighted_Assignment.Hungarian(A,u,v,p,w,vrblvl);
    for i in 1..u'last loop
      mincost := mincost + u(i);
    end loop;
    for i in 1..v'last loop
      mincost := mincost + v(i);
    end loop;
    match := Double_Weighted_Assignment.Row_Matching(p,A'last(1));
    put("Minimum : "); put(mincost,1,3,3);
    put(" at"); put(match); new_line;
    Check_Matching(A,match);
    val := Double_Weighted_Assignment.Value_Selection(A,match);
    if val = mincost
     then put_line("The value matches the minimum cost, okay.");
     else put_line("The value does not match the minimum cost, bug!");
    end if;
  end Run_Hungarian;

  procedure Run_Test_Example is

    A : Matrix(1..5,1..5);
    v1,v2 : double_float;

  begin
    -- 7 2 1 9 4
    A(1,1) := 7.0; A(1,2) := 2.0; A(1,3) := 1.0; A(1,4) := 9.0; A(1,5) := 4.0;
   -- 9 6 9 5 5
    A(2,1) := 9.0; A(2,2) := 6.0; A(2,3) := 9.0; A(2,4) := 5.0; A(2,5) := 5.0;
   -- 3 8 3 1 8
    A(3,1) := 3.0; A(3,2) := 8.0; A(3,3) := 3.0; A(3,4) := 1.0; A(3,5) := 8.0;
   -- 7 9 4 2 2
    A(4,1) := 7.0; A(4,2) := 9.0; A(4,3) := 4.0; A(4,4) := 2.0; A(4,5) := 2.0;
   -- 8 4 7 4 8
    A(5,1) := 8.0; A(5,2) := 4.0; A(5,3) := 7.0; A(5,4) := 4.0; A(5,5) := 8.0;
    put_line("The input matrix : "); put(A,1,5,0);
    put_line("-> running the brute force method ...");
    Run_Brute_Force(A,v1,1);
    put_line("-> running the Hungarian algorithm ...");
    Run_Hungarian(A,v2,1);
    if v1 = v2
     then put_line("Brute force output = Hungarian output, okay.");
     else put_line("Brute force output /= Hungarian output, bug!");
    end if;
  end Run_Test_Example;

  function Random_Matrix ( dim : integer32 ) return Matrix is

  -- DESCRIPTION :
  --   Returns a random matrix of dimension dim,
  --   with random integers in the range 1 to 9.

    res : Matrix(1..dim,1..dim);
    random_digit : integer32;

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        random_digit := Standard_Random_Numbers.Random(1,9);
        res(i,j) := double_float(random_digit);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Vector
             ( dim : integer32 )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a random vector of dimension dim,
  --   with random integers in the range 1 to 9.

    res : Standard_Floating_Vectors.Vector(1..dim);
    random_digit : integer32;

  begin
    for i in 1..dim loop
      random_digit := Standard_Random_Numbers.Random(1,9);
      res(i) := double_float(random_digit);
    end loop;
    return res;
  end Random_Vector;

  procedure Run_Random_Example ( dim : in integer32 ) is

    A : constant Matrix(1..dim,1..dim) := Random_Matrix(dim);
    v1,v2 : double_float;

  begin
    put_line("The input matrix : "); put(A,1,dim,0);
    put_line("-> running the brute force method ...");
    Run_Brute_Force(A,v1,1);
    put_line("-> running the Hungarian algorithm ...");
    Run_Hungarian(A,v2,1);
    if v1 = v2
     then put_line("Brute force output = Hungarian output, okay.");
     else put_line("Brute force output /= Hungarian output, bug!");
    end if;
  end Run_Random_Example;

  procedure Test_Cramer_Vector ( dim : in integer32 ) is

    A : constant Matrix(1..dim,1..dim) := Random_Matrix(dim);
    b : constant Standard_Floating_Vectors.Vector(1..dim)
      := Random_Vector(dim);
    c : Standard_Floating_Vectors.Vector(0..dim);
    m : Standard_Integer_VecVecs.VecVec(0..dim);

  begin
    put_line("A random matrix A : "); put(A,1,dim,0);
    put_line("A random vector b : "); put(b,0); new_line;
    for i in m'range loop
      m(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    Double_Weighted_Assignment.cramer_vector(A,b,c,m,1);
    for i in 0..dim loop
      put(i,1); put(" : "); put(c(i),0);
      put(" :"); put(m(i).all); new_line;
    end loop;
  end Test_Cramer_Vector;

  procedure Main is

    ans : character;
    dim : integer32 := 0;

  begin
    new_line;
    put_line("MENU to test the weighted assignment solvers :");
    put_line("  0. test a textbook example");
    put_line("  1. generate a random example to test");
    put_line("  2. compute a Cramer vector on random data");
    put("Type 0, 1, or 2 to select a test : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' =>
        put_line("-> testing textbook example ...");
        Run_Test_Example;
      when others =>
        put("-> give the dimension of the matrix : "); get(dim);
        if ans = '1'
         then Run_Random_Example(dim);
         else Test_Cramer_Vector(dim);
        end if;
    end case;
  end Main;

end Test_Weighted_Assignment;
