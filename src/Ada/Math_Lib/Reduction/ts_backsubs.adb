with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Double_Double_Vectors;
with Double_Double_Matrices;

procedure ts_backsubs is

-- DESCRIPTION :
--   Test on the accuracy of the back substitution on a random 
--   floating-point matrix of positive numbers in [-1, +1].
--   In standard double precision, problems may occur already at dimensions
--   as small as 40, while there is no problem in double double precision.

  function Random_Upper_Triangular_Matrix
             ( dim : integer32 ) return Standard_Floating_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a random dim-by-dim upper triangular matrix 
  --   with numbers in [0, 1].

    res : Standard_Floating_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        if i > j then
          res(i,j) := 0.0;
        else
          res(i,j) := Standard_Random_Numbers.Random;
         -- if res(i,j) < 0.0
         --  then res(i,j) := -res(i,j);
         -- end if;
        end if;
      end loop;
    end loop;
    return res;
  end Random_Upper_Triangular_Matrix;

  function Random_Upper_Triangular_Matrix
             ( dim : integer32 ) return Double_Double_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a random dim-by-dim upper triangular matrix 
  --   with numbers in [0, 1].

    res : Double_Double_Matrices.Matrix(1..dim,1..dim);

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        if i > j then
          res(i,j) := Double_Double_Numbers.Create(0.0);
        else
          res(i,j) := DoblDobl_Random_Numbers.Random;
         -- if res(i,j) < 0.0
         --  then res(i,j) := -res(i,j);
         -- end if;
        end if;
      end loop;
    end loop;
    return res;
  end Random_Upper_Triangular_Matrix;

  function Solve_Upper_Triangular_System
             ( mat : Standard_Floating_Matrices.Matrix;
               rhs : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Solves the square upper triangular system defined by the matrix mat
  --   and right hand vector in rhs.

    res : Standard_Floating_Vectors.Vector(rhs'range);

  begin
    for i in reverse res'range loop
      res(i) := rhs(i);
      for j in i+1..res'last loop
        res(i) := res(i) - mat(i,j)*res(j);
      end loop;
      res(i) := res(i)/mat(i,i);
    end loop;
    return res;
  end Solve_Upper_Triangular_System; 

  function Solve_Upper_Triangular_System
             ( mat : Double_Double_Matrices.Matrix;
               rhs : Double_Double_Vectors.Vector )
             return Double_Double_Vectors.Vector is

  -- DESCRIPTION :
  --   Solves the square upper triangular system defined by the matrix mat
  --   and right hand vector in rhs.

    res : Double_Double_Vectors.Vector(rhs'range);

  begin
    for i in reverse res'range loop
      res(i) := rhs(i);
      for j in i+1..res'last loop
        res(i) := res(i) - mat(i,j)*res(j);
      end loop;
      res(i) := res(i)/mat(i,i);
    end loop;
    return res;
  end Solve_Upper_Triangular_System; 

  procedure Standard_Back_Substitution ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random upper triangular matrix with positve numbers
  --   and a right hand side vector so the (1, 1, ..., 1) vector is
  --   an exact solution of the system.  Back substitution in plain
  --   double precision arithmetic is applied to solve the system.

    mat : constant Standard_Floating_Matrices.Matrix(1..dim,1..dim)
        := Random_Upper_Triangular_Matrix(dim);
    sol : constant Standard_Floating_Vectors.Vector(1..dim) := (1..dim => 1.0);
    use Standard_Floating_Matrices;
    rhs : Standard_Floating_Vectors.Vector(1..dim) := mat*sol;
    num : Standard_Floating_Vectors.Vector(1..dim) 
        := Solve_Upper_Triangular_System(mat,rhs);
    err : double_float;

  begin
   -- put_line("The solution :"); put_line(sol);
   -- put_line("The matrix :"); put(mat,3);
   -- put_line("The right hand side :"); put_line(rhs);
    put(dim,4);
   -- put_line("The first component of the computed solution vector : ");
    put("   ");
    put(num(num'first)); -- new_line;
   -- put("error : "); 
    put("   ");
    err := abs(1.0 - num(num'first));
    put(err,3); new_line;
  end Standard_Back_Substitution;

  procedure DoblDobl_Back_Substitution ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random upper triangular matrix with positve numbers
  --   and a right hand side vector so the (1, 1, ..., 1) vector is
  --   an exact solution of the system.  Back substitution in double
  --   double precision arithmetic is applied to solve the system.

    mat : constant Double_Double_Matrices.Matrix(1..dim,1..dim)
        := Random_Upper_Triangular_Matrix(dim);
    sol : constant Double_Double_Vectors.Vector(1..dim)
        := (1..dim => Double_Double_Numbers.Create(1.0));
    use Double_Double_Matrices;
    rhs : Double_Double_Vectors.Vector(1..dim) := mat*sol;
    num : Double_Double_Vectors.Vector(1..dim) 
        := Solve_Upper_Triangular_System(mat,rhs);
    err : double_double;

  begin
   -- put_line("The solution :"); put_line(sol);
   -- put_line("The matrix :"); put(mat,3);
   -- put_line("The right hand side :"); put_line(rhs);
    put(dim,4);
   -- put_line("The first component of the computed solution vector : ");
    put("   ");
    put(num(num'first)); -- new_line;
   -- put("error : "); 
    put("   ");
    err := abs(1.0 - num(num'first));
    put(err,3); new_line;
  end DoblDobl_Back_Substitution;

  procedure Main is

    dim : integer32 := 0;

  begin
    new_line;
    put_line("Testing the back substitution ...");
   -- put("Give the dimension : "); get(dim);
    dim := 8;
    put_line("  dim    first component    error  ");
    for k in 1..10 loop
     -- put("Dimension : "); put(dim,1); put_line(" :");
      Standard_Back_Substitution(dim);
      dim := dim + 8;
    end loop;
    dim := 8;
    put("  dim");
    put("                 first component");
    put_line("              error  ");
    for k in 1..10 loop
     -- put("Dimension : "); put(dim,1); put_line(" :");
      DoblDobl_Back_Substitution(dim);
      dim := dim + 64;
    end loop;
  end Main;

begin
  Main;
end ts_backsubs;
