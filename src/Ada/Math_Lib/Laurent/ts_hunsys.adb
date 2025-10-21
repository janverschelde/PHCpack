with Ada.Text_io;                       use Ada.Text_IO;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Standard_Floating_Matrices;        use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Double_Weighted_Assignment;

procedure ts_hunsys is

-- DESCRIPTION :
--   Compute the leading powers and coefficients of a linear
--   system of real power series.

  function Row_Min_Plus
             ( A : Matrix; x : Standard_Floating_Vectors.Vector )
             return Standard_Floating_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the minimum of the sum of the elements on each row of A,
  --   augmented in each column with the corresponding element in x,
  --   representing the leading powers of t^A*t^x.

    res : Standard_Floating_Vectors.Vector(A'range(1));
    val : double_float;

  begin
    for i in res'range loop
      res(i) := A(i,A'first(2)) + x(x'first);
      for j in A'first(2)+1..A'last(2) loop
        val := A(i,j) + x(j);
        if val < res(i)
         then res(i) := val;
        end if;
      end loop;
    end loop;
    return res;
  end Row_Min_Plus;

  procedure Test_Random_Input ( dim : in integer32 ) is

    n : constant natural32 := natural32(dim);
    A : Matrix(1..dim,1..dim)
      := Standard_Random_Matrices.Random_Matrix(n,n);
    x : Standard_Floating_Vectors.Vector(1..dim)
      := Standard_Random_Vectors.Random_Vector(1,dim);
    b : Standard_Floating_Vectors.Vector(1..dim);
    c : Standard_Floating_Vectors.Vector(0..dim);
    d : Standard_Floating_Vectors.Vector(1..dim);
    r : Standard_Floating_Vectors.Vector(1..dim);
    m : Standard_Integer_VecVecs.VecVec(0..dim);
    abc : Matrix(1..dim,1..dim+1);
    mv : Standard_Floating_Vectors.Vector(1..dim);
    idxmv : Standard_Integer_Vectors.Vector(1..2*dim);
    idxmv1,idxmv2 : Standard_Integer_Vectors.Vector(1..dim);
    fail : boolean;

  begin
    for i in 1..dim loop
      x(i) := 1.0 + abs(x(i));
      for j in A'range(2) loop
        if i = j
         then A(i,j) := abs(A(i,j));
         else A(i,j) := 1.0 + abs(A(i,j));
        end if;
      end loop;
    end loop;
    put("A random "); put(dim,1);
    put_line("-dimensional matrix :"); put(A,3);
    put("A random "); put(dim,1);
    put_line("-dimensional vector :"); put(x,3); new_line;
    b := Row_Min_Plus(A,x);
    put_line("The right hand side vector : "); put(b,3); new_line;
    put_line("-> computing the tropical Cramer vector ...");
    for i in m'range loop
      m(i) := new Standard_Integer_Vectors.Vector'(1..dim => 0);
    end loop;
    Double_Weighted_Assignment.cramer_vector(A,b,c,m,1);
    abc := Double_Weighted_Assignment.Abc_Matrix(A,b,c);
    put_line("-> checking if minimum is attained twice ...");
    put_line("A : "); put(A,1,A'last(1),3);
    put_line("b : "); put(b,3); new_line;
    put_line("c : "); put(c,3); new_line;
    put_line("Cramer vector added to A | b :");
    for i in abc'range(1) loop
      for j in abc'first(2)..abc'last(2)-1 loop
        put(abc(i,j),3);
      end loop;
      put(" | ");
      put(abc(i,abc'last(2)),3);
      new_line;
    end loop;
    Double_Weighted_Assignment.Abc_ArgMin(abc,1.0e-10,mv,idxmv,fail,1);
    if fail
     then put_line("Failed to reach minimum exactly twice!");
    end if;
    for i in 1..dim loop
      idxmv1(i) := idxmv(2*i-1);
      idxmv2(i) := idxmv(2*i);
    end loop;
    put("1st index :"); put(idxmv1); new_line;
    put("2nd index :"); put(idxmv2); new_line;
    for i in d'range loop
      d(i) := c(i) - c(0);
    end loop;
    put_line("The leading degrees : "); put(d,3); new_line;
    put_line("The original x : "); put(x,3); new_line;
    r := Row_Min_Plus(A,d);
    put_line("Comparing right hand sides :");
    for i in r'range loop
      put("b("); put(i,1); put(") :"); put(b(i));
      put(" - r("); put(i,1); put(") :"); put(r(i));
      put(" = "); put(b(i)-r(i),2,3,3); new_line;
    end loop;
  end Test_Random_Input;

  procedure Main is

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    Test_Random_Input(dim);
  end Main;

begin
  Main;
end ts_hunsys;
