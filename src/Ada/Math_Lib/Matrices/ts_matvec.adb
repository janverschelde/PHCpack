with Ada.Text_IO;                        use Ada.Text_IO;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Random_Vectors;
with Double_Complex_MatVecs;
with Double_Real_MatVecs;

procedure ts_matvec is

-- DESCRIPTION :
--   Tests instantiations of matrices of vectors.

  procedure Write ( A : in Double_Complex_MatVecs.MatVec ) is

  -- DESCRIPTION :
  --   Writes the matrix A.
 
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("["); put(i,1); put(",");
        put(j,1); put_line("] :");
        put_line(A(i,j));
      end loop;
    end loop;
  end Write;

  procedure Write ( A : in Double_Real_MatVecs.MatVec ) is

  -- DESCRIPTION :
  --   Writes the matrix A.
 
  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("["); put(i,1); put(",");
        put(j,1); put_line("] :");
        put_line(A(i,j));
      end loop;
    end loop;
  end Write;

  procedure Test_Double_Complex ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random complex matrix of dimension dim
  --   and writes the matrix.

    cmv : Double_Complex_MatVecs.MatVec(1..dim,1..dim);

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        declare
          v : constant Standard_Complex_Vectors.Vector(1..dim)
            := Standard_Random_Vectors.Random_Vector(1,dim);
        begin
          cmv(i,j) := new Standard_Complex_Vectors.Vector'(v);
        end;
      end loop;
    end loop;
    Write(cmv);
  end Test_Double_Complex;

  procedure Test_Double_Real ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random real matrix of dimension dim
  --   and writes the matrix.

    rmv : Double_Real_MatVecs.MatVec(1..dim,1..dim);

  begin
    for i in 1..dim loop
      for j in 1..dim loop
        declare
          v : constant Standard_Floating_Vectors.Vector(1..dim)
            := Standard_Random_Vectors.Random_Vector(1,dim);
        begin
          rmv(i,j) := new Standard_Floating_Vectors.Vector'(v);
        end;
      end loop;
    end loop;
    Write(rmv);
  end Test_Double_Real;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the dimension and then generates a random instance.

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    put_line("A random matrix of double complex vectors :");
    Test_Double_Complex(dim);
    put_line("A random matrix of double real vectors :");
    Test_Double_Real(dim);
  end Main;

begin
  Main;
end ts_matvec;
