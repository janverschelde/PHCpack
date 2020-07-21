with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Standard_Vector_Splitters;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;
with Standard_Matrix_Splitters;
with Standard_Complex_Singular_Values;  use Standard_Complex_Singular_Values;

procedure ts_perfdsvd is

-- DESCRIPTION :
--   Test the development of a better performing SVD implementation,
--   in double precision.

  procedure Test ( n,p : in integer32;
                   A : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in A is an n-by-p matrix.

    wrk : Standard_Complex_Matrices.Matrix(1..n,1..p) := A;
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    S : Standard_Complex_Vectors.Vector(1..mm);
    e : Standard_Complex_Vectors.Vector(1..p);
    U : Standard_Complex_Matrices.Matrix(1..n,1..n);
    V : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    SVD(wrk,n,p,S,e,U,V,job,info);
    put_line("The singular values : "); put_line(S);
    put_line("The matrix U :"); put(U,2);
    put_line("The matrix V :"); put(V,2);
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns
  --   and then generates a random complex matrix
  --   with that number of rows and columns.

    nrows,ncols : integer32 := 0;

  begin
    put("Give the number of rows : "); get(nrows);
    put("Give the number of columns : "); get(ncols);
    declare
      A : constant Standard_Complex_Matrices.Matrix(1..nrows,1..ncols)
        := Standard_Random_Matrices.Random_Matrix(1,nrows,1,ncols);
      rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
    begin
      put("A random "); put(nrows,1); put("-by-"); put(ncols,1);
      put_line(" complex matrix :");
      put(A,2);
      rvv := Standard_Vector_Splitters.Allocate(ncols,nrows,1,1);
      ivv := Standard_Vector_Splitters.Allocate(ncols,nrows,1,1);
      Standard_Matrix_Splitters.Complex_Parts(A,rvv,ivv);
      put_line("The real parts of the columns :"); put(rvv,2);
      put_line("The imaginary parts of the columns :"); put(ivv,2);
      Test(nrows,ncols,A);
    end;
  end Main;

begin
  Main;
end ts_perfdsvd;
