with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;       use Standard_Floating_VecVecs_io;
with Standard_Random_Vectors;
with Standard_Vector_Splitters;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;
with Standard_Matrix_Splitters;
with Standard_Complex_BLAS_Helpers;
with Standard_Complex_Singular_Values;   use Standard_Complex_Singular_Values;

procedure ts_perfdsvd is

-- DESCRIPTION :
--   Test the development of a better performing SVD implementation,
--   in double precision.

  function dznrm2 ( n : integer32;
                    xre : Standard_Floating_Vectors.Link_to_Vector;
                    xim : Standard_Floating_Vectors.Link_to_Vector;
                    ind,incx : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a vector given as a pair of vectors,
  --   as a vector of its real parts and a vector of its imaginary parts,
  --   starting at xre(ind), xim(ind), with n steps of increment incx.

  -- REQUIRED :
  --   xre'range = xim'range and xre'last >= ind + (n-1)*incx.

    ix : integer32;
    norm,scale,ssq,temp : double_float;

  begin
    if n < 1 or incx < 1 then
      norm := 0.0;
    else
      scale := 0.0;
      ssq := 1.0;
      ix := ind;
      while ix <= ind + (n-1)*incx loop
        if xre(ix) /= 0.0 then
          temp := abs(xre(ix));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2;
            scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        if xim(ix) /= 0.0 then
          temp := abs(xim(ix));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2;
            scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        ix := ix + incx;
      end loop;
      norm := scale*SQRT(ssq);
    end if;
    return norm;
  end dznrm2;

  function dznrm2 ( n : integer32;
                    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
                    row,col,incx : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Returns the Euclidean norm of a column of a matrix given
  --   as columns of real and imaginary parts.

  -- ON ENTRY :
  --   n        number of steps done in the matrix;
  --   rvv      real parts of the columns of the matrix;
  --   ivv      imaginary parts of the columns of the matrix;
  --   row      start row in the column col of the matrix;
  --   col      index to the column in the matrix;
  --   incx     increment in the matrix.

    ix : integer32;
    norm,scale,ssq,temp : double_float;
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if n < 1 or incx < 1 then
      norm := 0.0;
    else
      scale := 0.0; ssq := 1.0; ix := row;
      xre := rvv(col); xim := ivv(col);
      while ix <= row + (n-1)*incx loop
        if xre(ix) /= 0.0 then
          temp := abs(xre(ix));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2; scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        if xim(ix) /= 0.0 then
          temp := abs(xim(ix));
          if scale < temp then
            ssq := 1.0 + ssq*(scale/temp)**2; scale := temp;
          else
            ssq := ssq + (temp/scale)**2;
          end if;
        end if;
        ix := ix + incx;
      end loop;
      norm := scale*SQRT(ssq);
    end if;
    return norm;
  end dznrm2;

  procedure Test_Euclidean_Norm_of_Vector ( n : in integer32 ) is

  -- DESCRIPTION :
  --    Generates a random n-dimensional vector and compares
  --    the output of the dznrm2 on the complex vector with
  --    the output on the splitted vector.

    x : constant Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    xre,xim : Standard_Floating_Vectors.Link_to_Vector;
    nrm1,nrm2 : double_float;

  begin
    Standard_Vector_Splitters.Split_Complex(x,xre,xim);
    nrm1 := Standard_Complex_BLAS_Helpers.dznrm2(n,x,1,1);
    nrm2 := dznrm2(n,xre,xim,1,1);
    put("nrm1 : "); put(nrm1); new_line;
    put("nrm2 : "); put(nrm2); new_line;
  end Test_Euclidean_Norm_of_Vector;

  procedure Test_Euclidean_Norm_of_Column ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p complex matrix and compares the output of
  --   dznrm2 on the complex matrix with the output of the splitted matrix.

    A : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;
    nrm1,nrm2 : double_float;

  begin
    rvv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    ivv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(A,rvv,ivv);
    nrm1 := Standard_Complex_BLAS_Helpers.dznrm2(n,A,1,1,1);
    nrm2 := dznrm2(n,rvv,ivv,1,1,1);
    put("nrm1 : "); put(nrm1); new_line;
    put("nrm2 : "); put(nrm2); new_line;
  end Test_Euclidean_Norm_of_Column;

  procedure Test ( n,p : in integer32;
                   A : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given in A is an n-by-p matrix, computes its SVD.

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

  procedure Test_SVD ( n,p : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random n-by-p matrix and test the SVD.

    A : constant Standard_Complex_Matrices.Matrix(1..n,1..p)
      := Standard_Random_Matrices.Random_Matrix(1,n,1,p);
    rvv,ivv : Standard_Floating_VecVecs.Link_to_VecVec;

  begin
    put("A random "); put(n,1); put("-by-"); put(p,1);
    put_line(" complex matrix :"); put(A,2);
    rvv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    ivv := Standard_Vector_Splitters.Allocate(p,n,1,1);
    Standard_Matrix_Splitters.Complex_Parts(A,rvv,ivv);
    put_line("The real parts of the columns :"); put(rvv,2);
    put_line("The imaginary parts of the columns :"); put(ivv,2);
    Test(n,p,A);
  end Test_SVD;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the number of rows and columns
  --   and then generates a random complex matrix
  --   with that number of rows and columns.

    ans : character;

  begin
    new_line;
    put_line("MENU to test inlined SVD operations :");
    put_line("  1. Euclidean norm of a random vector; ");
    put_line("  2. Euclidean norm of a column of a random matrix; ");
    put_line("  3. SVD of a given matrix.");
    put("Type 1, 2, or 3 to select a test : ");
    Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' =>
        declare
          dim : integer32 := 0;
        begin
          put("Give the dimension : "); get(dim);
          Test_Euclidean_Norm_of_Vector(dim);
        end;
      when '2' | '3' =>
        declare
          nrows,ncols : integer32 := 0;
        begin
          put("Give the number of rows : "); get(nrows);
          put("Give the number of columns : "); get(ncols);
          if ans = '2'
           then Test_Euclidean_Norm_of_Column(nrows,ncols);
           else Test_SVD(nrows,ncols);
          end if;
        end;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_perfdsvd;
