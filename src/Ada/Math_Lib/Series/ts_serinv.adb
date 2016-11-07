with text_io;                           use text_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Singular_Values;
with Standard_Dense_Matrix_Series;

procedure ts_serinv is

-- DESCRIPTION :
--   Development of the inverse of a matrix series,
--   even if the rank of the leading matrix is not full.

  procedure Read ( s : in out Standard_Dense_Matrix_Series.Matrix;
                   n : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the n-by-n matrices for the coefficients in s,
  --   up to the degree of s.

    cff : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for i in 0..s.deg loop
      put("Reading the coefficient of degree ");
      put(i,1); put_line(" ...");
      put("-> give an "); put(n,1); put("-by-"); put(n,1);
      put_line(" matrix : "); get(cff);
      s.cff(i) := new Standard_Complex_Matrices.Matrix'(cff);
    end loop;
  end Read;

  procedure Write ( s : in Standard_Dense_Matrix_Series.Matrix ) is
  begin
    put("A series of degree "); put(s.deg,1); put_line(" :");
    for i in 0..s.deg loop
      put("-> coefficient of degree "); put(i,1); put_line(" :");
      put(s.cff(i).all,3);
    end loop;
  end Write;

  function Stack ( s : Standard_Dense_Matrix_Series.Matrix )
                 return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Stacks the coefficients in the series s into one big matrix
  --   to determine the column rank.
  --   The dimension of the matrix are as follows:
  --   number of rows : dimension of the matrices times the degree + 2,
  --   number of columns : dimension of the matrices times the degree + 1.

  -- REQUIRED : s.cff(0) is not null.

    dim : constant integer32 := s.cff(0)'last(1);
    rws : constant integer32 := (2*s.deg+1)*dim;
    cls : constant integer32 := (s.deg+1)*dim;
    res : Standard_Complex_Matrices.Matrix(1..rws,1..cls);
    mat : Standard_Complex_Matrices.Link_to_Matrix;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := create(0.0);
      end loop;
    end loop;
    for k in 0..s.deg loop
      mat := s.cff(k);
      for ll in 0..s.deg loop
        for i in 1..dim loop
          for j in 1..dim loop
            res((ll+k)*dim+i,ll*dim+j) := mat(i,j);
          end loop;
        end loop;
      end loop;
    end loop;
    return res;
  end Stack;

  function Rank ( A : Standard_Complex_Matrices.Matrix;
                  verbose : boolean := true ) 
                return integer32 is

  -- DESCRIPTION :
  --   Returns the numerical rank of the matrix A,
  --   based on the computation of a singular value decomposition of A.

    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    res : integer32 := 0;
    wrk : Standard_Complex_Matrices.Matrix(1..n,1..p) := A;
    mm : constant integer32 := Standard_Complex_Singular_Values.Min0(n+1,p);
    s : Standard_Complex_Vectors.Vector(1..mm);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Standard_Complex_Matrices.Matrix(1..n,1..n);
    v : Standard_Complex_Matrices.Matrix(1..p,1..p);
    job : constant integer32 := 11;
    info : integer32;

  begin
    Standard_Complex_Singular_Values.SVD(wrk,n,p,s,e,u,v,job,info);
    put_line("The singular values :"); put_line(s);
    res := Standard_Complex_Singular_Values.Rank(s);
    return res;
  end Rank;

  procedure Test ( s : in Standard_Dense_Matrix_Series.Matrix ) is

  -- DESCRIPTION :
  --   Tests the rank of the stacked coefficient matrix.

    dim : constant integer32 := s.cff(0)'last(1);
    rws : constant integer32 := (2*s.deg+1)*dim;
    cls : constant integer32 := (s.deg+1)*dim;
    cff : constant Standard_Complex_Matrices.Matrix(1..rws,1..cls)
        := Stack(s);
    rnk : integer32;
  
  begin
    put_line("The stacked coefficient matrix : "); put(cff,2);
    rnk := Rank(cff);
    put("The rank of the stacked matrix : "); put(rnk,1); new_line;
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the degree and the dimension of the matrices.
  
    s : Standard_Dense_Matrix_Series.Matrix;
    n : integer32 := 0;

  begin
    s.deg := 0;
    new_line;
    put_line("Reading a matrix series ...");
    put("-> give the degree of the series : "); get(s.deg);
    put("-> give the dimension of the matrices : "); get(n);
    Read(s,n);
    Write(s);
    Test(s);
  end Main;

begin
  Main;
end ts_serinv;
