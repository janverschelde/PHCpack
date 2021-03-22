with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Random_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecVecs;
with Standard_Laurent_Series;
with Random_Laurent_Series;             use Random_Laurent_Series;
with Standard_Linear_Laurent_Solvers;   use Standard_Linear_Laurent_Solvers;

procedure ts_laurmat is

-- DESCRIPTION :
--   A matrix of Laurent series is a tuple of
--   1) a matrix of leading exponents, and
--   2) 3-dimensional vector of vector of vectors
--   with the coefficients of the power series.

  procedure Allocate_Series_Coefficients
              ( dim,deg : in integer32;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Returns in cff the coefficients of dim series of degree deg,
  --   all equal to zero.

    res : Standard_Complex_VecVecs.VecVec(1..dim);
    zero : constant Standard_Complex_Numbers.Complex_Number
         := Standard_Complex_Numbers.Create(0.0);

  begin
    for i in 1..dim loop
      declare
        val : constant Standard_Complex_Vectors.Vector(0..deg)
            := (0..deg => zero);
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(val);
      end;
    end loop;
    cff := new Standard_Complex_VecVecs.VecVec'(res);
  end Allocate_Series_Coefficients;

  procedure Write ( e : in Standard_Integer_Matrices.Matrix;
                    c : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                    s : in string := "A" ) is

  -- DESCRIPTION :
  --   Writes the matrix of Laurent series, defined by
  --   the leading exponents in e and coefficients in c.

  -- REQUIRED :
  --   e'range(1) = c'range(1) and e'range(2) = c'range(2).

  begin
    for i in e'range(1) loop
      for j in e'range(2) loop
        put(s & "("); put(i,1); put(","); put(j,1); put_line(") :");
        Standard_Laurent_Series.Write(e(i,j),c(i)(j).all);
      end loop;
    end loop;
  end Write;

  procedure Write ( e : in Standard_Integer_Vectors.Vector;
                    c : in Standard_Complex_VecVecs.Link_to_VecVec;
                    s : in string := "v" ) is

  -- DESCRIPTION :
  --   Writes the vector of Laurent series, defined by
  --   the leading exponents in e and coefficients in c.

  -- REQUIRED : e'range = c'range.

  begin
    for i in e'range loop
      put(s & "("); put(i,1); put_line(") :");
      Standard_Laurent_Series.Write(e(i),c(i).all);
    end loop;
  end Write;

  procedure Test ( nrows,ncols,deg,low,upp : in integer32;
                   lower,upper : in boolean ) is

  -- DESCRIPTION :
  --   Generates a random linear system and runs a test.

  -- ON ENTRY :
  --   nrows    number of rows of the test matrix;
  --   ncols    number of columns of the test matrix;
  --   deg      degree of the series in the matrix;
  --   low      lower bound for the leading exponents;
  --   upp      upper bound for the leading exponents;
  --   lower    true if the test matrix is lower triangular;
  --   upper    true if the test matrix is upper triangular.

    nbrows : constant natural32 := natural32(nrows);
    nbcols : constant natural32 := natural32(ncols);
    Alead : Standard_Integer_Matrices.Matrix(1..nrows,1..ncols);
    xlead : Standard_Integer_Vectors.Vector(1..ncols);
    Acffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
    xcffs,ycffs : Standard_Complex_VecVecs.Link_to_VecVec;
    ylead : Standard_Integer_Vectors.Vector(1..ncols);
    blead : Standard_Integer_Vectors.Vector(1..nrows);
    bcffs : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
    Standard_Complex_VecVecVecs.Allocate(Acffs,1,nrows,1,ncols,0,deg);
    if not lower and not upper then
      Random_Matrix(nbrows,nbcols,low,upp,Alead,Acffs);
    elsif lower then
      Random_Lower_Matrix(nbrows,nbcols,low,upp,Alead,Acffs);
    else -- upper must be true
      Random_Upper_Matrix(nbrows,nbcols,low,upp,Alead,Acffs);
    end if;
    put_line("The matrix of leading exponents :"); put(Alead,1);
    put("A "); put(nrows,1); put("-by-"); put(ncols,1);
    put_line(" matrix of Laurent series : "); Write(Alead,Acffs);
    Random_Vector(ncols,deg,low,upp,xlead,xcffs);
    put("The vector of leading exponents :"); put(xlead,1); new_line;
    put("A "); put(ncols,1); put_line("-vector of Laurent series :");
    Write(xlead,xcffs,"x");
    Allocate_Series_Coefficients(nrows,deg,bcffs);
    Matrix_Vector_Product(deg,Alead,Acffs,xlead,xcffs,blead,bcffs);
    put_line("The product of the matrix with the vector :");
    Write(blead,bcffs,"b");
    Allocate_Series_Coefficients(ncols,deg,ycffs);
    if lower then
      Forward_Substitution(deg,Alead,Acffs,blead,bcffs,ylead,ycffs);
      put_line("The computed solution :");
      Write(ylead,ycffs,"y");
    elsif upper then
      Backward_Substitution(deg,Alead,Acffs,blead,bcffs,ylead,ycffs);
      put_line("The computed solution :");
      Write(ylead,ycffs,"y");
    end if;
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the dimensions and then generates a random matrix.

    deg,nrows,ncols,low,upp,seed : integer32 := 0;
    ans : character;
    lower,upper : boolean;

  begin
    new_line;
    put("Give the truncation degree : "); get(deg);
    put("Give the lower bound on the leading exponent : "); get(low);
    put("Give the upper bound on the leading exponent : "); get(upp);
    put("Give the number of rows : "); get(nrows);
    put("Give the number of columns : "); get(ncols);
    new_line;
    if nrows /= ncols then
      lower := false;
    else
      put("Lower triangular matrix ? (y/n) "); Ask_Yes_or_No(ans);
      lower := (ans = 'y');
      if lower then
        upper := false;
      else
        put("Upper triangular matrix ? (y/n) "); Ask_Yes_or_No(ans);
        upper := (ans = 'y');
      end if;
    end if;
    new_line;
    put("Fixed seed ? (y/n) "); Ask_Yes_or_No(ans);
    if ans /= 'y' then
      seed := Standard_Random_Numbers.Get_Seed;
    else
      put("Give the seed : "); get(seed);
      Standard_Random_Numbers.Set_Seed(natural32(seed));
    end if;
    Test(nrows,ncols,deg,low,upp,lower,upper);
    new_line;
    put("The seed used : "); put(seed,1); new_line;
  end Main;

begin
  Main;
end ts_laurmat;
