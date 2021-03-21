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

  procedure Matrix_Vector_Product
              ( d : in integer32;
                eA : in Standard_Integer_Matrices.Matrix;
                cA : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                ex : in Standard_Integer_Vectors.Vector;
                cx : in Standard_Complex_VecVecs.Link_to_VecVec;
                ey : out Standard_Integer_Vectors.Vector;
                cy : out Standard_Complex_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Returns the product of a matrix of Laurent series
  --   with a vector of Laurent series.

  -- REQUIRED :
  --   eA'range(1) = ey'range and eA'range(2) = ex'range,
  --   cA'range(1) = cy'range and cA'range(2) = cx'range.

  -- ON ENTRY :
  --   d        only coefficients in the range 0 to d are considered;
  --   eA       leading exponents of the series in the matrix A;
  --   cA       coefficients of the series in the matrix A;
  --   ex       leading exponents of the series in the vector x;
  --   cx       coefficients of the series in the vector x;
  --   cy       space allocated for the coefficients of the product.

  -- ON RETURN :
  --   ey       leading exponents of the series of the product;
  --   cy       leading coefficients of the series of the product.

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..d);

  begin
    for i in eA'range(1) loop
      declare
        cAi : constant Standard_Complex_VecVecs.Link_to_VecVec := cA(i);
        cyi : constant Standard_Complex_Vectors.Link_to_Vector := cy(i);
      begin -- initialize with first product instead of with zero
        Standard_Laurent_Series.Multiply
          (d,eA(i,eA'first(2)),ex(ex'first),cAi(cAi'first).all,
           cx(cx'first).all,ey(i),cyi.all);
        for j in eA'first(2)+1..eA'last(2) loop
          Standard_Laurent_Series.Multiply
            (d,eA(i,j),ex(j),cAi(j).all,cx(j).all,ze,zc);
          Standard_Laurent_Series.Add(d,ey(i),ze,cyi.all,zc,ewrk,cwrk);
          ey(i) := ewrk;
          for k in 0..d loop
            cyi(k) := cwrk(k);
          end loop;
        end loop;
      end;
    end loop;
  end Matrix_Vector_Product;

  procedure Forward_Substitution
              ( d : in integer32;
                eL : in Standard_Integer_Matrices.Matrix;
                cL : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                eb : in Standard_Integer_Vectors.Vector;
                cb : in Standard_Complex_VecVecs.Link_to_VecVec;
                ex : out Standard_Integer_Vectors.Vector;
                cx : out Standard_Complex_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Applies forward substitution to solve a lower triangular system
  --   with ones on the diagonal.

  -- REQUIRED :
  --   The matrix is square and all ranges are compatible.

  -- ON ENTRY :
  --   d        only coefficients in the range 0 to d are considered;
  --   eL       leading exponents in the lower triangular matrix L;
  --   cL       coefficients of the series in the matrix L;
  --   eb       leading exponents of the right hand side vector b;
  --   cb       coefficients of the series in the vector b;
  --   cx       space allocated for the coefficients of the solution.

  -- ON RETURN :
  --   ex       leading exponents of the series of the solution;
  --   cx       leading coefficients of the series of the solution.

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..d);

  begin
    for i in eb'range loop
      ex(i) := eb(i);
      declare
        cbi : constant Standard_Complex_Vectors.Link_to_Vector := cb(i);
        cLi : constant Standard_Complex_VecVecs.Link_to_VecVec := cL(i);
        cxi : constant Standard_Complex_Vectors.Link_to_Vector := cx(i);
      begin
        for k in 0..d loop
          cxi(k) := cbi(k);
        end loop;
        for j in ex'first..(i-1) loop
          Standard_Laurent_Series.Multiply
            (d,eL(i,j),ex(j),cLi(j).all,cx(j).all,ze,zc);
         -- put("L("); put(i,1); put(","); put(j,1);
         -- put(")*x("); put(j,1); put_line(") :");
         -- Standard_Laurent_Series.Write(ze,zc);
          Standard_Laurent_Series.Subtract(d,ex(i),ze,cxi.all,zc,ewrk,cwrk);
          ex(i) := ewrk;
          for k in 0..d loop
            cxi(k) := cwrk(k);
          end loop;
         -- put("x("); put(i,1); put_line(") after the update :");
         -- Standard_Laurent_Series.Write(ex(i),cxi.all);
        end loop;
      end;
    end loop;
  end Forward_Substitution;

  procedure Backward_Substitution
              ( d : in integer32;
                eU : in Standard_Integer_Matrices.Matrix;
                cU : in Standard_Complex_VecVecVecs.Link_to_VecVecVec;
                eb : in Standard_Integer_Vectors.Vector;
                cb : in Standard_Complex_VecVecs.Link_to_VecVec;
                ex : out Standard_Integer_Vectors.Vector;
                cx : out Standard_Complex_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Applies forward substitution to solve an upper triangular system
  --   with general, nonzero elements on the diagonal.

  -- REQUIRED :
  --   The matrix is square and all ranges are compatible.

  -- ON ENTRY :
  --   d        only coefficients in the range 0 to d are considered;
  --   eU       leading exponents in the upper triangular matrix U;
  --   cU       coefficients of the series in the matrix U;
  --   eb       leading exponents of the right hand side vector b;
  --   cb       coefficients of the series in the vector b;
  --   cx       space allocated for the coefficients of the solution.

  -- ON RETURN :
  --   ex       leading exponents of the series of the solution;
  --   cx       leading coefficients of the series of the solution.

    ze,ewrk : integer32;
    zc,cwrk : Standard_Complex_Vectors.Vector(0..d);
    cbi,cxi : Standard_Complex_Vectors.Link_to_Vector;
    cUi : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
    cbi := cb(cb'last); cxi := cx(cx'last); cUi := cU(cU'last);
    Standard_Laurent_Series.Divide
      (d,eb(eb'last),eU(eU'last(1),eU'last(2)),cbi.all,cUi(cUi'last).all,
       ex(ex'last),cxi.all,cwrk);
    for i in reverse eb'first..eb'last-1 loop
      ex(i) := eb(i);
      cbi := cb(i); cxi := cx(i); cUi := cU(i);
      for k in 0..d loop
        cxi(k) := cbi(k);
      end loop;
      for j in (i+1)..ex'last loop
        Standard_Laurent_Series.Multiply
          (d,eU(i,j),ex(j),cUi(j).all,cx(j).all,ze,zc);
       -- put("U("); put(i,1); put(","); put(j,1);
       -- put(")*x("); put(j,1); put_line(") :");
       -- Standard_Laurent_Series.Write(ze,zc);
        Standard_Laurent_Series.Subtract(d,ex(i),ze,cxi.all,zc,ewrk,cwrk);
        ex(i) := ewrk;
        for k in 0..d loop
          cxi(k) := cwrk(k);
        end loop;
       -- put("x("); put(i,1); put_line(") after the update :");
       -- Standard_Laurent_Series.Write(ex(i),cxi.all);
      end loop;
      Standard_Laurent_Series.Divide
        (d,ex(i),eU(i,i),cxi.all,cUi(i).all,ze,zc,cwrk);
      ex(i) := ze;
      for k in 0..d loop
        cxi(k) := zc(k);
      end loop;
    end loop;
  end Backward_Substitution;

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
