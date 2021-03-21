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
with Standard_Random_Vectors;
with Standard_Random_Matrices;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecVecs;
with Standard_Laurent_Series;

procedure ts_laurmat is

-- DESCRIPTION :
--   A matrix of Laurent series is a tuple of
--   1) a matrix of leading exponents, and
--   2) 3-dimensional vector of vector of vectors
--   with the coefficients of the power series.

  procedure Random_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

  -- DESCRIPTION :
  --   Given a fully allocated 3-dimensional v,
  --   fills it up with random complex numbers on the unit circle.

  begin
    for i in v'range loop
      declare
        vi : constant Standard_Complex_VecVecs.Link_to_VecVec := v(i);
      begin
        for j in vi'range loop
          declare
            vij : constant Standard_Complex_Vectors.Link_to_Vector := vi(j);
          begin
            for k in vij'range loop
              vij(k) := Standard_Random_Numbers.Random1;
            end loop;
          end;
        end loop;
      end;
    end loop;
  end Random_VecVecVec;

  procedure Random_Lower_VecVecVec
              ( v : in Standard_Complex_VecVecVecs.Link_to_VecVecVec ) is

  -- DESCRIPTION :
  --   Given a fully allocated 3-dimensional v,
  --   fills it up with random complex numbers on the unit circle,
  --   but only on the lower triangular part, below the diagonal.
  --   The diagonal elements are set to one.

  begin
    for i in v'range loop -- row i
      declare
        vi : constant Standard_Complex_VecVecs.Link_to_VecVec := v(i);
        vii : constant Standard_Complex_Vectors.Link_to_Vector := vi(i);
      begin
        for j in vi'first..(i-1) loop -- column j
          declare
            vij : constant Standard_Complex_Vectors.Link_to_Vector := vi(j);
          begin
            for k in vij'range loop
              vij(k) := Standard_Random_Numbers.Random1;
            end loop;
          end;
        end loop;
        vii(0) := Standard_Complex_Numbers.Create(1.0);
      end;
    end loop;
  end Random_Lower_VecVecVec;

  procedure Random_Series_Coefficients
              ( dim,deg : in integer32;
                cff : out Standard_Complex_VecVecs.Link_to_VecVec ) is

  -- DESCRIPTION :
  --   Returns in cff the coefficients of dim series of degree deg,
  --   randomly generated on the complex unit circle.

    res : Standard_Complex_VecVecs.VecVec(1..dim);

  begin
    for i in 1..dim loop
      declare
        val : constant Standard_Complex_Vectors.Vector(0..deg)
            := Standard_Random_Vectors.Random_Vector(0,deg);
      begin
        res(i) := new Standard_Complex_Vectors.Vector'(val);
      end;
    end loop;
    cff := new Standard_Complex_VecVecs.VecVec'(res);
  end Random_Series_Coefficients;

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
        cyi : Standard_Complex_Vectors.Link_to_Vector := cy(i);
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

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the dimensions and then generates a random matrix.

    deg,nrows,ncols,low,upp : integer32 := 0;
    ans : character;
    lower : boolean;

  begin
    new_line;
    put("Give the truncation degree : "); get(deg);
    put("Give the lower bound on the leading exponent : "); get(low);
    put("Give the upper bound on the leading exponent : "); get(upp);
    put("Give the number of rows : "); get(nrows);
    put("Give the number of columns : "); get(ncols);
    new_line;
    put("Lower triangular matrix ? (y/n) "); Ask_Yes_or_No(ans);
    lower := (ans = 'y');
    declare
      nbrows : constant natural32 := natural32(nrows);
      nbcols : constant natural32 := natural32(ncols);
      Alead : constant Standard_Integer_Matrices.Matrix(1..nrows,1..ncols)
            := Standard_Random_Matrices.Random_Matrix(nbrows,nbcols,low,upp);
      xlead : constant Standard_Integer_Vectors.Vector(1..ncols)
            := Standard_Random_Vectors.Random_Vector(1,ncols,low,upp);
      Acffs : Standard_Complex_VecVecVecs.Link_to_VecVecVec;
      xcffs : Standard_Complex_VecVecs.Link_to_VecVec;
      blead : Standard_Integer_Vectors.Vector(1..nrows);
      bcffs : Standard_Complex_VecVecs.Link_to_VecVec;
    begin
      put_line("The matrix of leading exponents :"); put(Alead,1);
      put("The vector of leading exponents :"); put(xlead,1); new_line;
      Standard_Complex_VecVecVecs.Allocate(Acffs,1,nrows,1,ncols,0,deg);
      if lower
       then Random_Lower_VecVecVec(Acffs);
       else Random_VecVecVec(Acffs);
      end if;
      put("A "); put(nrows,1); put("-by-"); put(ncols,1);
      put_line(" matrix of Laurent series : "); Write(Alead,Acffs);
      Random_Series_Coefficients(ncols,deg,xcffs);
      put("A "); put(ncols,1); put_line("-vector of Laurent series :");
      Write(xlead,xcffs,"x");
      Allocate_Series_Coefficients(nrows,deg,bcffs);
      Matrix_Vector_Product(deg,Alead,Acffs,xlead,xcffs,blead,bcffs);
      put_line("The product of the matrix with the vector :");
      Write(blead,bcffs,"b");
    end;
  end Main;

begin
  Main;
end ts_laurmat;
