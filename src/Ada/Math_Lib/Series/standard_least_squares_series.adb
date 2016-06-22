with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Dense_Series;               use Standard_Dense_Series;
with Standard_Dense_Series_Norms;         use Standard_Dense_Series_Norms;

package body Standard_Least_Squares_Series is

-- AUXILIARIES :

  function min0 ( a,b : integer32 ) return integer32 is

  -- DESCRIPTION : returns the minimum of a and b.

  begin
    if a <= b
     then return a;
     else return b;
    end if;
  end min0;

  function dmax1 ( a,b : double_float ) return double_float is

  -- DESCRIPTION : returns the maximum of a and b.

  begin
    if a >= b
     then return a;
     else return b;
    end if;
  end dmax1;

  function cdabs ( s : Series ) return double_float is

  -- DESCRIPTION :
  --   Computes the 2-norm of the series s.

    res : constant double_float := Two_Norm(s);

  begin
    return res;
  end cdabs;

  function csign ( a,b : Series ) return Series is

  -- DESCRIPTION : translated from
  --       csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2)) 

    use Standard_Complex_Numbers;

    res : Series;
    fac : constant double_float := cdabs(a)/cdabs(b);
    cff : constant Complex_Number := Create(fac);

  begin
    res := cff*b;
    return res;
  end csign;

  procedure zswap ( a : in out Standard_Dense_Series_Matrices.Matrix;
                    k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   Swaps the columns k1 and k2 in the matrix a.

    tmp : Series;

  begin
    for i in a'range(1) loop
      tmp := a(i,k1); a(i,k1) := a(i,k2); a(i,k2) := tmp;
    end loop;
  end zswap;

  procedure zcopy ( n,start : in integer32;
                    x : in Standard_Dense_Series_Vectors.Vector;
                    y : out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Copies the n entries of x to y, beginning at start.

  begin
    for i in start..start+n-1 loop
      y(i) := x(i);
    end loop;
  end zcopy;
 
  function znrm2 ( a : Standard_Dense_Series_Matrices.Matrix;
                   row,col : integer32 ) return double_float is

  -- DESCRIPTION :
  --   Computes the 2-norm of the vector in the column col of the matrix,
  --   starting at the given row.

    sum : Series := Create(0.0);

  begin
    for i in row..a'last(1) loop
      sum := sum + Conjugate(a(i,col))*a(i,col);
    end loop;
    return Max_Norm(sum); -- SQRT(REAL_PART(sum));
  end znrm2;

  function zdot ( a : Standard_Dense_Series_Matrices.Matrix;
                  row,k1,k2 : integer32 ) return Series is

  -- DESCRIPTION :
  --   Returns the inner product of the vectors in the columns k1 and k2,
  --   starting at the given row.

    res : Series := Create(0.0);

  begin
    for i in row..a'last(1) loop
      res := res + Conjugate(a(i,k1))*a(i,k2);
    end loop;
    return res;
  end zdot;

  function zdotc ( row : integer32;
                   x : Standard_Dense_Series_Matrices.Matrix;
                   y : Standard_Dense_Series_Vectors.Vector )
                 return Series is

  -- DESCRIPTION :
  --   Dot product of two vectors : x(row..x'last(1),row)*y(row..y'last).

    res : Series := Create(0.0);

  begin
    for i in row..y'last loop
      res := res + Conjugate(x(i,row))*y(i);
    end loop;
    return res;
  end zdotc;

  procedure zaxpy ( a : in out Standard_Dense_Series_Matrices.Matrix; 
                    f : in Series; row,k1,k2 : in integer32 ) is

  -- DESCRIPTION :
  --   The column k2 is added with f times the column k1, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,k2) := a(i,k2) + f*a(i,k1);
    end loop;
  end zaxpy;

  procedure zscal ( a : in out Standard_Dense_Series_Matrices.Matrix;
                    f : in Series; row,col : in integer32 ) is

  -- DESCRIPTION :
  --   Multiplies the column col of the matrix with f, starting at row.

  begin
    for i in row..a'last(1) loop
      a(i,col) := f*a(i,col);
    end loop;
  end zscal;

  procedure zaxpy ( n,row,col : in integer32; f : in Series;
                    x : in Standard_Dense_Series_Matrices.Matrix;
                    y : in out Standard_Dense_Series_Vectors.Vector ) is

  -- DESCRIPTION :
  --   y(i) := y(i) + f*x(i,col) for i in row..row+n-1.

  begin
    for i in row..row+n-1 loop
      y(i) := y(i) + f*x(i,col);
    end loop;
  end zaxpy;

-- TARGET PROCEDURES :

  procedure QRD ( x : in out Standard_Dense_Series_Matrices.Matrix;
                  qraux : in out Standard_Dense_Series_Vectors.Vector;
                  jpvt : in out Standard_Integer_Vectors.Vector;
                  piv : in boolean ) is
  begin
    null;
  end QRD;

end Standard_Least_Squares_Series;
