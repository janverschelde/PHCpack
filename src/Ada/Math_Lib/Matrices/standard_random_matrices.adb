with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Floating_QR_Least_Squares; use Standard_Floating_QR_Least_Squares;
with Standard_Complex_QR_Least_Squares;  use Standard_Complex_QR_Least_Squares;

package body Standard_Random_Matrices is

  function Random_Matrix
             ( dim : natural32; prb : double_float := 0.5 )
             return Boolean_Matrices.Matrix is

    idm : constant integer32 := integer32(dim);
    res : Boolean_Matrices.Matrix(1..idm,1..idm);

  begin
    for i in 1..idm loop
      for j in 1..idm loop
        res(i,j) := Standard_Random_Numbers.Random(prb);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix
             ( rows,cols : natural32; prb : double_float := 0.5 )
             return Boolean_Matrices.Matrix is

    nbrows : constant integer32 := integer32(rows);
    nbcols : constant integer32 := integer32(cols);
    res : Boolean_Matrices.Matrix(1..nbrows,1..nbcols);

  begin
    for i in 1..nbrows loop
      for j in 1..nbcols loop
        res(i,j) := Standard_Random_Numbers.Random(prb);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m : natural32; low,upp : integer32 )
                         return Standard_Integer_Matrices.Matrix is

    res : Standard_Integer_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(low,upp);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m : natural32; low,upp : integer64 )
                         return Standard_Integer64_Matrices.Matrix is

    res : Standard_Integer64_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(low,upp);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m : natural32 )
                         return Standard_Floating_Matrices.Matrix is

    res : Standard_Floating_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : Standard_Floating_Matrices.Matrix )
                         return Standard_Floating_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : Standard_Floating_Matrices.Matrix(1..n,1..m);
    wrk : Standard_Floating_Matrices.Matrix(1..n,1..m);
    bas : Standard_Floating_Matrices.Matrix(1..n,1..n);
    qra : Standard_Floating_Vectors.Vector(1..n) := (1..n => 0.0);
    pvt : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0); 

  begin
    wrk := mat;
    QRD(wrk,qra,pvt,false);
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        bas(i,j) := wrk(i,j);
      end loop;
      for j in m+1..n loop
        bas(i,j) := 0.0;
      end loop;
    end loop;
    Basis(bas,mat); 
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := bas(i,j);
      end loop;
    end loop;
    return res;
  end Orthogonalize;

  function Random_Orthogonal_Matrix
               ( n,m : natural32 ) return Standard_Floating_Matrices.Matrix is

    res : constant Standard_Floating_Matrices.Matrix
                     (1..integer32(n),1..integer32(m))
        := Orthogonalize(Random_Matrix(n,m));

  begin
    return res;
  end Random_Orthogonal_Matrix;

  function Random_Matrix ( n,m : natural32 )
                         return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random1;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( low_row,upp_row,low_col,upp_col : integer32 )
                         return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns a random matrix with row range low_row..upp_row,
  --   and column range low_col..upp_col.

    use Standard_Complex_Matrices;
    res : Matrix(low_row..upp_row,low_col..upp_col);

  begin
    for i in low_row..upp_row loop
      for j in low_col..upp_col loop
        res(i,j) := Random1;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : Standard_Complex_Matrices.Matrix )
                         return Standard_Complex_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : Standard_Complex_Matrices.Matrix(1..n,1..m);
    wrk : Standard_Complex_Matrices.Matrix(1..n,1..m);
    bas : Standard_Complex_Matrices.Matrix(1..n,1..n);
    qra : Standard_Complex_Vectors.Vector(1..n) := (1..n => Create(0.0));
    pvt : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0); 

  begin
    wrk := mat;
    QRD(wrk,qra,pvt,false);
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        bas(i,j) := wrk(i,j);
      end loop;
      for j in m+1..n loop
        bas(i,j) := Create(0.0);
      end loop;
    end loop;
    Basis(bas,mat); 
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := bas(i,j);
      end loop;
    end loop;
    return res;
  end Orthogonalize;

  function Random_Orthogonal_Matrix
                ( n,m : natural32 ) return Standard_Complex_Matrices.Matrix is

    res : constant Standard_Complex_Matrices.Matrix
                     (1..integer32(n),1..integer32(m))
        := Orthogonalize(Random_Matrix(n,m));

  begin
    return res;
  end Random_Orthogonal_Matrix;

end Standard_Random_Matrices;
