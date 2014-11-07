with Standard_Random_Numbers;
with Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Random_Numbers;            use Multprec_Random_Numbers;
with Standard_Integer_Vectors;
with Multprec_Floating_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Floating_QR_Least_Squares; use Multprec_Floating_QR_Least_Squares;
with Multprec_Complex_QR_Least_Squares;  use Multprec_Complex_QR_Least_Squares;

package body Multprec_Random_Matrices is

  function Random_Matrix ( n,m : natural32; low,upp : integer32 )
                         return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..integer32(n),1..integer32(m));
    rnd : integer32;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        rnd := Standard_Random_Numbers.Random(low,upp);
        res(i,j) := Multprec_Integer_Numbers.Create(rnd);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Integer_Matrices.Matrix is

    res : Multprec_Integer_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Integer64_Matrices.Matrix is

    res : Multprec_Integer64_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Floating_Matrices.Matrix is

    res : Multprec_Floating_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : Multprec_Floating_Matrices.Matrix )
                         return Multprec_Floating_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : Multprec_Floating_Matrices.Matrix(1..n,1..m);
    wrk : Multprec_Floating_Matrices.Matrix(1..n,1..m);
    bas : Multprec_Floating_Matrices.Matrix(1..n,1..n);
    zero : Floating_Number := Create(0.0);
    qra : Multprec_Floating_Vectors.Vector(1..n);
    pvt : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0); 

  begin
    for i in 1..n loop
      Copy(zero,qra(i));
    end loop;
    Multprec_Floating_Matrices.Copy(mat,wrk);
    QRD(wrk,qra,pvt,false);
    Multprec_Floating_Vectors.Clear(qra);
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        Copy(wrk(i,j),bas(i,j));
      end loop;
      for j in m+1..n loop
        Copy(zero,bas(i,j));
      end loop;
    end loop;
    Basis(bas,mat); 
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := bas(i,j);
      end loop;
    end loop;
    Multprec_Floating_Numbers.Clear(zero);
    Multprec_Floating_Matrices.Clear(wrk);
    return res;
  end Orthogonalize;

  function Random_Orthogonal_Matrix
             ( n,m,sz : natural32 ) return Multprec_Floating_Matrices.Matrix is

    res : Multprec_Floating_Matrices.Matrix(1..integer32(n),1..integer32(m));
    rnd : Multprec_Floating_Matrices.Matrix(1..integer32(n),1..integer32(m))
        := Random_Matrix(n,m,sz);

  begin
    res := Orthogonalize(rnd);
    Multprec_Floating_Matrices.Clear(rnd);
    return res;
  end Random_Orthogonal_Matrix;

  function Random_Matrix ( n,m,sz : natural32 )
                         return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random(sz);
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : Multprec_Complex_Matrices.Matrix )
                         return Multprec_Complex_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    wrk : Multprec_Complex_Matrices.Matrix(1..n,1..m);
    bas : Multprec_Complex_Matrices.Matrix(1..n,1..n);
    mpflt_zero : Floating_Number := Create(0.0);
    cmplx_zero : Complex_Number := Create(mpflt_zero);
    qra : Multprec_Complex_Vectors.Vector(1..n);
    pvt : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0); 

  begin
    for i in 1..n loop
      Copy(cmplx_zero,qra(i));
    end loop;
    Multprec_Complex_Matrices.Copy(mat,wrk);
    QRD(wrk,qra,pvt,false);
    Multprec_Complex_Vectors.Clear(qra);
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        Copy(wrk(i,j),bas(i,j));
      end loop;
      for j in m+1..n loop
        Copy(cmplx_zero,bas(i,j));
      end loop;
    end loop;
    Basis(bas,mat); 
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := bas(i,j);
      end loop;
    end loop;
    Multprec_Floating_Numbers.Clear(mpflt_zero);
    Multprec_Complex_Numbers.Clear(cmplx_zero);
    Multprec_Complex_Matrices.Clear(wrk);
    return res;
  end Orthogonalize;

  function Random_Orthogonal_Matrix
             ( n,m,sz : natural32 ) return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));
    rnd : Multprec_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m))
        := Random_Matrix(n,m,sz);

  begin
    res := Orthogonalize(rnd);
    Multprec_Complex_Matrices.Clear(rnd);
    return res;
  end Random_Orthogonal_Matrix;

end Multprec_Random_Matrices;
