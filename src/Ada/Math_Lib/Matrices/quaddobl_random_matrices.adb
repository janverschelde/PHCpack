with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Integer_Vectors;
with Quad_Double_Vectors;
with QuadDobl_Complex_Vectors;
with Quad_Double_QR_Least_Squares;       use Quad_Double_QR_Least_Squares;
with QuadDobl_Complex_QR_Least_Squares;  use QuadDobl_Complex_QR_Least_Squares;

package body QuadDobl_Random_Matrices is

  function Random_Matrix ( n,m : natural32 )
                         return Quad_Double_Matrices.Matrix is

    res : Quad_Double_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Random_Numbers.Random;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : Quad_Double_Matrices.Matrix )
                         return Quad_Double_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : Quad_Double_Matrices.Matrix(1..n,1..m);
    wrk : Quad_Double_Matrices.Matrix(1..n,1..m);
    bas : Quad_Double_Matrices.Matrix(1..n,1..n);
    zero : constant quad_double := create(0.0);
    qra : Quad_Double_Vectors.Vector(1..n) := (1..n => zero);
    pvt : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0); 

  begin
    wrk := mat;
    QRD(wrk,qra,pvt,false);
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        bas(i,j) := wrk(i,j);
      end loop;
      for j in m+1..n loop
        bas(i,j) := zero;
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
                ( n,m : natural32 ) return Quad_Double_Matrices.Matrix is

    res : constant Quad_Double_Matrices.Matrix
                     (1..integer32(n),1..integer32(m))
        := Orthogonalize(Random_Matrix(n,m));

  begin
    return res;
  end Random_Orthogonal_Matrix;

  function Random_Matrix ( n,m : natural32 )
                         return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Random_Numbers.Random1;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : QuadDobl_Complex_Matrices.Matrix )
                         return QuadDobl_Complex_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);
    wrk : QuadDobl_Complex_Matrices.Matrix(1..n,1..m);
    bas : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant quad_double := create(0.0);
    qra : QuadDobl_Complex_Vectors.Vector(1..n) := (1..n => Create(zero));
    pvt : Standard_Integer_Vectors.Vector(1..n) := (1..n => 0); 

  begin
    wrk := mat;
    QRD(wrk,qra,pvt,false);
    for i in wrk'range(1) loop
      for j in wrk'range(2) loop
        bas(i,j) := wrk(i,j);
      end loop;
      for j in m+1..n loop
        bas(i,j) := Create(zero);
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
             ( n,m : natural32 ) return QuadDobl_Complex_Matrices.Matrix is

    res : constant QuadDobl_Complex_Matrices.Matrix
                     (1..integer32(n),1..integer32(m))
        := Orthogonalize(Random_Matrix(n,m));

  begin
    return res;
  end Random_Orthogonal_Matrix;

end QuadDobl_Random_Matrices;
