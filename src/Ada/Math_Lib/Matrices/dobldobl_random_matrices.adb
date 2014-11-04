with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with Standard_Integer_Vectors;
with Double_Double_Vectors;
with DoblDobl_Complex_Vectors;
with Double_Double_QR_Least_Squares;     use Double_Double_QR_Least_Squares;
with DoblDobl_Complex_QR_Least_Squares;  use DoblDobl_Complex_QR_Least_Squares;

package body DoblDobl_Random_Matrices is

  function Random_Matrix ( n,m : natural32 )
                         return Double_Double_Matrices.Matrix is

    res : Double_Double_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DoblDobl_Random_Numbers.Random;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : Double_Double_Matrices.Matrix )
                         return Double_Double_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : Double_Double_Matrices.Matrix(1..n,1..m);
    wrk : Double_Double_Matrices.Matrix(1..n,1..m);
    bas : Double_Double_Matrices.Matrix(1..n,1..n);
    zero : constant double_double := create(0.0);
    qra : Double_Double_Vectors.Vector(1..n) := (1..n => zero);
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
                ( n,m : natural32 ) return Double_Double_Matrices.Matrix is

    res : constant Double_Double_Matrices.Matrix
                     (1..integer32(n),1..integer32(m))
        := Orthogonalize(Random_Matrix(n,m));

  begin
    return res;
  end Random_Orthogonal_Matrix;

  function Random_Matrix ( n,m : natural32 )
                         return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DoblDobl_Random_Numbers.Random1;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Orthogonalize ( mat : DoblDobl_Complex_Matrices.Matrix )
                         return DoblDobl_Complex_Matrices.Matrix is

    n : constant integer32 := mat'length(1);
    m : constant integer32 := mat'length(2);
    res : DoblDobl_Complex_Matrices.Matrix(1..n,1..m);
    wrk : DoblDobl_Complex_Matrices.Matrix(1..n,1..m);
    bas : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    zero : constant double_double := create(0.0);
    qra : DoblDobl_Complex_Vectors.Vector(1..n) := (1..n => Create(zero));
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
             ( n,m : natural32 ) return DoblDobl_Complex_Matrices.Matrix is

    res : constant DoblDobl_Complex_Matrices.Matrix
                     (1..integer32(n),1..integer32(m))
        := Orthogonalize(Random_Matrix(n,m));

  begin
    return res;
  end Random_Orthogonal_Matrix;

end DoblDobl_Random_Matrices;
