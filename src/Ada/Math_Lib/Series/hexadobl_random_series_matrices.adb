with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Hexa_Double_Numbers;               use Hexa_Double_Numbers;
with HexaDobl_Complex_Numbers;
with Standard_Random_Numbers;
with HexaDobl_Random_Numbers;
with HexaDobl_Complex_Matrices;
with HexaDobl_Complex_Random_Series;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_Random_Series_Vectors;

package body HexaDobl_Random_Series_Matrices is

  function Random_Series_VecVec
             ( vvfirst,vvlast,first,last,degree : integer32 ) return VecVec is

    res : VecVec(vvfirst..vvlast);

  begin
    for i in res'range loop
      declare
        v : constant HexaDobl_Complex_Series_Vectors.Vector(first..last)
          := HexaDobl_Random_Series_Vectors.Random_Series_Vector
               (first,last,degree);
      begin
        res(i) := new HexaDobl_Complex_Series_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Random_Series_VecVec;

  function Random_Series_Matrix
             ( rowfirst,rowlast,columnfirst,columnlast,degree : integer32 )
             return HexaDobl_Complex_Series_Matrices.Matrix is

    use HexaDobl_Complex_Series_Matrices;

    res : Matrix(rowfirst..rowlast,columnfirst..columnlast);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := HexaDobl_Complex_Random_Series.Random_Series(degree);
      end loop;
    end loop;
    return res;
  end Random_Series_Matrix;

  function Random_Matrix_Series
             ( deg,dim,lower,upper : integer32 )
             return HexaDobl_Complex_Matrix_Series.Matrix is

    res : HexaDobl_Complex_Matrix_Series.Matrix(deg);

  begin
    for d in 0..res.deg loop
      declare
        mat : HexaDobl_Complex_Matrices.Matrix(1..dim,1..dim);
        rnd : integer32;
        nbr : hexa_double;
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            rnd := Standard_Random_Numbers.Random(lower,upper);
            nbr := Create(double_float(rnd));
            mat(i,j) := HexaDobl_Complex_Numbers.Create(nbr);
          end loop;
        end loop;
        res.cff(d) := new HexaDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Random_Matrix_Series;

  function Random_Matrix_Series
             ( deg,dim : integer32 )
             return HexaDobl_Complex_Matrix_Series.Matrix is

    res : HexaDobl_Complex_Matrix_Series.Matrix(deg);

  begin
    for d in 0..res.deg loop
      declare
        mat : HexaDobl_Complex_Matrices.Matrix(1..dim,1..dim);
      begin
        for i in 1..dim loop
          for j in 1..dim loop
            mat(i,j) := HexaDobl_Random_Numbers.Random1;
          end loop;
        end loop;
        res.cff(d) := new HexaDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Random_Matrix_Series;

end HexaDobl_Random_Series_Matrices;
