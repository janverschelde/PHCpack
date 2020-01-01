with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;

package body DoblDobl_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : constant DoblDobl_Complex_Vectors.Vector(0..degree)
        := DoblDobl_Random_Vectors.Random_Vector(0,degree);

  begin
    return DoblDobl_Dense_Series.Create(cff);
  end Random_Series;

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return DoblDobl_Dense_Series_Vectors.Vector is

    res : DoblDobl_Dense_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return DoblDobl_Dense_Vector_Series.Vector is

    rnd : constant DoblDobl_Dense_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant DoblDobl_Dense_Vector_Series.Vector
        := DoblDobl_Dense_Vector_Series.Create(rnd); 

  begin
    return res;
  end Random_Vector_Series;

  function Random_Series_VecVec
             ( vvfirst,vvlast,first,last,degree : integer32 ) return VecVec is

    res : VecVec(vvfirst..vvlast);

  begin
    for i in res'range loop
      declare
        v : constant DoblDobl_Dense_Series_Vectors.Vector(first..last)
          := Random_Series_Vector(first,last,degree);
      begin
        res(i) := new DoblDobl_Dense_Series_Vectors.Vector'(v);
      end;
    end loop;
    return res;
  end Random_Series_VecVec;

  function Random_Series_Matrix
             ( rowfirst,rowlast,columnfirst,columnlast,degree : integer32 )
             return Matrix is


    res : Matrix(rowfirst..rowlast,columnfirst..columnlast);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Random_Series(degree);
      end loop;
    end loop;
    return res;
  end Random_Series_Matrix;

end DoblDobl_Random_Series;
