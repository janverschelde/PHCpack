with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;

package body DoblDobl_Random_Series is

  function Random_Series ( degree : integer32 ) return Series is

    cff : DoblDobl_Complex_Vectors.Vector(0..degree)
        := DoblDobl_Random_Vectors.Random_Vector(0,degree);

  begin
    return DoblDobl_Dense_Series.Create(cff);
  end Random_Series;

  function Random_Series_Vector
             ( first,last,degree : integer32 ) return Vector is

    res : Vector(first..last);

  begin
    for k in res'range loop
      res(k) := Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Series_VecVec
             ( vvfirst,vvlast,first,last,degree : integer32 ) return VecVec is

    res : VecVec(vvfirst..vvlast);

  begin
    for i in res'range loop
      declare
        v : constant Vector(first..last)
          := Random_Series_Vector(first,last,degree);
      begin
        res(i) := new Vector'(v);
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
