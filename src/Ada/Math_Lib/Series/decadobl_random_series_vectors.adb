with DecaDobl_Complex_Random_Series;

package body DecaDobl_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return DecaDobl_Complex_Series_Vectors.Vector is

    res : DecaDobl_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := DecaDobl_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return DecaDobl_Complex_Vector_Series.Vector is

    rnd : DecaDobl_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant DecaDobl_Complex_Vector_Series.Vector
        := DecaDobl_Complex_Vector_Series.Create(rnd); 

  begin
    DecaDobl_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end DecaDobl_Random_Series_Vectors;
