with DoblDobl_Complex_Random_Series;

package body DoblDobl_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return DoblDobl_Complex_Series_Vectors.Vector is

    res : DoblDobl_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := DoblDobl_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return DoblDobl_Complex_Vector_Series.Vector is

    rnd : DoblDobl_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant DoblDobl_Complex_Vector_Series.Vector
        := DoblDobl_Complex_Vector_Series.Create(rnd); 

  begin
    DoblDobl_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end DoblDobl_Random_Series_Vectors;
