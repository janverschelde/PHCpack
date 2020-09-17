with PentDobl_Complex_Random_Series;

package body PentDobl_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return PentDobl_Complex_Series_Vectors.Vector is

    res : PentDobl_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := PentDobl_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return PentDobl_Complex_Vector_Series.Vector is

    rnd : PentDobl_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant PentDobl_Complex_Vector_Series.Vector
        := PentDobl_Complex_Vector_Series.Create(rnd); 

  begin
    PentDobl_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end PentDobl_Random_Series_Vectors;
