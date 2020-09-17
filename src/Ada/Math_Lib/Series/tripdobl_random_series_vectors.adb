with TripDobl_Complex_Random_Series;

package body TripDobl_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return TripDobl_Complex_Series_Vectors.Vector is

    res : TripDobl_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := TripDobl_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return TripDobl_Complex_Vector_Series.Vector is

    rnd : TripDobl_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant TripDobl_Complex_Vector_Series.Vector
        := TripDobl_Complex_Vector_Series.Create(rnd); 

  begin
    TripDobl_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end TripDobl_Random_Series_Vectors;
