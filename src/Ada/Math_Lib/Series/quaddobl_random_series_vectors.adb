with QuadDobl_Complex_Random_Series;

package body QuadDobl_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return QuadDobl_Complex_Series_Vectors.Vector is

    res : QuadDobl_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := QuadDobl_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return QuadDobl_Complex_Vector_Series.Vector is

    rnd : QuadDobl_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant QuadDobl_Complex_Vector_Series.Vector
        := QuadDobl_Complex_Vector_Series.Create(rnd); 

  begin
    QuadDobl_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end QuadDobl_Random_Series_Vectors;
