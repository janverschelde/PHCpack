with Standard_Complex_Random_Series;

package body Standard_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return Standard_Complex_Series_Vectors.Vector is

    res : Standard_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := Standard_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return Standard_Complex_Vector_Series.Vector is

    rnd : Standard_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant Standard_Complex_Vector_Series.Vector
        := Standard_Complex_Vector_Series.Create(rnd); 

  begin
    Standard_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end Standard_Random_Series_Vectors;
