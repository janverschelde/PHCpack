with Random_Series_Generators;

package body Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return Standard_Dense_Series2_Vectors.Vector is

    res : Standard_Dense_Series2_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := Random_Series_Generators.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return Standard_Dense_Vector_Series2.Vector is

    rnd : Standard_Dense_Series2_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : Standard_Dense_Vector_Series2.Vector
        := Standard_Dense_Vector_Series2.Create(rnd); 

  begin
    Standard_Dense_Series2_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end Random_Series_Vectors;
