with Standard_Complex_Vectors;
with Standard_Random_Vectors;

package body Standard_Random_Series is

  function Random_Series ( order : integer32 ) return Series is

    cff : Standard_Complex_Vectors.Vector(0..order)
        := Standard_Random_Vectors.Random_Vector(0,order);

  begin
    return Create(cff);
  end Random_Series;

  function Random_Series_Vector
             ( first,last,order : integer32 )
             return Standard_Dense_Series_Vectors.Vector is

    res : Standard_Dense_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := Random_Series(order);
    end loop;
    return res;
  end Random_Series_Vector;

end Standard_Random_Series;
