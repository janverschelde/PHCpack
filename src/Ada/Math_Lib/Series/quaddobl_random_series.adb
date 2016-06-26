with QuadDobl_Complex_Vectors;
with QuadDobl_Random_Vectors;

package body QuadDobl_Random_Series is

  function Random_Series ( order : integer32 ) return Series is

    cff : QuadDobl_Complex_Vectors.Vector(0..order)
        := QuadDobl_Random_Vectors.Random_Vector(0,order);

  begin
    return QuadDobl_Dense_Series.Create(cff);
  end Random_Series;

  function Random_Series_Vector
             ( first,last,order : integer32 ) return Vector is

    res : Vector(first..last);

  begin
    for k in res'range loop
      res(k) := Random_Series(order);
    end loop;
    return res;
  end Random_Series_Vector;

end QuadDobl_Random_Series;
