with DoblDobl_Complex_Vectors;
with DoblDobl_Random_Vectors;

package body DoblDobl_Random_Series is

  function Random_Series ( order : integer32 ) return Series is

    cff : DoblDobl_Complex_Vectors.Vector(0..order)
        := DoblDobl_Random_Vectors.Random_Vector(0,order);

  begin
    return DoblDobl_Dense_Series.Create(cff);
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

end DoblDobl_Random_Series;
