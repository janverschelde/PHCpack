with HexaDobl_Complex_Random_Series;

package body HexaDobl_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return HexaDobl_Complex_Series_Vectors.Vector is

    res : HexaDobl_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := HexaDobl_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return HexaDobl_Complex_Vector_Series.Vector is

    rnd : HexaDobl_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant HexaDobl_Complex_Vector_Series.Vector
        := HexaDobl_Complex_Vector_Series.Create(rnd); 

  begin
    HexaDobl_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end HexaDobl_Random_Series_Vectors;
