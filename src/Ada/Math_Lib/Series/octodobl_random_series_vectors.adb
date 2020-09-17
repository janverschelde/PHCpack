with OctoDobl_Complex_Random_Series;

package body OctoDobl_Random_Series_Vectors is

  function Random_Series_Vector
             ( first,last,degree : integer32 )
             return OctoDobl_Complex_Series_Vectors.Vector is

    res : OctoDobl_Complex_Series_Vectors.Vector(first..last);

  begin
    for k in res'range loop
      res(k) := OctoDobl_Complex_Random_Series.Random_Series(degree);
    end loop;
    return res;
  end Random_Series_Vector;

  function Random_Vector_Series
             ( first,last,degree : integer32 )
             return OctoDobl_Complex_Vector_Series.Vector is

    rnd : OctoDobl_Complex_Series_Vectors.Vector(first..last)
        := Random_Series_Vector(first,last,degree);
    res : constant OctoDobl_Complex_Vector_Series.Vector
        := OctoDobl_Complex_Vector_Series.Create(rnd); 

  begin
    OctoDobl_Complex_Series_Vectors.Clear(rnd);
    return res;
  end Random_Vector_Series;

end OctoDobl_Random_Series_Vectors;
