with Standard_Complex_Vectors;

package body Standard_Dense_Vector_Series is

-- CONSTRUCTORS :

  function Create ( v : Standard_Dense_Series_Vectors.Vector )
                  return Standard_Dense_Vector_Series.Vector is

    res : Standard_Dense_Vector_Series.Vector;
    dim : constant integer32 := v'last;
    deg : constant integer32 := v(v'first).deg;

  begin
    res.deg := deg;
    for i in 0..res.deg loop
      res.cff(i) := new Standard_Complex_Vectors.Vector(1..dim);
    end loop;
    for i in v'range loop
      for j in 0..v(i).deg loop
        res.cff(j)(i) := v(i).cff(j);
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( v : Standard_Dense_Vector_Series.Vector )
                  return Standard_Dense_Series_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    res : Standard_Dense_Series_Vectors.Vector(1..dim);

  begin
    for i in res'range loop
      res(i).deg := v.deg;
    end loop;
    for i in 0..v.deg loop
      for j in v.cff'range loop
        res(j).cff(i) := v.cff(i)(j);
      end loop;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( v : in out Vector ) is
  begin
    for i in 0..v.deg loop
      Standard_Complex_Vectors.Clear(v.cff(i));
    end loop;
    v.deg := -1;
  end Clear;

end Standard_Dense_Vector_Series;
