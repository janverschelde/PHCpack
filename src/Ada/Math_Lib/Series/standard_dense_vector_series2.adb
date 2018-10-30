with unchecked_deallocation;

package body Standard_Dense_Vector_Series2 is

-- CONSTRUCTORS :

  function Create ( v : Standard_Dense_Series2_Vectors.Vector )
                  return Standard_Dense_Vector_Series2.Vector is

    dim : constant integer32 := v'last;
    deg : constant integer32 := v(v'first).deg;
    res : Standard_Dense_Vector_Series2.Vector(deg);

  begin
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

  function Create ( v : Standard_Dense_Vector_Series2.Vector )
                  return Standard_Dense_Series2_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    res : Standard_Dense_Series2_Vectors.Vector(1..dim);

  begin
    for i in res'range loop
      res(i) := Standard_Dense_Series2.Create(0,v.deg);
    end loop;
    for i in 0..v.deg loop
      for j in res'range loop
        res(j).cff(i) := v.cff(i)(j);
      end loop;
    end loop;
    return res;
  end Create;

-- EVALUATORS :

  function Eval ( v : Standard_Dense_Vector_Series2.Vector;
                  t : double_float )
                return Standard_Complex_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    cff : Standard_Complex_Vectors.Link_to_Vector := v.cff(v.deg);
    res : Standard_Complex_Vectors.Vector(1..dim) := cff.all;

  begin
    for i in reverse 0..(v.deg-1) loop
      cff := v.cff(i);
      for j in 1..dim loop
        res(j) := res(j)*t + cff(j);
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( v : Standard_Dense_Vector_Series2.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    cff : Standard_Complex_Vectors.Link_to_Vector := v.cff(v.deg);
    res : Standard_Complex_Vectors.Vector(1..dim) := cff.all;

  begin
    for i in reverse 0..(v.deg-1) loop
      cff := v.cff(i);
      for j in 1..dim loop
        res(j) := res(j)*t + cff(j);
      end loop;
    end loop;
    return res;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( v : in out Standard_Dense_Vector_Series2.Vector ) is
  begin
    for i in 0..v.deg loop
      Standard_Complex_Vectors.Clear(v.cff(i));
    end loop;
  end Clear;

  procedure Clear
              ( v : in out Standard_Dense_Vector_Series2.Link_to_Vector ) is

    procedure free is
      new unchecked_deallocation
            (Standard_Dense_Vector_Series2.Vector,
             Standard_Dense_Vector_Series2.Link_to_Vector);


  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

end Standard_Dense_Vector_Series2;
