with unchecked_deallocation;
with QuadDobl_Complex_Series;

package body QuadDobl_Complex_Vector_Series is

-- CONSTRUCTORS :

  function Create ( v : QuadDobl_Complex_Series_Vectors.Vector )
                  return QuadDobl_Complex_Vector_Series.Vector is

    dim : constant integer32 := v'last;
    deg : constant integer32 := v(v'first).deg;
    res : QuadDobl_Complex_Vector_Series.Vector(deg);

  begin
    for i in 0..res.deg loop
      res.cff(i) := new QuadDobl_Complex_Vectors.Vector(1..dim);
    end loop;
    for i in v'range loop
      for j in 0..v(i).deg loop
        res.cff(j)(i) := v(i).cff(j);
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( v : QuadDobl_Complex_Vector_Series.Vector )
                  return QuadDobl_Complex_Series_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    res : QuadDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    for i in res'range loop
      res(i) := QuadDobl_Complex_Series.Create(0,v.deg);
    end loop;
    for i in 0..v.deg loop
      for j in res'range loop
        res(j).cff(i) := v.cff(i)(j);
      end loop;
    end loop;
    return res;
  end Create;

-- EVALUATORS :

  function Eval ( v : QuadDobl_Complex_Vector_Series.Vector;
                  t : quad_double )
                return QuadDobl_Complex_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    cff : QuadDobl_Complex_Vectors.Link_to_Vector := v.cff(v.deg);
    res : QuadDobl_Complex_Vectors.Vector(1..dim) := cff.all;

  begin
    for i in reverse 0..(v.deg-1) loop
      cff := v.cff(i);
      for j in 1..dim loop
        res(j) := res(j)*t + cff(j);
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( v : QuadDobl_Complex_Vector_Series.Vector;
                  t : Complex_Number )
                return QuadDobl_Complex_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    cff : QuadDobl_Complex_Vectors.Link_to_Vector := v.cff(v.deg);
    res : QuadDobl_Complex_Vectors.Vector(1..dim) := cff.all;

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

  procedure Clear ( v : in out QuadDobl_Complex_Vector_Series.Vector ) is
  begin
    for i in 0..v.deg loop
      QuadDobl_Complex_Vectors.Clear(v.cff(i));
    end loop;
  end Clear;

  procedure Clear
              ( v : in out QuadDobl_Complex_Vector_Series.Link_to_Vector ) is

    procedure free is
      new unchecked_deallocation
            (QuadDobl_Complex_Vector_Series.Vector,
             QuadDobl_Complex_Vector_Series.Link_to_Vector);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

end QuadDobl_Complex_Vector_Series;
