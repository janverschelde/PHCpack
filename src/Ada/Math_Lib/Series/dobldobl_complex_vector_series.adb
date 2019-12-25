with unchecked_deallocation;
with DoblDobl_Complex_Series;

package body DoblDobl_Complex_Vector_Series is

-- CONSTRUCTORS :

  function Create ( v : DoblDobl_Complex_Series_Vectors.Vector )
                  return DoblDobl_Complex_Vector_Series.Vector is

    dim : constant integer32 := v'last;
    deg : constant integer32 := v(v'first).deg;
    res : DoblDobl_Complex_Vector_Series.Vector(deg);

  begin
    for i in 0..res.deg loop
      res.cff(i) := new DoblDobl_Complex_Vectors.Vector(1..dim);
    end loop;
    for i in v'range loop
      for j in 0..v(i).deg loop
        res.cff(j)(i) := v(i).cff(j);
      end loop;
    end loop;
    return res;
  end Create;

  function Create ( v : DoblDobl_Complex_Vector_Series.Vector )
                  return DoblDobl_Complex_Series_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    res : DoblDobl_Complex_Series_Vectors.Vector(1..dim);

  begin
    for i in res'range loop
      res(i) := DoblDobl_Complex_Series.Create(0,v.deg);
    end loop;
    for i in 0..v.deg loop
      for j in res'range loop
        res(j).cff(i) := v.cff(i)(j);
      end loop;
    end loop;
    return res;
  end Create;

-- EVALUATORS :

  function Eval ( v : DoblDobl_Complex_Vector_Series.Vector;
                  t : double_double )
                return DoblDobl_Complex_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    cff : DoblDobl_Complex_Vectors.Link_to_Vector := v.cff(v.deg);
    res : DoblDobl_Complex_Vectors.Vector(1..dim) := cff.all;

  begin
    for i in reverse 0..(v.deg-1) loop
      cff := v.cff(i);
      for j in 1..dim loop
        res(j) := res(j)*t + cff(j);
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( v : DoblDobl_Complex_Vector_Series.Vector;
                  t : Complex_Number )
                return DoblDobl_Complex_Vectors.Vector is

    dim : constant integer32 := v.cff(v.cff'first)'last;
    cff : DoblDobl_Complex_Vectors.Link_to_Vector := v.cff(v.deg);
    res : DoblDobl_Complex_Vectors.Vector(1..dim) := cff.all;

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

  procedure Clear ( v : in out DoblDobl_Complex_Vector_Series.Vector ) is
  begin
    for i in 0..v.deg loop
      DoblDobl_Complex_Vectors.Clear(v.cff(i));
    end loop;
  end Clear;

  procedure Clear
              ( v : in out DoblDobl_Complex_Vector_Series.Link_to_Vector ) is

    procedure free is
      new unchecked_deallocation
            (DoblDobl_Complex_Vector_Series.Vector,
             DoblDobl_Complex_Vector_Series.Link_to_Vector);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

end DoblDobl_Complex_Vector_Series;
