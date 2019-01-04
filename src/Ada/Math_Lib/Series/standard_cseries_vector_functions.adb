with Standard_Complex_Series;
with Standard_Complex_Series_Functions;

package body Standard_CSeries_Vector_Functions is

  function Eval ( v : Standard_Complex_Series_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := Standard_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : Standard_Complex_Series_Vectors.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := Standard_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : Standard_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := Standard_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Eval ( v : Standard_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := Standard_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Shift ( v : Standard_Complex_Series_Vectors.Vector;
                   c : double_float )
                 return Standard_Complex_Series_Vectors.Vector is

    res : Standard_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : Standard_Complex_Series_Vectors.Vector;
                   c : Complex_Number )
                 return Standard_Complex_Series_Vectors.Vector is

    res : Standard_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : Standard_Complex_Series_VecVecs.VecVec;
                   c : double_float )
                 return Standard_Complex_Series_VecVecs.VecVec is

    res : Standard_Complex_Series_VecVecs.VecVec(v'range);

    use Standard_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null then
        declare
          svi : constant Vector := Shift(v(i).all,c);
        begin
          res(i) := new Vector'(svi); 
        end;
      end if;
    end loop;
    return res;
  end Shift;

  function Shift ( v : Standard_Complex_Series_VecVecs.VecVec;
                   c : Complex_Number )
                 return Standard_Complex_Series_VecVecs.VecVec is

    res : Standard_Complex_Series_VecVecs.VecVec(v'range);

    use Standard_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null then
        declare
          svi : constant Vector := Shift(v(i).all,c);
        begin
          res(i) := new Vector'(svi); 
        end;
      end if;
    end loop;
    return res;
  end Shift;

  procedure Shift ( v : in out Standard_Complex_Series_Vectors.Vector;
                    c : in double_float ) is
  begin
    for i in v'range loop
      Standard_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out Standard_Complex_Series_Vectors.Vector;
                    c : in Complex_Number ) is
  begin
    for i in v'range loop
      Standard_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out Standard_Complex_Series_VecVecs.VecVec;
                    c : in double_float ) is

    use Standard_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  procedure Shift ( v : in out Standard_Complex_Series_VecVecs.VecVec;
                    c : in Complex_Number ) is

    use Standard_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  function Make_Deep_Copy
             ( v : Standard_Complex_Series_Vectors.Vector )
             return Standard_Complex_Series_Vectors.Vector is

    res : Standard_Complex_Series_Vectors.Vector(v'range);

    use Standard_Complex_Series;

  begin
    for i in v'range loop
      res(i) := new Series'(Create(v(i).all,v(i).deg));
    end loop;
    return res;
  end Make_Deep_Copy;

  function Make_Deep_Copy
             ( v : Standard_Complex_Series_VecVecs.VecVec )
             return Standard_Complex_Series_VecVecs.VecVec is

    res : Standard_Complex_Series_VecVecs.VecVec(v'range);

    use Standard_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null then
        declare
          cp : constant Standard_Complex_Series_Vectors.Vector(v(i)'range)
             := Make_Deep_Copy(v(i).all);
        begin
          res(i) := new Standard_Complex_Series_Vectors.Vector'(cp);
        end;
      end if;
    end loop;
    return res;
  end Make_Deep_Copy;

  procedure Deep_Clear
              ( v : in out Standard_Complex_Series_Vectors.Vector ) is

    use Standard_Complex_Series;

  begin
    for i in v'range loop
      if v(i) /= null
       then Standard_Complex_Series.Clear(v(i));
      end if;
    end loop;
  end Deep_Clear;

  procedure Deep_Clear
              ( v : in out Standard_Complex_Series_VecVecs.VecVec ) is

    use Standard_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Deep_Clear(v(i).all);
      end if;
    end loop;
    Standard_Complex_Series_VecVecs.Clear(v);
  end Deep_Clear;

end Standard_CSeries_Vector_Functions;
