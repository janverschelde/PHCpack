with OctoDobl_Complex_Series;
with OctoDobl_Complex_Series_Functions;

package body OctoDobl_CSeries_Vector_Functions is

  function Eval ( v : OctoDobl_Complex_Series_Vectors.Vector;
                  t : octo_double )
                return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := OctoDobl_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : OctoDobl_Complex_Series_Vectors.Vector;
                  t : Complex_Number )
                return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := OctoDobl_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : OctoDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : octo_double )
                return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := OctoDobl_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Eval ( v : OctoDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := OctoDobl_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Shift ( v : OctoDobl_Complex_Series_Vectors.Vector;
                   c : octo_double )
                 return OctoDobl_Complex_Series_Vectors.Vector is

    res : OctoDobl_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := OctoDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : OctoDobl_Complex_Series_Vectors.Vector;
                   c : Complex_Number )
                 return OctoDobl_Complex_Series_Vectors.Vector is

    res : OctoDobl_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := OctoDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : OctoDobl_Complex_Series_VecVecs.VecVec;
                   c : octo_double )
                 return OctoDobl_Complex_Series_VecVecs.VecVec is

    res : OctoDobl_Complex_Series_VecVecs.VecVec(v'range);

    use OctoDobl_Complex_Series_Vectors;

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

  function Shift ( v : OctoDobl_Complex_Series_VecVecs.VecVec;
                   c : Complex_Number )
                 return OctoDobl_Complex_Series_VecVecs.VecVec is

    res : OctoDobl_Complex_Series_VecVecs.VecVec(v'range);

    use OctoDobl_Complex_Series_Vectors;

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

  procedure Shift ( v : in out OctoDobl_Complex_Series_Vectors.Vector;
                    c : in octo_double ) is
  begin
    for i in v'range loop
      OctoDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out OctoDobl_Complex_Series_Vectors.Vector;
                    c : in Complex_Number ) is
  begin
    for i in v'range loop
      OctoDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out OctoDobl_Complex_Series_VecVecs.VecVec;
                    c : in octo_double ) is

    use OctoDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  procedure Shift ( v : in out OctoDobl_Complex_Series_VecVecs.VecVec;
                    c : in Complex_Number ) is

    use OctoDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  function Make_Deep_Copy
             ( v : OctoDobl_Complex_Series_Vectors.Vector )
             return OctoDobl_Complex_Series_Vectors.Vector is

    res : OctoDobl_Complex_Series_Vectors.Vector(v'range);

    use OctoDobl_Complex_Series;

  begin
    for i in v'range loop
      res(i) := new Series'(Create(v(i).all,v(i).deg));
    end loop;
    return res;
  end Make_Deep_Copy;

  function Make_Deep_Copy
             ( v : OctoDobl_Complex_Series_VecVecs.VecVec )
             return OctoDobl_Complex_Series_VecVecs.VecVec is

    res : OctoDobl_Complex_Series_VecVecs.VecVec(v'range);

    use OctoDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null then
        declare
          cp : constant OctoDobl_Complex_Series_Vectors.Vector(v(i)'range)
             := Make_Deep_Copy(v(i).all);
        begin
          res(i) := new OctoDobl_Complex_Series_Vectors.Vector'(cp);
        end;
      end if;
    end loop;
    return res;
  end Make_Deep_Copy;

  procedure Deep_Clear
              ( v : in out OctoDobl_Complex_Series_Vectors.Vector ) is

    use OctoDobl_Complex_Series;

  begin
    for i in v'range loop
      if v(i) /= null
       then OctoDobl_Complex_Series.Clear(v(i));
      end if;
    end loop;
  end Deep_Clear;

  procedure Deep_Clear
              ( v : in out OctoDobl_Complex_Series_VecVecs.VecVec ) is

    use OctoDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Deep_Clear(v(i).all);
      end if;
    end loop;
    OctoDobl_Complex_Series_VecVecs.Clear(v);
  end Deep_Clear;

end OctoDobl_CSeries_Vector_Functions;
