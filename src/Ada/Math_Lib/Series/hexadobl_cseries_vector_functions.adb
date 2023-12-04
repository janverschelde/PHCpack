with HexaDobl_Complex_Series;
with HexaDobl_Complex_Series_Functions;

package body HexaDobl_CSeries_Vector_Functions is

  function Eval ( v : HexaDobl_Complex_Series_Vectors.Vector;
                  t : hexa_double )
                return HexaDobl_Complex_Vectors.Vector is

    res : HexaDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := HexaDobl_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : HexaDobl_Complex_Series_Vectors.Vector;
                  t : Complex_Number )
                return HexaDobl_Complex_Vectors.Vector is

    res : HexaDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := HexaDobl_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : HexaDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : hexa_double )
                return HexaDobl_Complex_Vectors.Vector is

    res : HexaDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := HexaDobl_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Eval ( v : HexaDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return HexaDobl_Complex_Vectors.Vector is

    res : HexaDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := HexaDobl_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Shift ( v : HexaDobl_Complex_Series_Vectors.Vector;
                   c : hexa_double )
                 return HexaDobl_Complex_Series_Vectors.Vector is

    res : HexaDobl_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : HexaDobl_Complex_Series_Vectors.Vector;
                   c : Complex_Number )
                 return HexaDobl_Complex_Series_Vectors.Vector is

    res : HexaDobl_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := HexaDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : HexaDobl_Complex_Series_VecVecs.VecVec;
                   c : hexa_double )
                 return HexaDobl_Complex_Series_VecVecs.VecVec is

    res : HexaDobl_Complex_Series_VecVecs.VecVec(v'range);

    use HexaDobl_Complex_Series_Vectors;

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

  function Shift ( v : HexaDobl_Complex_Series_VecVecs.VecVec;
                   c : Complex_Number )
                 return HexaDobl_Complex_Series_VecVecs.VecVec is

    res : HexaDobl_Complex_Series_VecVecs.VecVec(v'range);

    use HexaDobl_Complex_Series_Vectors;

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

  procedure Shift ( v : in out HexaDobl_Complex_Series_Vectors.Vector;
                    c : in hexa_double ) is
  begin
    for i in v'range loop
      HexaDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out HexaDobl_Complex_Series_Vectors.Vector;
                    c : in Complex_Number ) is
  begin
    for i in v'range loop
      HexaDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out HexaDobl_Complex_Series_VecVecs.VecVec;
                    c : in hexa_double ) is

    use HexaDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  procedure Shift ( v : in out HexaDobl_Complex_Series_VecVecs.VecVec;
                    c : in Complex_Number ) is

    use HexaDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  function Make_Deep_Copy
             ( v : HexaDobl_Complex_Series_Vectors.Vector )
             return HexaDobl_Complex_Series_Vectors.Vector is

    res : HexaDobl_Complex_Series_Vectors.Vector(v'range);

    use HexaDobl_Complex_Series;

  begin
    for i in v'range loop
      res(i) := new Series'(Create(v(i).all,v(i).deg));
    end loop;
    return res;
  end Make_Deep_Copy;

  function Make_Deep_Copy
             ( v : HexaDobl_Complex_Series_VecVecs.VecVec )
             return HexaDobl_Complex_Series_VecVecs.VecVec is

    res : HexaDobl_Complex_Series_VecVecs.VecVec(v'range);

    use HexaDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null then
        declare
          cp : constant HexaDobl_Complex_Series_Vectors.Vector(v(i)'range)
             := Make_Deep_Copy(v(i).all);
        begin
          res(i) := new HexaDobl_Complex_Series_Vectors.Vector'(cp);
        end;
      end if;
    end loop;
    return res;
  end Make_Deep_Copy;

  procedure Deep_Clear
              ( v : in out HexaDobl_Complex_Series_Vectors.Vector ) is

    use HexaDobl_Complex_Series;

  begin
    for i in v'range loop
      if v(i) /= null
       then HexaDobl_Complex_Series.Clear(v(i));
      end if;
    end loop;
  end Deep_Clear;

  procedure Deep_Clear
              ( v : in out HexaDobl_Complex_Series_VecVecs.VecVec ) is

    use HexaDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Deep_Clear(v(i).all);
      end if;
    end loop;
    HexaDobl_Complex_Series_VecVecs.Clear(v);
  end Deep_Clear;

end HexaDobl_CSeries_Vector_Functions;
