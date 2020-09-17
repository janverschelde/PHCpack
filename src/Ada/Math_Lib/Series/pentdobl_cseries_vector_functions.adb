with PentDobl_Complex_Series;
with PentDobl_Complex_Series_Functions;

package body PentDobl_CSeries_Vector_Functions is

  function Eval ( v : PentDobl_Complex_Series_Vectors.Vector;
                  t : penta_double )
                return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := PentDobl_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : PentDobl_Complex_Series_Vectors.Vector;
                  t : Complex_Number )
                return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := PentDobl_Complex_Series_Functions.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : PentDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : penta_double )
                return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := PentDobl_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Eval ( v : PentDobl_Complex_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := PentDobl_Complex_Series_Functions.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Shift ( v : PentDobl_Complex_Series_Vectors.Vector;
                   c : penta_double )
                 return PentDobl_Complex_Series_Vectors.Vector is

    res : PentDobl_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := PentDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : PentDobl_Complex_Series_Vectors.Vector;
                   c : Complex_Number )
                 return PentDobl_Complex_Series_Vectors.Vector is

    res : PentDobl_Complex_Series_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := PentDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
    return res;
  end Shift;

  function Shift ( v : PentDobl_Complex_Series_VecVecs.VecVec;
                   c : penta_double )
                 return PentDobl_Complex_Series_VecVecs.VecVec is

    res : PentDobl_Complex_Series_VecVecs.VecVec(v'range);

    use PentDobl_Complex_Series_Vectors;

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

  function Shift ( v : PentDobl_Complex_Series_VecVecs.VecVec;
                   c : Complex_Number )
                 return PentDobl_Complex_Series_VecVecs.VecVec is

    res : PentDobl_Complex_Series_VecVecs.VecVec(v'range);

    use PentDobl_Complex_Series_Vectors;

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

  procedure Shift ( v : in out PentDobl_Complex_Series_Vectors.Vector;
                    c : in penta_double ) is
  begin
    for i in v'range loop
      PentDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out PentDobl_Complex_Series_Vectors.Vector;
                    c : in Complex_Number ) is
  begin
    for i in v'range loop
      PentDobl_Complex_Series_Functions.Shift(v(i),c);
    end loop;
  end Shift;

  procedure Shift ( v : in out PentDobl_Complex_Series_VecVecs.VecVec;
                    c : in penta_double ) is

    use PentDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  procedure Shift ( v : in out PentDobl_Complex_Series_VecVecs.VecVec;
                    c : in Complex_Number ) is

    use PentDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Shift(v(i).all,c);
      end if;
    end loop;
  end Shift;

  function Make_Deep_Copy
             ( v : PentDobl_Complex_Series_Vectors.Vector )
             return PentDobl_Complex_Series_Vectors.Vector is

    res : PentDobl_Complex_Series_Vectors.Vector(v'range);

    use PentDobl_Complex_Series;

  begin
    for i in v'range loop
      res(i) := new Series'(Create(v(i).all,v(i).deg));
    end loop;
    return res;
  end Make_Deep_Copy;

  function Make_Deep_Copy
             ( v : PentDobl_Complex_Series_VecVecs.VecVec )
             return PentDobl_Complex_Series_VecVecs.VecVec is

    res : PentDobl_Complex_Series_VecVecs.VecVec(v'range);

    use PentDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null then
        declare
          cp : constant PentDobl_Complex_Series_Vectors.Vector(v(i)'range)
             := Make_Deep_Copy(v(i).all);
        begin
          res(i) := new PentDobl_Complex_Series_Vectors.Vector'(cp);
        end;
      end if;
    end loop;
    return res;
  end Make_Deep_Copy;

  procedure Deep_Clear
              ( v : in out PentDobl_Complex_Series_Vectors.Vector ) is

    use PentDobl_Complex_Series;

  begin
    for i in v'range loop
      if v(i) /= null
       then PentDobl_Complex_Series.Clear(v(i));
      end if;
    end loop;
  end Deep_Clear;

  procedure Deep_Clear
              ( v : in out PentDobl_Complex_Series_VecVecs.VecVec ) is

    use PentDobl_Complex_Series_Vectors;

  begin
    for i in v'range loop
      if v(i) /= null
       then Deep_Clear(v(i).all);
      end if;
    end loop;
    PentDobl_Complex_Series_VecVecs.Clear(v);
  end Deep_Clear;

end PentDobl_CSeries_Vector_Functions;
