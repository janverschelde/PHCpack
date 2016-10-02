with QuadDobl_Dense_Series;

package body QuadDobl_Series_Vector_Functions is

  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  t : quad_double )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := QuadDobl_Dense_Series.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  t : Complex_Number )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := QuadDobl_Dense_Series.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : quad_double )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := QuadDobl_Dense_Series.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Eval ( v : QuadDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := QuadDobl_Dense_Series.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

end QuadDobl_Series_Vector_Functions;
