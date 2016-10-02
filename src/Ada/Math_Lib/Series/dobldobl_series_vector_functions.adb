with DoblDobl_Dense_Series;

package body DoblDobl_Series_Vector_Functions is

  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  t : double_double )
                return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := DoblDobl_Dense_Series.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  t : Complex_Number )
                return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := DoblDobl_Dense_Series.Eval(v(k),t);
    end loop;
    return res;
  end Eval;

  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : double_double )
                return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := DoblDobl_Dense_Series.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

  function Eval ( v : DoblDobl_Dense_Series_Vectors.Vector;
                  w : Standard_Integer_Vectors.Vector;
                  t : Complex_Number )
                return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for k in v'range loop
      res(k) := DoblDobl_Dense_Series.Eval(v(k),t,w(k),w(w'last));
    end loop;
    return res;
  end Eval;

end DoblDobl_Series_Vector_Functions;
