with QuadDobl_Complex_Numbers_cv;        use QuadDobl_Complex_Numbers_cv;

package body QuadDobl_Complex_Vectors_cv is

  function Standard_to_QuadDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_to_QuadDobl_Complex(v(i));
    end loop;
    return res;
  end Standard_to_QuadDobl_Complex;

  function Multprec_to_QuadDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Multprec_to_QuadDobl_Complex(v(i));
    end loop;
    return res;
  end Multprec_to_QuadDobl_Complex;

  function QuadDobl_Complex_to_Standard
             ( v : QuadDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := QuadDobl_Complex_to_Standard(v(i));
    end loop;
    return res;
  end QuadDobl_Complex_to_Standard;

  function QuadDobl_Complex_to_Multprec
             ( v : QuadDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := QuadDobl_Complex_to_Multprec(v(i));
    end loop;
    return res;
  end QuadDobl_Complex_to_Multprec;

end QuadDobl_Complex_Vectors_cv;
