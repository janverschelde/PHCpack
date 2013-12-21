with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;

package body DoblDobl_Complex_Vectors_cv is

  function Standard_to_DoblDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_to_DoblDobl_Complex(v(i));
    end loop;
    return res;
  end Standard_to_DoblDobl_Complex;

  function Multprec_to_DoblDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Multprec_to_DoblDobl_Complex(v(i));
    end loop;
    return res;
  end Multprec_to_DoblDobl_Complex;

  function DoblDobl_Complex_to_Standard
             ( v : DoblDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := DoblDobl_Complex_to_Standard(v(i));
    end loop;
    return res;
  end DoblDobl_Complex_to_Standard;

  function DoblDobl_Complex_to_Multprec
             ( v : DoblDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := DoblDobl_Complex_to_Multprec(v(i));
    end loop;
    return res;
  end DoblDobl_Complex_to_Multprec;

end DoblDobl_Complex_Vectors_cv;
