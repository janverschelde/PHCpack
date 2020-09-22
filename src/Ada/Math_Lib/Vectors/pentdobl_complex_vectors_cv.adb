with PentDobl_Complex_Numbers_cv;        use PentDobl_Complex_Numbers_cv;

package body PentDobl_Complex_Vectors_cv is

  function Standard_to_PentDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_to_PentDobl_Complex(v(i));
    end loop;
    return res;
  end Standard_to_PentDobl_Complex;

  function Multprec_to_PentDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return PentDobl_Complex_Vectors.Vector is

    res : PentDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Multprec_to_PentDobl_Complex(v(i));
    end loop;
    return res;
  end Multprec_to_PentDobl_Complex;

  function PentDobl_Complex_to_Standard
             ( v : PentDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := PentDobl_Complex_to_Standard(v(i));
    end loop;
    return res;
  end PentDobl_Complex_to_Standard;

  function PentDobl_Complex_to_Multprec
             ( v : PentDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := PentDobl_Complex_to_Multprec(v(i));
    end loop;
    return res;
  end PentDobl_Complex_to_Multprec;

end PentDobl_Complex_Vectors_cv;
