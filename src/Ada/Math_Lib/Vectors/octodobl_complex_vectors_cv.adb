with OctoDobl_Complex_Numbers_cv;        use OctoDobl_Complex_Numbers_cv;

package body OctoDobl_Complex_Vectors_cv is

  function Standard_to_OctoDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_to_OctoDobl_Complex(v(i));
    end loop;
    return res;
  end Standard_to_OctoDobl_Complex;

  function Multprec_to_OctoDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return OctoDobl_Complex_Vectors.Vector is

    res : OctoDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Multprec_to_OctoDobl_Complex(v(i));
    end loop;
    return res;
  end Multprec_to_OctoDobl_Complex;

  function OctoDobl_Complex_to_Standard
             ( v : OctoDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := OctoDobl_Complex_to_Standard(v(i));
    end loop;
    return res;
  end OctoDobl_Complex_to_Standard;

  function OctoDobl_Complex_to_Multprec
             ( v : OctoDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := OctoDobl_Complex_to_Multprec(v(i));
    end loop;
    return res;
  end OctoDobl_Complex_to_Multprec;

end OctoDobl_Complex_Vectors_cv;
