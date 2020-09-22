with TripDobl_Complex_Numbers_cv;        use TripDobl_Complex_Numbers_cv;

package body TripDobl_Complex_Vectors_cv is

  function Standard_to_TripDobl_Complex
             ( v : Standard_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector is

    res : TripDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Standard_to_TripDobl_Complex(v(i));
    end loop;
    return res;
  end Standard_to_TripDobl_Complex;

  function Multprec_to_TripDobl_Complex
             ( v : Multprec_Complex_Vectors.Vector )
             return TripDobl_Complex_Vectors.Vector is

    res : TripDobl_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := Multprec_to_TripDobl_Complex(v(i));
    end loop;
    return res;
  end Multprec_to_TripDobl_Complex;

  function TripDobl_Complex_to_Standard
             ( v : TripDobl_Complex_Vectors.Vector )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := TripDobl_Complex_to_Standard(v(i));
    end loop;
    return res;
  end TripDobl_Complex_to_Standard;

  function TripDobl_Complex_to_Multprec
             ( v : TripDobl_Complex_Vectors.Vector )
             return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      res(i) := TripDobl_Complex_to_Multprec(v(i));
    end loop;
    return res;
  end TripDobl_Complex_to_Multprec;

end TripDobl_Complex_Vectors_cv;
