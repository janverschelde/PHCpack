with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;

package body Series_Coefficient_Vectors is

  function Standard_Series_Coefficients
             ( s : Standard_Complex_Series_Vectors.Vector )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant Standard_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end Standard_Series_Coefficients;

  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Complex_Series_Vectors.Vector )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant DoblDobl_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end DoblDobl_Series_Coefficients;

  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Series_Vectors.Vector )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(s'range);

  begin
    for k in s'range loop
      declare
        cff : constant QuadDobl_Complex_Vectors.Vector(0..s(k).deg)
            := s(k).cff(0..s(k).deg);
      begin
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    return res;
  end QuadDobl_Series_Coefficients;

end Series_Coefficient_Vectors;
