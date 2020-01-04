with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;

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

  function Standard_Series_Coefficients
             ( s : Standard_Complex_Vector_Series.Vector )
             return Standard_Complex_VecVecs.VecVec is

    res : Standard_Complex_VecVecs.VecVec(0..s.deg);

  begin
    for k in s.cff'range loop
      declare
        cf : constant Standard_Complex_Vectors.Link_to_Vector := s.cff(k);
        cp : constant Standard_Complex_Vectors.Vector(cf'range) := cf.all;
      begin
        res(k) := new Standard_Complex_Vectors.Vector'(cp);
      end;
    end loop;
    return res;
  end Standard_Series_Coefficients;

  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Complex_Vector_Series.Vector )
             return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(0..s.deg);

  begin
    for k in s.cff'range loop
      declare
        cf : constant DoblDobl_Complex_Vectors.Link_to_Vector := s.cff(k);
        cp : constant DoblDobl_Complex_Vectors.Vector(cf'range) := cf.all;
      begin
        res(k) := new DoblDobl_Complex_Vectors.Vector'(cp);
      end;
    end loop;
    return res;
  end DoblDobl_Series_Coefficients;

  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Vector_Series.Vector )
             return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(0..s.deg);

  begin
    for k in s.cff'range loop
      declare
        cf : constant QuadDobl_Complex_Vectors.Link_to_Vector := s.cff(k);
        cp : constant QuadDobl_Complex_Vectors.Vector(cf'range) := cf.all;
      begin
        res(k) := new QuadDobl_Complex_Vectors.Vector'(cp);
      end;
    end loop;
    return res;
  end QuadDobl_Series_Coefficients;

  function Standard_Series_Coefficients
             ( s : Standard_Complex_Matrix_Series.Matrix )
             return Standard_Complex_VecMats.VecMat is

    res : Standard_Complex_VecMats.VecMat(0..s.deg);

    use Standard_Complex_Matrices;

  begin
    for k in s.cff'range loop
      declare
        cf : constant Link_to_Matrix := s.cff(k);
        cp : constant Matrix(cf'range(1),cf'range(2)) := cf.all;
      begin
        res(k) := new Matrix'(cp);
      end;
    end loop;
    return res;
  end Standard_Series_Coefficients;

  function DoblDobl_Series_Coefficients
             ( s : DoblDobl_Complex_Matrix_Series.Matrix )
             return DoblDobl_Complex_VecMats.VecMat is

    res : DoblDobl_Complex_VecMats.VecMat(0..s.deg);

    use DoblDobl_Complex_Matrices;

  begin
    for k in s.cff'range loop
      declare
        cf : constant Link_to_Matrix := s.cff(k);
        cp : constant Matrix(cf'range(1),cf'range(2)) := cf.all;
      begin
        res(k) := new Matrix'(cp);
      end;
    end loop;
    return res;
  end DoblDobl_Series_Coefficients;

  function QuadDobl_Series_Coefficients
             ( s : QuadDobl_Complex_Matrix_Series.Matrix )
             return QuadDobl_Complex_VecMats.VecMat is

    res : QuadDobl_Complex_VecMats.VecMat(0..s.deg);

    use QuadDobl_Complex_Matrices;

  begin
    for k in s.cff'range loop
      declare
        cf : constant Link_to_Matrix := s.cff(k);
        cp : constant Matrix(cf'range(1),cf'range(2)) := cf.all;
      begin
        res(k) := new Matrix'(cp);
      end;
    end loop;
    return res;
  end QuadDobl_Series_Coefficients;

end Series_Coefficient_Vectors;
