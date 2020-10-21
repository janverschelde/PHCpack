with DecaDobl_Complex_Numbers_cv;

package body DecaDobl_Complex_Matrices_cv is

  function to_octo_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return OctoDobl_Complex_Matrices.Matrix is

    res : OctoDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use DecaDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DecaDobl_Complex_to_OctoDobl(A(i,j));
      end loop;
    end loop;
    return res;
  end to_octo_double;

  function to_penta_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return PentDobl_Complex_Matrices.Matrix is

    res : PentDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use DecaDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DecaDobl_Complex_to_PentDobl(A(i,j));
      end loop;
    end loop;
    return res;
  end to_penta_double;

  function to_quad_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use DecaDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DecaDobl_Complex_to_QuadDobl(A(i,j));
      end loop;
    end loop;
    return res;
  end to_quad_double;

  function to_triple_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return TripDobl_Complex_Matrices.Matrix is

    res : TripDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use DecaDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DecaDobl_Complex_to_TripDobl(A(i,j));
      end loop;
    end loop;
    return res;
  end to_triple_double;

  function to_double_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use DecaDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DecaDobl_Complex_to_DoblDobl(A(i,j));
      end loop;
    end loop;
    return res;
  end to_double_double;

  function to_double
             ( A : DecaDobl_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use DecaDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := DecaDobl_Complex_to_Standard(A(i,j));
      end loop;
    end loop;
    return res;
  end to_double;

  function to_octo_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return OctoDobl_Complex_VecMats.VecMat is

    res : OctoDobl_Complex_VecMats.VecMat(A'range);
    lnk : DecaDobl_Complex_Matrices.Link_to_Matrix;

  begin
    for i in A'range loop
      lnk := A(i);
      declare
        B : constant OctoDobl_Complex_Matrices.Matrix(lnk'range(1),lnk'range(2))
          := to_octo_double(lnk.all);
      begin
        res(i) := new OctoDobl_Complex_Matrices.Matrix'(B);
      end;
    end loop;
    return res;
  end to_octo_double;

  function to_penta_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return PentDobl_Complex_VecMats.VecMat is

    res : PentDobl_Complex_VecMats.VecMat(A'range);
    lnk : DecaDobl_Complex_Matrices.Link_to_Matrix;

  begin
    for i in A'range loop
      lnk := A(i);
      declare
        B : constant PentDobl_Complex_Matrices.Matrix(lnk'range(1),lnk'range(2))
          := to_penta_double(lnk.all);
      begin
        res(i) := new PentDobl_Complex_Matrices.Matrix'(B);
      end;
    end loop;
    return res;
  end to_penta_double;

  function to_quad_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return QuadDobl_Complex_VecMats.VecMat is

    res : QuadDobl_Complex_VecMats.VecMat(A'range);
    lnk : DecaDobl_Complex_Matrices.Link_to_Matrix;

  begin
    for i in A'range loop
      lnk := A(i);
      declare
        B : constant QuadDobl_Complex_Matrices.Matrix(lnk'range(1),lnk'range(2))
          := to_quad_double(lnk.all);
      begin
        res(i) := new QuadDobl_Complex_Matrices.Matrix'(B);
      end;
    end loop;
    return res;
  end to_quad_double;

  function to_triple_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return TripDobl_Complex_VecMats.VecMat is

    res : TripDobl_Complex_VecMats.VecMat(A'range);
    lnk : DecaDobl_Complex_Matrices.Link_to_Matrix;

  begin
    for i in A'range loop
      lnk := A(i);
      declare
        B : constant TripDobl_Complex_Matrices.Matrix(lnk'range(1),lnk'range(2))
          := to_triple_double(lnk.all);
      begin
        res(i) := new TripDobl_Complex_Matrices.Matrix'(B);
      end;
    end loop;
    return res;
  end to_triple_double;

  function to_double_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return DoblDobl_Complex_VecMats.VecMat is

    res : DoblDobl_Complex_VecMats.VecMat(A'range);
    lnk : DecaDobl_Complex_Matrices.Link_to_Matrix;

  begin
    for i in A'range loop
      lnk := A(i);
      declare
        B : constant DoblDobl_Complex_Matrices.Matrix(lnk'range(1),lnk'range(2))
          := to_double_double(lnk.all);
      begin
        res(i) := new DoblDobl_Complex_Matrices.Matrix'(B);
      end;
    end loop;
    return res;
  end to_double_double;

  function to_double
             ( A : DecaDobl_Complex_VecMats.VecMat )
             return Standard_Complex_VecMats.VecMat is

    res : Standard_Complex_VecMats.VecMat(A'range);
    lnk : DecaDobl_Complex_Matrices.Link_to_Matrix;

  begin
    for i in A'range loop
      lnk := A(i);
      declare
        B : constant Standard_Complex_Matrices.Matrix(lnk'range(1),lnk'range(2))
          := to_double(lnk.all);
      begin
        res(i) := new Standard_Complex_Matrices.Matrix'(B);
      end;
    end loop;
    return res;
  end to_double;

end DecaDobl_Complex_Matrices_cv;
