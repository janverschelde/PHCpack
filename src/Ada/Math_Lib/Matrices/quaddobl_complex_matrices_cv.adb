with QuadDobl_Complex_Numbers_cv;

package body QuadDobl_Complex_Matrices_cv is

  function to_triple_double
             ( A : QuadDobl_Complex_Matrices.Matrix )
             return TripDobl_Complex_Matrices.Matrix is

    res : TripDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use QuadDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Complex_to_TripDobl(A(i,j));
      end loop;
    end loop;
    return res;
  end to_triple_double;

  function to_double_double
             ( A : QuadDobl_Complex_Matrices.Matrix )
             return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use QuadDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Complex_to_DoblDobl(A(i,j));
      end loop;
    end loop;
    return res;
  end to_double_double;

  function to_double
             ( A : QuadDobl_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(A'range(1),A'range(2));

    use QuadDobl_Complex_Numbers_cv;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Complex_to_Standard(A(i,j));
      end loop;
    end loop;
    return res;
  end to_double;

  function to_triple_double
             ( A : QuadDobl_Complex_VecMats.VecMat )
             return TripDobl_Complex_VecMats.VecMat is

    res : TripDobl_Complex_VecMats.VecMat(A'range);
    lnk : QuadDobl_Complex_Matrices.Link_to_Matrix;

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
             ( A : QuadDobl_Complex_VecMats.VecMat )
             return DoblDobl_Complex_VecMats.VecMat is

    res : DoblDobl_Complex_VecMats.VecMat(A'range);
    lnk : QuadDobl_Complex_Matrices.Link_to_Matrix;

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
             ( A : QuadDobl_Complex_VecMats.VecMat )
             return Standard_Complex_VecMats.VecMat is

    res : Standard_Complex_VecMats.VecMat(A'range);
    lnk : QuadDobl_Complex_Matrices.Link_to_Matrix;

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

end QuadDobl_Complex_Matrices_cv;
