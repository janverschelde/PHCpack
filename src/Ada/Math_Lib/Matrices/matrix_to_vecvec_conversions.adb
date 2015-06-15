with Multprec_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Multprec_Complex_Vectors;

package body Matrix_to_VecVec_Conversions is

  function mat2vv ( A : Standard_Complex_Matrices.Matrix )
                  return Standard_Complex_VecVecs.VecVec is


    res : Standard_Complex_VecVecs.VecVec(A'range(2));
    col : Standard_Complex_Vectors.Vector(A'range(1));

  begin
    for k in A'range(2) loop
      for i in col'range loop
        col(i) := A(i,k); 
      end loop;
      res(k) := new Standard_Complex_Vectors.Vector'(col);
    end loop;
    return res;
  end mat2vv;

  function mat2vv ( A : DoblDobl_Complex_Matrices.Matrix )
                  return DoblDobl_Complex_VecVecs.VecVec is

    res : DoblDobl_Complex_VecVecs.VecVec(A'range(2));
    col : DoblDobl_Complex_Vectors.Vector(A'range(1));

  begin
    for k in A'range(2) loop
      for i in col'range loop
        col(i) := A(i,k); 
      end loop;
      res(k) := new DoblDobl_Complex_Vectors.Vector'(col);
    end loop;
    return res;
  end mat2vv;

  function mat2vv ( A : QuadDobl_Complex_Matrices.Matrix )
                  return QuadDobl_Complex_VecVecs.VecVec is

    res : QuadDobl_Complex_VecVecs.VecVec(A'range(2));
    col : QuadDobl_Complex_Vectors.Vector(A'range(1));

  begin
    for k in A'range(2) loop
      for i in col'range loop
        col(i) := A(i,k); 
      end loop;
      res(k) := new QuadDobl_Complex_Vectors.Vector'(col);
    end loop;
    return res;
  end mat2vv;

  function mat2vv ( A : Multprec_Complex_Matrices.Matrix )
                  return Multprec_Complex_VecVecs.VecVec is

    res : Multprec_Complex_VecVecs.VecVec(A'range(2));

  begin
    for k in A'range(2) loop
      declare
        col : Multprec_Complex_Vectors.Vector(A'range(1));
      begin
        for i in col'range loop
          Multprec_Complex_Numbers.Copy(A(i,k),col(i)); 
        end loop;
        res(k) := new Multprec_Complex_Vectors.Vector'(col);
      end;
    end loop;
    return res;
  end mat2vv;

end Matrix_to_VecVec_Conversions;
