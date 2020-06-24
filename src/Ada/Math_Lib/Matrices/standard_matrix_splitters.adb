with Standard_Complex_Numbers;
with Standard_Floating_Vectors;

package body Standard_Matrix_Splitters is

  procedure Complex_Parts
              ( mat : in Standard_Complex_Matrices.Matrix;
                rvv,ivv : in Standard_Floating_VecVecs.Link_to_VecVec ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in mat'range(2) loop -- split k-th column of the matrix
      rlnk := rvv(k);
      ilnk := ivv(k);
      for i in mat'range(1) loop
        rlnk(i) := Standard_Complex_Numbers.REAL_PART(mat(i,k));
        ilnk(i) := Standard_Complex_Numbers.IMAG_PART(mat(i,k));
      end loop;
    end loop;
  end Complex_Parts;

  procedure Complex_Merge
              ( rvv,ivv : in Standard_Floating_VecVecs.Link_to_VecVec;
                mat : out Standard_Complex_Matrices.Matrix ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for k in rvv'range loop
      rlnk := rvv(k);
      ilnk := ivv(k);
      for i in rlnk'range loop
        mat(i,k) := Standard_Complex_Numbers.Create(rlnk(i),ilnk(i));
      end loop;
    end loop;
  end Complex_Merge;

end Standard_Matrix_Splitters;
