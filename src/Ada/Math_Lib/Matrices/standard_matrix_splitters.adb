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

  procedure Split_Rows
              ( A : in Standard_Complex_Matrices.Link_to_Matrix;
                rArows : in Standard_Floating_VecVecs.Link_to_VecVec;
                iArows : in Standard_Floating_VecVecs.Link_to_VecVec ) is

    rlnk,ilnk : Standard_Floating_Vectors.Link_to_Vector;

  begin
    for i in A'range(1) loop
      rlnk := rArows(i);
      ilnk := iArows(i);
      for j in A'range(2) loop
        rlnk(j) := Standard_Complex_Numbers.REAL_PART(A(i,j));
        ilnk(j) := Standard_Complex_Numbers.IMAG_PART(A(i,j));
      end loop;
    end loop;
  end Split_Rows;

  procedure Split_Rows
              ( vm : in Standard_Complex_VecMats.VecMat;
                rv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec;
                iv : in Standard_Floating_VecVecVecs.Link_to_VecVecVec ) is
  begin
    for k in rv'range loop
      Split_Rows(vm(k),rv(k),iv(k));
    end loop;
  end Split_Rows;

end Standard_Matrix_Splitters;
