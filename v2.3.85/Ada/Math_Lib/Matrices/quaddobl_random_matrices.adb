with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Random_Numbers;

package body QuadDobl_Random_Matrices is

  function Random_Matrix ( n,m : natural32 )
                         return Quad_Double_Matrices.Matrix is

    res : Quad_Double_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Random_Numbers.Random;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

  function Random_Matrix ( n,m : natural32 )
                         return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(1..integer32(n),1..integer32(m));

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := QuadDobl_Random_Numbers.Random1;
      end loop;
    end loop;
    return res;
  end Random_Matrix;

end QuadDobl_Random_Matrices;
