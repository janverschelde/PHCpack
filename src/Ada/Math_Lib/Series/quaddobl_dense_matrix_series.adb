with QuadDobl_Complex_Matrices;

package body QuadDobl_Dense_Matrix_Series is

-- CONSTRUCTORS :

  function Create ( A : QuadDobl_Dense_Series_Matrices.Matrix )
                  return QuadDobl_Dense_Matrix_Series.Matrix is

    res : QuadDobl_Dense_Matrix_Series.Matrix;
    deg : constant integer32 := A(A'first(1),A'first(2)).deg;

  begin
    res.deg := deg;
    for k in 0..res.deg loop
      declare
        mat : QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));
      begin
        for i in A'range(1) loop
          for j in A'range(2) loop
            mat(i,j) := A(i,j).cff(k);
          end loop;
        end loop;
        res.cff(k) := new QuadDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Create;

  function Create ( A : QuadDobl_Dense_Matrix_Series.Matrix )
                  return QuadDobl_Dense_Series_Matrices.Matrix is

    mat : QuadDobl_Complex_Matrices.Link_to_Matrix := A.cff(A.cff'first);
    rows : constant integer32 := mat'last(1);
    cols : constant integer32 := mat'last(2);
    res : QuadDobl_Dense_Series_Matrices.Matrix(1..rows,1..cols);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j).deg := A.deg;
        for k in 0..A.deg loop
          mat := A.cff(k);
          res(i,j).cff(k) := mat(i,j);
        end loop;
      end loop;
    end loop;
    return res;
  end Create;

-- DESTRUCTOR :

  procedure Clear ( A : in out Matrix ) is
  begin
    for i in 0..A.deg loop
      QuadDobl_Complex_Matrices.Clear(A.cff(i));
    end loop;
    A.deg := -1;
  end Clear;

end QuadDobl_Dense_Matrix_Series;
