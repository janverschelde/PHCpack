with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;

package body DoblDobl_Dense_Matrix_Series is

-- CONSTRUCTORS :

  function Create ( A : DoblDobl_Dense_Series_Matrices.Matrix )
                  return DoblDobl_Dense_Matrix_Series.Matrix is

    res : DoblDobl_Dense_Matrix_Series.Matrix;
    deg : constant integer32 := A(A'first(1),A'first(2)).deg;

  begin
    res.deg := deg;
    for k in 0..res.deg loop
      declare
        mat : DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));
      begin
        for i in A'range(1) loop
          for j in A'range(2) loop
            mat(i,j) := A(i,j).cff(k);
          end loop;
        end loop;
        res.cff(k) := new DoblDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Create;

  function Create ( A : DoblDobl_Dense_Matrix_Series.Matrix )
                  return DoblDobl_Dense_Series_Matrices.Matrix is

    mat : DoblDobl_Complex_Matrices.Link_to_Matrix := A.cff(A.cff'first);
    rows : constant integer32 := mat'last(1);
    cols : constant integer32 := mat'last(2);
    res : DoblDobl_Dense_Series_Matrices.Matrix(1..rows,1..cols);

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

-- MULTIPLIER :

  function Multiply
             ( mat : DoblDobl_Dense_Matrix_Series.Matrix;
               vec : DoblDobl_Dense_Vector_Series.Vector )
             return DoblDobl_Dense_Vector_Series.Vector is

    res : DoblDobl_Dense_Vector_Series.Vector;
    deg : constant integer32 := mat.deg + vec.deg;
    dim : constant integer32 := vec.cff(0)'last;
    xdg : integer32;

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

  begin
    if deg > DoblDobl_Dense_Series.max_deg
     then res.deg := DoblDobl_Dense_Series.max_deg;
     else res.deg := deg;
    end if;
    for k in 0..mat.deg loop
      declare
        acc : DoblDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(0).all*vec.cff(k).all;
      begin
        for i in 1..k loop
          acc := acc + mat.cff(i).all*vec.cff(k-i).all;
        end loop;
        res.cff(k) := new DoblDobl_Complex_Vectors.Vector'(acc);
      end;
    end loop;
    xdg := mat.deg+1;
    for k in 1..mat.deg loop -- computed extended degree terms
      declare
        acc : DoblDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(k).all*vec.cff(xdg-k).all;
      begin
        for i in (k+1)..mat.deg loop
          acc := acc + mat.cff(i).all*vec.cff(xdg-i).all;
        end loop;
        res.cff(xdg) := new DoblDobl_Complex_Vectors.Vector'(acc);
      end;
      xdg := xdg + 1;
    end loop;
    return res;
  end Multiply;

-- DESTRUCTOR :

  procedure Clear ( A : in out Matrix ) is
  begin
    for i in 0..A.deg loop
      DoblDobl_Complex_Matrices.Clear(A.cff(i));
    end loop;
    A.deg := -1;
  end Clear;

end DoblDobl_Dense_Matrix_Series;
