with OctoDobl_Complex_Vectors;
with OctoDobl_Complex_Matrices;

package body OctoDobl_Complex_Matrix_Series is

-- CONSTRUCTORS :

  function Create ( A : OctoDobl_Complex_Series_Matrices.Matrix )
                  return OctoDobl_Complex_Matrix_Series.Matrix is

    deg : constant integer32 := A(A'first(1),A'first(2)).deg;
    res : OctoDobl_Complex_Matrix_Series.Matrix(deg);

  begin
    for k in 0..deg loop
      declare
        mat : OctoDobl_Complex_Matrices.Matrix(A'range(1),A'range(2));
      begin
        for i in A'range(1) loop
          for j in A'range(2) loop
            mat(i,j) := A(i,j).cff(k);
          end loop;
        end loop;
        res.cff(k) := new OctoDobl_Complex_Matrices.Matrix'(mat);
      end;
    end loop;
    return res;
  end Create;

  function Create ( A : OctoDobl_Complex_Matrix_Series.Matrix )
                  return OctoDobl_Complex_Series_Matrices.Matrix is

    mat : OctoDobl_Complex_Matrices.Link_to_Matrix := A.cff(A.cff'first);
    rows : constant integer32 := mat'last(1);
    cols : constant integer32 := mat'last(2);
    res : OctoDobl_Complex_Series_Matrices.Matrix(1..rows,1..cols);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        declare
          s : OctoDobl_Complex_Series.Series(A.deg);
        begin
          for k in 0..A.deg loop
            mat := A.cff(k);
            s.cff(k) := mat(i,j);
          end loop;
          res(i,j) := new OctoDobl_Complex_Series.Series'(s);
        end;
      end loop;
    end loop;
    return res;
  end Create;

-- MULTIPLIER :

  function Multiply
             ( mat : OctoDobl_Complex_Matrix_Series.Matrix;
               vec : OctoDobl_Complex_Vector_Series.Vector )
             return OctoDobl_Complex_Vector_Series.Vector is

    deg : constant integer32 := mat.deg + vec.deg;
    res : OctoDobl_Complex_Vector_Series.Vector(deg);
    dim : constant integer32 := vec.cff(0)'last;
    xdg : integer32;

    use OctoDobl_Complex_Vectors;
    use OctoDobl_Complex_Matrices;

  begin
    for k in 0..mat.deg loop
      declare
        acc : OctoDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(0).all*vec.cff(k).all;
      begin
        for i in 1..k loop
          acc := acc + mat.cff(i).all*vec.cff(k-i).all;
        end loop;
        res.cff(k) := new OctoDobl_Complex_Vectors.Vector'(acc);
      end;
    end loop;
    xdg := mat.deg+1;
    for k in 1..mat.deg loop -- computed extended degree terms
      declare
        acc : OctoDobl_Complex_Vectors.Vector(1..dim)
            := mat.cff(k).all*vec.cff(xdg-k).all;
      begin
        for i in (k+1)..mat.deg loop
          acc := acc + mat.cff(i).all*vec.cff(xdg-i).all;
        end loop;
        res.cff(xdg) := new OctoDobl_Complex_Vectors.Vector'(acc);
      end;
      xdg := xdg + 1;
    end loop;
    return res;
  end Multiply;

-- DESTRUCTOR :

  procedure Clear ( A : in out Matrix ) is
  begin
    for i in 0..A.deg loop
      OctoDobl_Complex_Matrices.Clear(A.cff(i));
    end loop;
  end Clear;

end OctoDobl_Complex_Matrix_Series;
