with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Plane_Representations is

  function Localize ( locmap : Standard_Natural_Matrices.Matrix;
                      plamat : Standard_Complex_Matrices.Matrix )
                    return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(plamat'range(1),plamat'range(2));
    tol : constant double_float := 10.0**(-10);
    done_j : boolean;

  begin
    for j in locmap'range(2) loop
      done_j := false;
      for i in locmap'range(1) loop
        if locmap(i,j) = 1 then
          if AbsVal(plamat(i,j)) > tol then
            for k in plamat'range(1) loop
              if AbsVal(plamat(k,j)) > tol
               then res(k,j) := plamat(k,j)/plamat(i,j);
               else res(k,j) := plamat(k,j);
              end if;
            end loop;
            res(i,j) := Create(1.0);
          end if;
          done_j := true;
        else
          res(i,j) := plamat(i,j);
        end if;
        exit when done_j;
      end loop;
    end loop;
    return res;
  end Localize;

  function Vector_Rep ( plamat : Standard_Complex_Matrices.Matrix )
                      return Standard_Complex_Vectors.Vector is 

    dim : constant integer32 := plamat'length(1)*plamat'length(2);
    res : Standard_Complex_Vectors.Vector(1..dim);
    cnt : integer32 := 0;

  begin
    for i in plamat'range(1) loop
      for j in plamat'range(2) loop
        cnt := cnt + 1;
        res(cnt) := plamat(i,j);
      end loop;
    end loop;
    return res;
  end Vector_Rep;

  function Vector_Rep ( locmap : Standard_Natural_Matrices.Matrix;
                        plamat : Standard_Complex_Matrices.Matrix )
                      return Standard_Complex_Vectors.Vector is

    dim : constant integer32 := plamat'length(1)*plamat'length(2);
    res : Standard_Complex_Vectors.Vector(1..dim);
    cnt : integer32 := 0;

  begin
    for i in plamat'range(1) loop
      for j in plamat'range(2) loop
        if locmap(i,j) = 2 then
          cnt := cnt + 1;
          res(cnt) := plamat(i,j);
        end if;
      end loop;
    end loop;
    return res(1..cnt);
  end Vector_Rep;

  function Matrix_Rep ( locmap : Standard_Natural_Matrices.Matrix;
                        plavec : Standard_Complex_Vectors.Vector )
                      return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(locmap'range(1),locmap'range(2));
    cnt : integer32 := 0;

  begin
    for i in locmap'range(1) loop
      for j in locmap'range(2) loop
        if locmap(i,j) = 0 then
          res(i,j) := Create(0.0);
        elsif locmap(i,j) = 1 then
          res(i,j) := Create(1.0);
        else
          cnt := cnt + 1;
          res(i,j) := plavec(cnt);
        end if;
      end loop;
    end loop;
    return res;
  end Matrix_Rep;

end Plane_Representations;
