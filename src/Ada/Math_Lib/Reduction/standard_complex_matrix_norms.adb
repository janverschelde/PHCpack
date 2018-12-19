with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Standard_Complex_Matrix_Norms is

  function Max_Norm ( m : Matrix ) return double_float is

    res : double_float := 0.0;
    nrm : double_float;

  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        nrm := AbsVal(m(i,j));
        if nrm > res
         then res := nrm;
        end if;
      end loop;
    end loop;
    return res;
  end Max_Norm;

end Standard_Complex_Matrix_Norms;
