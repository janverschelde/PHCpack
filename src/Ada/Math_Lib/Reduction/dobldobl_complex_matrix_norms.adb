with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;

package body DoblDobl_Complex_Matrix_Norms is

  function Max_Norm ( m : Matrix ) return double_double is

    res : double_double := create(0.0);
    nrm : double_double;

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

end DoblDobl_Complex_Matrix_Norms;
