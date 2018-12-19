with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;

package body QuadDobl_Complex_Matrix_Norms is

  function Max_Norm ( m : Matrix ) return quad_double is

    res : quad_double := create(0.0);
    nrm : quad_double;

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

end QuadDobl_Complex_Matrix_Norms;
