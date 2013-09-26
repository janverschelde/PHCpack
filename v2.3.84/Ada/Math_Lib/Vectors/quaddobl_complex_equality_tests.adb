with Quad_Double_Numbers;                use Quad_Double_Numbers;

package body QuadDobl_Complex_Equality_Tests is

  function Equal ( x,y : Complex_Number; tol : double_float ) return boolean is

    dif : constant Complex_Number := x-y;
    absdif : constant quad_double := AbsVal(dif);
    res : constant boolean := (hihi_part(absdif) < tol);

  begin
    return res;
  end Equal;

  function Equal ( x,y : Vector; tol : double_float ) return boolean is
  begin
    for i in x'range loop
      if not Equal(x(i),y(i),tol)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

end QuadDobl_Complex_Equality_Tests;
