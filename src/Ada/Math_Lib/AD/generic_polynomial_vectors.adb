with unchecked_deallocation;

package body Generic_Polynomial_Vectors is

  procedure Clear ( v : in out Polynomial_Vector ) is
  begin
    for i in v'range loop
      Polynomials.Clear(v(i));
    end loop;
  end Clear;

  procedure Clear ( v : in out Link_to_Polynomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Polynomial_Vector,Link_to_Polynomial_Vector);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

end Generic_Polynomial_Vectors;
