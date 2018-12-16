with unchecked_deallocation;

package body Generic_Monomial_Vectors is

  procedure Clear ( v : in out Monomial_Vector ) is
  begin
    for i in v'range loop
      Monomials.Clear(v(i));
    end loop;
  end Clear;

  procedure Clear ( v : in out Link_to_Monomial_Vector ) is

    procedure free is
      new unchecked_deallocation(Monomial_Vector,Link_to_Monomial_Vector);

  begin
    if v /= null then
      Clear(v.all);
      free(v);
    end if;
  end Clear;

end Generic_Monomial_Vectors;
