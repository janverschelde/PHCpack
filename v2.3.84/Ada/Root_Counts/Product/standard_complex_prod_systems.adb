with unchecked_deallocation;

package body Standard_Complex_Prod_Systems is

  function Expand ( p : Prod_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Expand(p(i));
    end loop;
    return res;
  end Expand;

-- DESTRUCTORS :

  procedure free is new unchecked_deallocation(Prod_Sys,Link_to_Prod_Sys);

  procedure Shallow_Clear ( p : in out Prod_Sys ) is
  begin
    for i in p'range loop
      Shallow_Clear(p(i));
    end loop;
  end Shallow_Clear;

  procedure Shallow_Clear ( p : in out Link_to_Prod_Sys ) is
  begin
    if p /= null
     then Shallow_Clear(p.all);
          free(p);
    end if;
  end Shallow_Clear;

  procedure Deep_Clear ( p : in out Prod_Sys ) is
  begin
    for i in p'range loop
      Deep_Clear(p(i));
    end loop;
  end Deep_Clear;

  procedure Deep_Clear ( p : in out Link_to_Prod_Sys ) is
  begin
    if p /= null
     then Deep_Clear(p.all); free(p);
    end if;
  end Deep_Clear;

end Standard_Complex_Prod_Systems; 
