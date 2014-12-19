with unchecked_deallocation;

package body Generic_Laur_Poly_Systems is 

-- COPYING :

  procedure Copy ( p : in Laur_Sys; q : in out Laur_Sys ) is
  begin
    for i in p'range loop
      Copy(p(i),q(i));
    end loop;
  end Copy;

-- SELECTOR :

  function Size_of_Support ( p : Laur_Sys ) return natural32 is

    res,size : natural32 := 0;

  begin
    for i in p'range loop
      size := Size_of_Support(p(i));
      if size > res
       then res := size;
      end if;
    end loop;
    return size;
  end Size_of_Support;

-- ARITHMETIC OPERATIONS :

  function "+" ( p,q : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := p(i)+q(i);
    end loop;
    return res;
  end "+";

  function "-" ( p,q : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := p(i)-q(i);
    end loop;
    return res;
  end "-";

  function "-" ( p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := -p(i);
    end loop;
    return res;
  end "-";

  function "*" ( a : number; p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := a*p(k);
    end loop;
    return res;
  end "*";

  function "*" ( p : Laur_Sys; a : number ) return Laur_Sys is
  begin
    return a*p;
  end "*";

  procedure Add ( p : in out Laur_Sys; q : in Laur_Sys ) is
  begin
    for i in p'range loop
      Add(p(i),q(i));
    end loop;
  end Add;

  procedure Sub ( p : in out Laur_Sys; q : in Laur_Sys ) is
  begin
    for i in p'range loop
      Sub(p(i),q(i));
    end loop;
  end Sub;

  procedure Min ( p : in out Laur_Sys ) is
  begin
    for i in p'range loop
      Min(p(i));
    end loop;
  end Min;

  procedure Mul ( p : in out Laur_Sys; a : in number ) is
  begin
    for k in p'range loop
      Mul(p(k),a);
    end loop;
  end Mul;

  function Diff ( p : Laur_Sys; i : integer32 ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for j in p'range loop
      res(j) := Diff(p(j),i);
    end loop;
    return res;
  end Diff;

  procedure Diff ( p : in out Laur_Sys; i : in integer32 ) is
  begin
    for j in p'range loop
      Diff(p(j),i);
    end loop;
  end Diff;

-- DESTRUCTORS :

  procedure Clear ( p : in out Laur_Sys ) is
  begin
    for i in p'range loop
      Clear(p(i));
    end loop;
  end Clear;

  procedure Clear ( p : in out Link_to_Laur_Sys ) is

    procedure free is new unchecked_deallocation(Laur_Sys,Link_to_Laur_Sys);

  begin
    if p /= null
     then Clear(p.all);
    end if;
    free(p);
  end Clear;

end Generic_Laur_Poly_Systems;
