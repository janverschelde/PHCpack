with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Complex_Numbers;            use Standard_Complex_Numbers;
with Standard_Integer_Vectors;

package body Standard_Monomial_Map_Substitutors is

  function Filter ( p : Poly; tol : double_float ) return Poly is

    res : Poly := Null_Poly;

    procedure Filter_Term ( t : in Term; continue : out boolean ) is
    begin
      if AbsVal(t.cf) > tol
       then Add(res,t);
      end if;
      continue := true;
    end Filter_Term;
    procedure Filter_Terms is new Visiting_Iterator(Filter_Term);

  begin
    Filter_Terms(p);
    return res;
  end Filter;

  function Filter ( p : Laur_Sys; tol : double_float ) return Laur_Sys is

    res : Laur_Sys(p'range);
    f : Poly;
    cnt : integer32 := p'first - 1;

  begin
    for i in p'range loop
      f := Filter(p(i),tol);
      if f /= Null_Poly then
        cnt := cnt + 1;
        res(cnt) := f;
      end if;
    end loop;
    return res(res'first..cnt);
  end Filter;

  function Subs ( t : Term; map : Monomial_Map ) return Term is

    res : Term;

  begin
    res.cf := t.cf;
    for i in t.dg'range loop     -- compute (map.c)**t.dg
      if t.dg(i) < 0 then
        for j in 1..(-t.dg(i)) loop
          res.cf := res.cf/map.c(i);
        end loop;
      elsif t.dg(i) > 0 then
        for j in 1..t.dg(i) loop
          res.cf := res.cf*map.c(i);
        end loop;
      end if;
    end loop;
    res.dg := new Standard_Integer_Vectors.Vector'(1..map.d => 0);
    for i in 1..map.d loop
      for j in t.dg'range loop
        res.dg(i) := res.dg(i) + t.dg(j)*map.v(j)(i);
      end loop;
    end loop;
    return res;
  end Subs;

  function Subs ( p : Poly; map : Monomial_Map ) return Poly is

    res : Poly := Null_Poly;

    procedure Subs_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Subs(t,map);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Subs_Term;
    procedure Subs_Terms is new Visiting_Iterator(Subs_Term);

  begin
    Subs_Terms(p);
    return res;
  end Subs;

end Standard_Monomial_Map_Substitutors;
