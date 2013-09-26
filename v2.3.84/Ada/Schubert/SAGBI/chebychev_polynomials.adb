with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Chebychev_Polynomials is

  function Create ( k : natural32 ) return Vector is
  begin
    if k = 0 then
      declare
        res : Vector(0..0);
      begin
        res(0) := 1.0;
        return res;
      end;
    elsif k = 1 then
      declare
        res : Vector(0..1);
      begin
        res(0) := 0.0;
        res(1) := 1.0;
        return res;
      end;
    else
      declare
        res : Vector(0..integer32(k));
        tk1 : constant Vector := Create(k-1);
        tk2 : constant Vector := Create(k-2);
      begin
        for i in tk1'range loop
          res(i+1) := 2.0*tk1(i);
        end loop;
        res(0) := 0.0;
        for i in tk2'range loop
          res(i) := res(i) - tk2(i);
        end loop;
        return res;
      end;
    end if;
  end Create;

  function Eval ( k : natural32; x : double_float ) return double_float is

    res : constant double_float := COS(double_float(k)*ARCCOS(x));

  begin
    return res;
  end Eval;

  function Eval ( p : Vector; x : double_float ) return double_float is

    res : double_float := p(0);
    pow : double_float := x;

  begin
    for i in 1..p'last loop
      res := res + p(i)*pow;
      pow := pow*x;
    end loop;
    return res;
  end Eval;

  function Diff ( p : Vector ) return Vector is

    res : Vector(0..p'last-1);

  begin
    for i in res'range loop
      res(i) := double_float(i+1)*p(i+1);
    end loop;
    return res;
  end Diff;

  function Diff ( p : Vector; k : natural32 ) return Vector is
  begin
    if k = 0 then
      return p;
    elsif k = 1 then
      return Diff(p);
    else
      return Diff(Diff(p),k-1);
    end if;
  end Diff;

  function Int ( p : Vector ) return Vector is

    res : Vector(0..p'last+1);

  begin
    res(0) := 0.0;
    for i in p'range loop
      res(i+1) := p(i)/double_float(i+1);
    end loop;
    return res;
  end Int;

  function Int ( p : Vector; k : natural32 ) return Vector is
  begin
    if k = 0 then
      return p;
    elsif k = 1 then
      return Int(p);
    else
      return Int(Int(p),k-1);
    end if;
  end Int;

end Chebychev_Polynomials;
