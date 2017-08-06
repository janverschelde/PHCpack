with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;

package body Polynomial_Drops is

  function Drop ( t : Standard_Complex_Polynomials.Term; k : integer32 )
                return Standard_Complex_Polynomials.Term is
  begin
    if k < t.dg'first or k > t.dg'last then
      return t;
    else
      declare
        r : Standard_Complex_Polynomials.Term;
      begin
        r.cf := t.cf;
        r.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
        for i in t.dg'first..(k-1) loop
          r.dg(i) := t.dg(i);
        end loop;
        for i in (k+1)..t.dg'last loop
          r.dg(i-1) := t.dg(i);
        end loop;
        return r;
      end;
    end if;
  end Drop;

  function Drop ( t : Standard_Complex_Laurentials.Term; k : integer32 )
                return Standard_Complex_Laurentials.Term is
  begin
    if k < t.dg'first or k > t.dg'last then
      return t;
    else
      declare
        r : Standard_Complex_Laurentials.Term;
      begin
        r.cf := t.cf;
        r.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-1);
        for i in t.dg'first..(k-1) loop
          r.dg(i) := t.dg(i);
        end loop;
        for i in (k+1)..t.dg'last loop
          r.dg(i-1) := t.dg(i);
        end loop;
        return r;
      end;
    end if;
  end Drop;

  function Drop ( t : DoblDobl_Complex_Polynomials.Term; k : integer32 )
                return DoblDobl_Complex_Polynomials.Term is
  begin
    if k < t.dg'first or k > t.dg'last then
      return t;
    else
      declare
        r : DoblDobl_Complex_Polynomials.Term;
      begin
        r.cf := t.cf;
        r.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
        for i in t.dg'first..(k-1) loop
          r.dg(i) := t.dg(i);
        end loop;
        for i in (k+1)..t.dg'last loop
          r.dg(i-1) := t.dg(i);
        end loop;
        return r;
      end;
    end if;
  end Drop;

  function Drop ( t : DoblDobl_Complex_Laurentials.Term; k : integer32 )
                return DoblDobl_Complex_Laurentials.Term is
  begin
    if k < t.dg'first or k > t.dg'last then
      return t;
    else
      declare
        r : DoblDobl_Complex_Laurentials.Term;
      begin
        r.cf := t.cf;
        r.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-1);
        for i in t.dg'first..(k-1) loop
          r.dg(i) := t.dg(i);
        end loop;
        for i in (k+1)..t.dg'last loop
          r.dg(i-1) := t.dg(i);
        end loop;
        return r;
      end;
    end if;
  end Drop;

  function Drop ( t : QuadDobl_Complex_Polynomials.Term; k : integer32 )
                return QuadDobl_Complex_Polynomials.Term is
  begin
    if k < t.dg'first or k > t.dg'last then
      return t;
    else
      declare
        r : QuadDobl_Complex_Polynomials.Term;
      begin
        r.cf := t.cf;
        r.dg := new Standard_Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
        for i in t.dg'first..(k-1) loop
          r.dg(i) := t.dg(i);
        end loop;
        for i in (k+1)..t.dg'last loop
          r.dg(i-1) := t.dg(i);
        end loop;
        return r;
      end;
    end if;
  end Drop;

  function Drop ( t : QuadDobl_Complex_Laurentials.Term; k : integer32 )
                return QuadDobl_Complex_Laurentials.Term is
  begin
    if k < t.dg'first or k > t.dg'last then
      return t;
    else
      declare
        r : QuadDobl_Complex_Laurentials.Term;
      begin
        r.cf := t.cf;
        r.dg := new Standard_Integer_Vectors.Vector(t.dg'first..t.dg'last-1);
        for i in t.dg'first..(k-1) loop
          r.dg(i) := t.dg(i);
        end loop;
        for i in (k+1)..t.dg'last loop
          r.dg(i-1) := t.dg(i);
        end loop;
        return r;
      end;
    end if;
  end Drop;

  function Drop ( p : Standard_Complex_Polynomials.Poly; k : integer32 )
                return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Drop_Term ( t : in Term; continue : out boolean ) is

      r : Term;

    begin
      if t.dg(k) = 0 then
        r := Drop(t,k);
        Add(res,r);
        Clear(r);
      end if;
      continue := true;
    end Drop_Term;
    procedure Drop_Terms is new Visiting_Iterator(Drop_Term);

  begin
    Drop_Terms(p);
    return res;
  end Drop;

  function Drop ( p : Standard_Complex_Laurentials.Poly; k : integer32 )
                return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Drop_Term ( t : in Term; continue : out boolean ) is

      r : Term;

    begin
      if t.dg(k) = 0 then
        r := Drop(t,k);
        Add(res,r);
        Clear(r);
      end if;
      continue := true;
    end Drop_Term;
    procedure Drop_Terms is new Visiting_Iterator(Drop_Term);

  begin
    Drop_Terms(p);
    return res;
  end Drop;

  function Drop ( p : DoblDobl_Complex_Polynomials.Poly; k : integer32 )
                return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Drop_Term ( t : in Term; continue : out boolean ) is

      r : Term;

    begin
      if t.dg(k) = 0 then
        r := Drop(t,k);
        Add(res,r);
        Clear(r);
      end if;
      continue := true;
    end Drop_Term;
    procedure Drop_Terms is new Visiting_Iterator(Drop_Term);

  begin
    Drop_Terms(p);
    return res;
  end Drop;

  function Drop ( p : DoblDobl_Complex_Laurentials.Poly; k : integer32 )
                return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Drop_Term ( t : in Term; continue : out boolean ) is

      r : Term;

    begin
      if t.dg(k) = 0 then
        r := Drop(t,k);
        Add(res,r);
        Clear(r);
      end if;
      continue := true;
    end Drop_Term;
    procedure Drop_Terms is new Visiting_Iterator(Drop_Term);

  begin
    Drop_Terms(p);
    return res;
  end Drop;

  function Drop ( p : QuadDobl_Complex_Polynomials.Poly; k : integer32 )
                return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;
    res : Poly := Null_Poly;

    procedure Drop_Term ( t : in Term; continue : out boolean ) is

      r : Term;

    begin
      if t.dg(k) = 0 then
        r := Drop(t,k);
        Add(res,r);
        Clear(r);
      end if;
      continue := true;
    end Drop_Term;
    procedure Drop_Terms is new Visiting_Iterator(Drop_Term);

  begin
    Drop_Terms(p);
    return res;
  end Drop;

  function Drop ( p : QuadDobl_Complex_Laurentials.Poly; k : integer32 )
                return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Drop_Term ( t : in Term; continue : out boolean ) is

      r : Term;

    begin
      if t.dg(k) = 0 then
        r := Drop(t,k);
        Add(res,r);
        Clear(r);
      end if;
      continue := true;
    end Drop_Term;
    procedure Drop_Terms is new Visiting_Iterator(Drop_Term);

  begin
    Drop_Terms(p);
    return res;
  end Drop;

  function Drop ( p : Standard_Complex_Poly_Systems.Poly_Sys; k : integer32 )
                return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Drop(p(i),k);
    end loop;
    return res;
  end Drop;

  function Drop ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys; k : integer32 )
                return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Drop(p(i),k);
    end loop;
    return res;
  end Drop;

  function Drop ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys; k : integer32 )
                return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Drop(p(i),k);
    end loop;
    return res;
  end Drop;

  function Drop ( p : Standard_Complex_Laur_Systems.Laur_Sys; k : integer32 )
                return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Drop(p(i),k);
    end loop;
    return res;
  end Drop;

  function Drop ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys; k : integer32 )
                return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Drop(p(i),k);
    end loop;
    return res;
  end Drop;

  function Drop ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys; k : integer32 )
                return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Drop(p(i),k);
    end loop;
    return res;
  end Drop;

  function Remove_Variable
             ( p : Standard_Complex_Laurentials.Poly; k : integer32 )
             return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Drop_Term ( t : in Term; continue : out boolean ) is

      r : Term;

    begin
      r := Drop(t,k);
      Add(res,r);
      Clear(r);
      continue := true;
    end Drop_Term;
    procedure Drop_Terms is new Visiting_Iterator(Drop_Term);

  begin
    Drop_Terms(p);
    return res;
  end Remove_Variable;

  function Remove_Variable
             ( p : Standard_Complex_Laur_Systems.Laur_Sys; k : integer32 )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Variable(p(i),k);
    end loop;
    return res;
  end Remove_Variable;

end Polynomial_Drops;
