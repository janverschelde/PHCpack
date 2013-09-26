package body Standard_Initial_Forms is

  function Degree ( t : Standard_Complex_Polynomials.Term; v : Vector )
                  return integer32 is

    res : integer32 := 0;
 
  begin
    for i in v'range loop
      res := res + integer32(t.dg(i))*v(i);
    end loop;
    return res;
  end Degree;

  function Degree ( t : Standard_Complex_Laurentials.Term; v : Vector )
                  return integer32 is

    res : integer32 := 0;
 
  begin
    for i in v'range loop
      res := res + t.dg(i)*v(i);
    end loop;
    return res;
  end Degree;

  function Degree ( p : Standard_Complex_Polynomials.Poly; v : Vector )
                  return integer32 is

    use Standard_Complex_Polynomials;

    res : integer32 := 0;
    first : boolean := true;

    procedure Degree_Term ( t : Term; continue : out boolean ) is

      dt : constant integer32 := Degree(t,v);

    begin
      if first then
        res := dt;
        first := false;
      elsif dt < res then
        res := dt;
      end if;
      continue := true;
    end Degree_Term;
    procedure Degree_Terms is new Visiting_Iterator(Degree_Term);

  begin
    Degree_Terms(p);
    return res;
  end Degree;

  function Degree ( p : Standard_Complex_Laurentials.Poly; v : Vector ) 
                  return integer32 is

    res : integer32 := 0;
    first : boolean := true;

    use Standard_Complex_Laurentials;

    procedure Degree_Term ( t : Term; continue : out boolean ) is

      dt : constant integer32 := Degree(t,v);

    begin
      if first then
        res := dt;
        first := false;
      elsif dt < res then
        res := dt;
      end if;
      continue := true;
    end Degree_Term;
    procedure Degree_Terms is new Visiting_Iterator(Degree_Term);

  begin
    Degree_Terms(p);
    return res;
  end Degree;

  function Form ( p : Standard_Complex_Polynomials.Poly;
                  v : Vector; d : integer32 )
                return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Select_Term ( t : Term; continue : out boolean ) is
    begin
      if Degree(t,v) = d
       then Add(res,t);
      end if;
      continue := true;
    end Select_Term;
    procedure Select_Terms is new Visiting_Iterator(Select_Term);

  begin
    Select_Terms(p);
    return res;
  end Form;

  function Form ( p : Standard_Complex_Laurentials.Poly;
                  v : Vector; d : integer32 )
                return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Select_Term ( t : Term; continue : out boolean ) is
    begin
      if Degree(t,v) = d
       then Add(res,t);
      end if;
      continue := true;
    end Select_Term;
    procedure Select_Terms is new Visiting_Iterator(Select_Term);

  begin
    Select_Terms(p);
    return res;
  end Form;

  function Initial ( p : Standard_Complex_Polynomials.Poly; v : Vector )
                   return Standard_Complex_Polynomials.Poly is

    d : constant integer32 := Degree(p,v);

  begin
    return Form(p,v,d);
  end Initial;

  function Initial ( p : Standard_Complex_Laurentials.Poly; v : Vector )
                   return Standard_Complex_Laurentials.Poly is

    d : constant integer32 := Degree(p,v);

  begin
    return Form(p,v,d);
  end Initial;

  function Initial ( p : Poly_Sys; v : Vector ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Initial(p(i),v);
    end loop;
    return res;
  end Initial;

  function Initial ( s : Laur_Sys; v : Vector ) return Laur_Sys is

    res : Laur_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Initial(s(i),v);
    end loop;
    return res;
  end Initial;

  function Transform ( t : Standard_Complex_Laurentials.Term; m : Matrix ) 
                     return Standard_Complex_Laurentials.Term is

    use Standard_Complex_Laurentials;
    res : Term;
    e : constant Vector(t.dg'range) := m*t.dg.all;

  begin
    res.cf := t.cf;
    res.dg := new Vector'(e);
    return res;
  end Transform;

  function Transform ( p : Standard_Complex_Laurentials.Poly; m : Matrix )
                     return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Transform_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Transform(t,m);

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Transform_Term;
    procedure Transform_Terms is new Visiting_Iterator(Transform_Term);

  begin
    Transform_Terms(p);
    return res;
  end Transform;

  function Transform ( s : Laur_Sys; m : Matrix ) return Laur_Sys is

    res : Laur_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Transform(s(i),m);
    end loop;
    return res;
  end Transform;

  function Eliminate ( t : Standard_Complex_Laurentials.Term; k : integer32 )
                     return Standard_Complex_Laurentials.Term is

    use Standard_Complex_Laurentials;

    res : Term;

  begin
    res.cf := t.cf;
    res.dg := new Vector(t.dg'first..t.dg'last-1);
    for i in t.dg'first..k-1 loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in k+1..t.dg'last loop
      res.dg(i-1) := t.dg(i);
    end loop;
    return res;
  end Eliminate;

  function Eliminate ( p : Standard_Complex_Laurentials.Poly; k : integer32 )
                     return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Eliminate_Term ( t : in Term; continue : out boolean ) is

      nt : Term := Eliminate(t,k);

    begin
      Add(res,nt);
      Clear(nt);
      continue := true;
    end Eliminate_Term;
    procedure Eliminate_Terms is new Visiting_Iterator(Eliminate_Term);

  begin
    Eliminate_Terms(p);
    return res;
  end Eliminate;

  function Eliminate ( s : Laur_Sys; k : integer32 ) return Laur_Sys is

    res : Laur_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Eliminate(s(i),k);
    end loop;
    return res;
  end Eliminate;

end Standard_Initial_Forms;
