with unchecked_deallocation;

package body Exponent_Vectors is

-- CREATORS :

  function Create ( p : Standard_Complex_Laurentials.Poly )
                  return Standard_Integer_VecVecs.VecVec is

    use Standard_Complex_Laurentials;
    res : Standard_Integer_VecVecs.VecVec(1..integer32(Number_of_Terms(p)));
    ind : integer32 := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Standard_Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : Standard_Complex_Polynomials.Poly )
                  return Standard_Integer_VecVecs.VecVec is

    use Standard_Complex_Polynomials;
    res : Standard_Integer_VecVecs.VecVec(1..integer32(Number_of_Terms(p)));
    ind : integer32 := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Standard_Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := integer32(t.dg(i));
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Laurentials.Poly )
                  return Standard_Integer_VecVecs.VecVec is

    use DoblDobl_Complex_Laurentials;
    res : Standard_Integer_VecVecs.VecVec(1..integer32(Number_of_Terms(p)));
    ind : integer32 := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Standard_Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Polynomials.Poly )
                  return Standard_Integer_VecVecs.VecVec is

    use DoblDobl_Complex_Polynomials;
    res : Standard_Integer_VecVecs.VecVec(1..integer32(Number_of_Terms(p)));
    ind : integer32 := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Standard_Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := integer32(t.dg(i));
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Laurentials.Poly )
                  return Standard_Integer_VecVecs.VecVec is

    use QuadDobl_Complex_Laurentials;
    res : Standard_Integer_VecVecs.VecVec(1..integer32(Number_of_Terms(p)));
    ind : integer32 := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Standard_Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Polynomials.Poly )
                  return Standard_Integer_VecVecs.VecVec is

    use QuadDobl_Complex_Polynomials;
    res : Standard_Integer_VecVecs.VecVec(1..integer32(Number_of_Terms(p)));
    ind : integer32 := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Standard_Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := integer32(t.dg(i));
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                  return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Standard_Integer_VecVecs.VecVec := Create(p(i));
      begin
        res(i) := new Standard_Integer_VecVecs.VecVec(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new VecVec'(Create(p(i)));
    end loop;
    return res;
  end Create;

  function Create ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Standard_Integer_VecVecs.VecVec := Create(p(i));
      begin
        res(i) := new Standard_Integer_VecVecs.VecVec(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new VecVec'(Create(p(i)));
    end loop;
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                  return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Standard_Integer_VecVecs.VecVec := Create(p(i));
      begin
        res(i) := new Standard_Integer_VecVecs.VecVec(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new VecVec'(Create(p(i)));
    end loop;
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Standard_Integer_VecVecs.VecVec := Create(p(i));
      begin
        res(i) := new Standard_Integer_VecVecs.VecVec(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new VecVec'(Create(p(i)));
    end loop;
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                  return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Standard_Integer_VecVecs.VecVec := Create(p(i));
      begin
        res(i) := new Standard_Integer_VecVecs.VecVec(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new VecVec'(Create(p(i)));
    end loop;
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Standard_Integer_VecVecs.VecVec := Create(p(i));
      begin
        res(i) := new Standard_Integer_VecVecs.VecVec(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new VecVec'(Create(p(i)));
    end loop;
    return res;
  end Create;

-- SELECTOR :

  function Position ( ev : Standard_Integer_VecVecs.VecVec;
                      v : Standard_Integer_Vectors.Vector )
                    return integer32 is
  begin
    for i in ev'range loop
      if Standard_Integer_Vectors.Equal(ev(i).all,v)
       then return i;
      end if;
    end loop;
    return ev'last+1;
  end Position;

-- EVALUATORS :

  function Eval ( e : Standard_Integer_Vectors.Vector;
                  c : Standard_Complex_Numbers.Complex_Number;
                  x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;
    res : Complex_Number := c;

  begin
    for i in e'range loop
      for j in 1..e(i) loop
        res := res*x(i);
      end loop;
      for j in 1..-e(i) loop
        res := res/x(i);
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( ev : Standard_Integer_VecVecs.VecVec;
                  c,x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;
    res : Complex_Number := Eval(ev(ev'first).all,c(c'first),x);

  begin
    for i in c'first+1..c'last loop
      res := res + Eval(ev(i).all,c(i),x);
    end loop;
    return res;
  end Eval;

  function Eval ( ev : Exponent_Vectors_Array;
                  c : Standard_Complex_VecVecs.VecVec;
                  x : Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    for i in ev'range loop
      res(i) := Eval(ev(i).all,c(i).all,x);
    end loop;
    return res;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( v : in out Exponent_Vectors_Array ) is
  begin
    for i in v'range loop
      Standard_Integer_VecVecs.Deep_Clear(v(i));
    end loop;
  end Clear;

  procedure Clear ( v : in out Link_to_Exponent_Vectors_Array ) is

    procedure free is
      new unchecked_deallocation(Exponent_Vectors_Array,
                                 Link_to_Exponent_Vectors_Array);

  begin
    if v /= null
     then Clear(v.all); free(v);
    end if;
  end Clear;

end Exponent_Vectors;
