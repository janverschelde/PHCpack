with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Integer32_Vectors_Utilities;        use Integer32_Vectors_Utilities;

package body Transforming_Laurent_Systems is

  function Initial_Link_to_Vector
             ( p : Standard_Complex_Laurentials.Poly )
             return Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the initial degrees of the polynomial p.

    use Standard_Complex_Laurentials;
    init : Link_to_Vector;

    procedure Init_Term ( t : in Term; cont : out boolean ) is
    begin
      init := new Standard_Integer_Vectors.Vector'(t.dg.all);
      cont := false;
    end Init_Term;
    procedure Initial_Term is new Visiting_Iterator (Init_Term);

  begin
    Initial_Term(p);
    return init;
  end Initial_Link_to_Vector;

  function Initial_Link_to_Vector
             ( p : DoblDobl_Complex_Laurentials.Poly )
             return Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the initial degrees of the polynomial p.

    use DoblDobl_Complex_Laurentials;
    init : Link_to_Vector;

    procedure Init_Term ( t : in Term; cont : out boolean ) is
    begin
      init := new Standard_Integer_Vectors.Vector'(t.dg.all);
      cont := false;
    end Init_Term;
    procedure Initial_Term is new Visiting_Iterator (Init_Term);

  begin
    Initial_Term(p);
    return init;
  end Initial_Link_to_Vector;

  function Initial_Link_to_Vector
             ( p : QuadDobl_Complex_Laurentials.Poly )
             return Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the initial degrees of the polynomial p.

    use QuadDobl_Complex_Laurentials;
    init : Link_to_Vector;

    procedure Init_Term ( t : in Term; cont : out boolean ) is
    begin
      init := new Standard_Integer_Vectors.Vector'(t.dg.all);
      cont := false;
    end Init_Term;
    procedure Initial_Term is new Visiting_Iterator (Init_Term);

  begin
    Initial_Term(p);
    return init;
  end Initial_Link_to_Vector;

  procedure Shift ( p : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;
    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in out Term; cont : out boolean ) is
    begin
      Sub(Link_to_Vector(t.dg),init);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Changing_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
  end Shift;

  function Shift ( p : Standard_Complex_Laurentials.Poly )
                 return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;
    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := t.cf;
      rt.dg := t.dg - Degrees(init);
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Visiting_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
    return res;
  end Shift;

  procedure Shift ( p : in out DoblDobl_Complex_Laurentials.Poly ) is

    use DoblDobl_Complex_Laurentials;
    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in out Term; cont : out boolean ) is
    begin
      Sub(Link_to_Vector(t.dg),init);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Changing_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
  end Shift;

  function Shift ( p : DoblDobl_Complex_Laurentials.Poly )
                 return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Laurentials;
    res : Poly := Null_Poly;
    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := t.cf;
      rt.dg := t.dg - Degrees(init);
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Visiting_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
    return res;
  end Shift;

  procedure Shift ( p : in out QuadDobl_Complex_Laurentials.Poly ) is

    use QuadDobl_Complex_Laurentials;
    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in out Term; cont : out boolean ) is
    begin
      Sub(Link_to_Vector(t.dg),init);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Changing_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
  end Shift;

  function Shift ( p : QuadDobl_Complex_Laurentials.Poly )
                 return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;
    res : Poly := Null_Poly;
    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := t.cf;
      rt.dg := t.dg - Degrees(init);
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Visiting_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
    return res;
  end Shift;

  procedure Shift ( L : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for k in L'range loop
      Shift(l(k));
    end loop;
  end Shift;

  function Shift ( L : Standard_Complex_Laur_Systems.Laur_Sys )
                 return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for k in L'range loop
      res(k) := Shift(L(k));
    end loop;
    return res;
  end Shift;

  procedure Shift ( L : in out DoblDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for k in L'range loop
      Shift(l(k));
    end loop;
  end Shift;

  function Shift ( L : DoblDobl_Complex_Laur_Systems.Laur_Sys )
                 return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for k in L'range loop
      res(k) := Shift(L(k));
    end loop;
    return res;
  end Shift;

  procedure Shift ( L : in out QuadDobl_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for k in L'range loop
      Shift(l(k));
    end loop;
  end Shift;

  function Shift ( L : QuadDobl_Complex_Laur_Systems.Laur_Sys )
                 return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for k in L'range loop
      res(k) := Shift(L(k));
    end loop;
    return res;
  end Shift;

  procedure Transform ( t : in Transfo;
                        p : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    procedure Transform_Term ( tt : in out Term; cont : out boolean ) is
    begin
      Apply(t,Link_to_Vector(tt.dg));
      cont := true;
    end Transform_Term;
    procedure Transform_Terms is new Changing_Iterator (Transform_Term);

  begin
    Transform_Terms(p);
  end Transform;

  function Transform ( t : Transfo;
                       p : Standard_Complex_Laurentials.Poly )
                     return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly;

  begin
    Copy(p,res);
    Transform(t,res);
    return res;
  end Transform;

  function Transform2 ( t : Transfo;
                        p : Standard_Complex_Laurentials.Poly )
                      return Standard_Complex_Laurentials.Poly is

  -- IMPORTANT : This function might change the term order !

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Transform_Term ( tt : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := tt.cf;
      rt.dg := Degrees(t*Link_to_Vector(tt.dg));
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Transform_Term;
    procedure Transform_Terms is new Visiting_Iterator (Transform_Term);

  begin
    Transform_Terms(p);
    return res;
  end Transform2;

  procedure Transform ( t : in Transfo;
                        L : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in L'range loop
      Transform(t,L(i));
    end loop;
  end Transform;

  function Transform ( t : Transfo;
                       L : Standard_Complex_Laur_Systems.Laur_Sys )
                     return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in L'range loop
      res(i) := Transform(t,L(i));
    end loop;
    return res;
  end Transform;

  function Maximal_Support ( p : Standard_Complex_Laurentials.Poly;
                             v : Vector ) return integer32 is

    use Standard_Complex_Laurentials;
    res : integer32;
    first : boolean := true;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is

      sp : integer32 := t.dg.all*v;

    begin
      if first 
       then res := sp; first := false;
       elsif sp > res
           then res := sp;
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator (Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Maximal_Support;

  function Maximal_Support ( p : Standard_Complex_Laurentials.Poly;
                             v : Link_to_Vector ) return integer32 is
  begin
    return Maximal_Support(p,v.all);
  end Maximal_Support;

  procedure Face ( i,m : in integer32;
                   p : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    procedure Face_Term ( t : in out Term; cont : out boolean ) is
    begin
      if t.dg(i) /= m
       then t.cf := Create(0.0);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Changing_Iterator(Face_Term);

  begin
    Face_Terms(p);
  end Face;

  function Face ( i,m : integer32;
                  p : Standard_Complex_Laurentials.Poly )
                return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly;

  begin
    Copy(p,res);
    Face(i,m,res);
    return res;
  end Face;

  function Face2 ( i,m : integer32;
                   p : Standard_Complex_Laurentials.Poly )
                 return Standard_Complex_Laurentials.Poly is

  -- IMPORTANT : This function might change the term order !

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Face_Term ( t : in Term; cont : out boolean ) is
    begin
      if t.dg(i) = m
       then Add(res,t);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Visiting_Iterator(Face_Term);

  begin
    Face_Terms(p);
    return res;
  end Face2;

  procedure Face ( i,m : in integer32;
                   L : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for j in L'range loop
      Face(i,m,l(j));
    end loop;
  end Face;

  function Face ( i,m : integer32;
                  L : Standard_Complex_Laur_Systems.Laur_Sys )
                return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for j in L'range loop
      res(j) := Face(i,m,L(j));
    end loop;
    return res;
  end Face;

  procedure Face ( v : in Vector; m : in integer32;
                   p : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    procedure Face_Term ( t : in out Term; cont : out boolean ) is
    begin
      if t.dg.all*v /= m
       then t.cf := Create(0.0);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Changing_Iterator(Face_Term);

  begin
    Face_Terms(p);
  end Face;

  function Face ( v : Vector; m : integer32;
                  p : Standard_Complex_Laurentials.Poly )
                return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly;

  begin
    Copy(p,res);
    Face(v,m,res);
    return res;
  end Face;

  function Face2 ( v : Vector; m : integer32;
                   p : Standard_Complex_Laurentials.Poly )
                 return Standard_Complex_Laurentials.Poly is

  -- IMPORTANT : This procedure might change the term order !

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Face_Term ( t : in Term; cont : out boolean ) is
    begin
      if t.dg.all*v = m
       then Add(res,t);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Visiting_Iterator(Face_Term);

  begin
    Face_Terms(p);
    return res;
  end Face2;

  procedure Face ( v,m : in Vector;
                   L : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for i in L'range loop
      Face(v,m(i),L(i));
    end loop;
  end Face;

  function Face ( v,m : Vector;
                  L : Standard_Complex_Laur_Systems.Laur_Sys )
                return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for i in L'range loop
      res(i) := Face(v,m(i),L(i));
    end loop;
    return res;
  end Face;

  procedure Reduce ( i : in integer32;
                     p : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    procedure Reduce_Term ( t : in out Term; cont : out boolean ) is
    begin
      Reduce(Link_to_Vector(t.dg),i);
      cont := true;
    end Reduce_Term;
    procedure Reduce_Terms is new Changing_Iterator(Reduce_Term);

  begin
    Reduce_Terms(p);
  end Reduce;

  function Reduce ( i : integer32; p : Standard_Complex_Laurentials.Poly )
                  return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly;

  begin
    Copy(p,res);
    Reduce(i,res);
    return res;
  end Reduce;

  function Reduce2 ( i : integer32;
                     p : Standard_Complex_Laurentials.Poly )
                   return Standard_Complex_Laurentials.Poly is

  -- IMPORTANT : This function might change the term order !

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Reduce_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := t.cf;
      rt.dg := Degrees(Reduce(Link_to_Vector(t.dg),i));
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Reduce_Term;
    procedure Reduce_Terms is new Visiting_Iterator(Reduce_Term);

  begin
    Reduce_Terms(p);
    return res;
  end Reduce2;

  procedure Reduce ( i : in integer32;
                     L : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for j in L'range loop
      Reduce(i,L(j));
    end loop;
  end Reduce;

  function Reduce ( i : integer32;
                    L : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for j in L'range loop
      res(j) := Reduce(i,L(j));
    end loop;
    return res;
  end Reduce;

  procedure Insert ( i,d : in integer32;
                     p : in out Standard_Complex_Laurentials.Poly ) is

    use Standard_Complex_Laurentials;

    procedure Insert_Term ( t : in out Term; cont : out boolean ) is
    begin
      Insert(Link_to_Vector(t.dg),i,d);
      cont := true;
    end Insert_Term;
    procedure Insert_Terms is new Changing_Iterator(Insert_Term);

  begin
    Insert_Terms(p);
  end Insert;

  function Insert ( i,d : integer32;
                    p : Standard_Complex_Laurentials.Poly )
                  return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;
    res : Poly;

  begin
    Copy(p,res);
    Insert(i,d,res);
    return res;
  end Insert;

  function Insert2 ( i,d : integer32;
                     p : Standard_Complex_Laurentials.Poly )
                   return Standard_Complex_Laurentials.Poly is

  -- IMPORTANT : This function might change the term order !

    use Standard_Complex_Laurentials;
    res : Poly := Null_Poly;

    procedure Insert_Term ( t : in Term; cont : out boolean ) is

      rt : Term;

    begin
      rt.cf := t.cf;
      rt.dg := Degrees(Insert(Link_to_Vector(t.dg),i,d));
      Add(res,rt);
      Clear(rt);
      cont := true;
    end Insert_Term;
    procedure Insert_Terms is new Visiting_Iterator (Insert_Term);

  begin
    Insert_Terms(p);
    return res;
  end Insert2;

  procedure Insert ( i,d : in integer32;
                     L : in out Standard_Complex_Laur_Systems.Laur_Sys ) is
  begin
    for j in L'range loop
      Insert(i,d,l(j));
    end loop;
  end Insert;

  function Insert ( i,d : integer32;
                    L : Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(L'range);

  begin
    for j in L'range loop
      res(j) := Insert(i,d,l(j));
    end loop;
    return res;
  end Insert;

end Transforming_Laurent_Systems;
