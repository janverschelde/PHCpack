with Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;

package body Standard_to_Multprec_Convertors is

-- AUXILIARIES :

  function Convert ( d : Standard_Complex_Polynomials.Degrees )
                   return Multprec_Complex_Polynomials.Degrees is

    res : Multprec_Complex_Polynomials.Degrees;

  begin
    res := new Standard_Natural_Vectors.Vector(d'range);
    for i in res'range loop
      res(i) := d(i);
    end loop;
    return res;
  end Convert;

  function Convert ( d : Standard_Complex_Laurentials.Degrees )
                   return Multprec_Complex_Laurentials.Degrees is

    res : Multprec_Complex_Laurentials.Degrees;

  begin
    res := new Standard_Integer_Vectors.Vector(d'range);
    for i in res'range loop
      res(i) := d(i);
    end loop;
    return res;
  end Convert;

  function Convert ( t : Standard_Complex_Polynomials.Term ) 
                   return Multprec_Complex_Polynomials.Term is

    res : Multprec_Complex_Polynomials.Term;

  begin
    res.cf := Create(t.cf);
    res.dg := Convert(t.dg);
    return res;
  end Convert;

  function Convert ( t : Standard_Complex_Laurentials.Term ) 
                   return Multprec_Complex_Laurentials.Term is

    res : Multprec_Complex_Laurentials.Term;

  begin
    res.cf := Create(t.cf);
    res.dg := Convert(t.dg);
    return res;
  end Convert;

-- TARGET FUNCTIONS :

  function Convert ( p : Standard_Complex_Polynomials.Poly )
                   return Multprec_Complex_Polynomials.Poly is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    procedure Convert_Term ( t : in Standard_Complex_Polynomials.Term;
                             continue : out boolean ) is

      ct : Multprec_Complex_Polynomials.Term := Convert(t);

    begin
      Multprec_Complex_Polynomials.Add(res,ct);
      Multprec_Complex_Polynomials.Clear(ct);
      continue := true;
    end Convert_Term;
    procedure Convert_Terms is 
      new Standard_Complex_Polynomials.Visiting_Iterator(Convert_Term);

  begin
    Convert_Terms(p);
    return res;
  end Convert;

  function Convert ( p : Standard_Complex_Laurentials.Poly )
                   return Multprec_Complex_Laurentials.Poly is

    res : Multprec_Complex_Laurentials.Poly
        := Multprec_Complex_Laurentials.Null_Poly;

    procedure Convert_Term ( t : in Standard_Complex_Laurentials.Term;
                             continue : out boolean ) is

      ct : Multprec_Complex_Laurentials.Term := Convert(t);

    begin
      Multprec_Complex_Laurentials.Add(res,ct);
      Multprec_Complex_Laurentials.Clear(ct);
      continue := true;
    end Convert_Term;
    procedure Convert_Terms is 
      new Standard_Complex_Laurentials.Visiting_Iterator(Convert_Term);

  begin
    Convert_Terms(p);
    return res;
  end Convert;

  function Convert ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                   return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert(p(i));
    end loop;
    return res;
  end Convert;

  function Convert ( p : Standard_Complex_Laur_Systems.Laur_Sys )
                   return Multprec_Complex_Laur_Systems.Laur_Sys is

    res : Multprec_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Convert(p(i));
    end loop;
    return res;
  end Convert;

  procedure Set_Size ( p : in out Multprec_Complex_Polynomials.Poly;
                       size : in natural32 ) is

    res : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Polynomials.Null_Poly;

    procedure Set_Size_Term ( t : in Multprec_Complex_Polynomials.Term;
                              cont : out boolean ) is

      nt : Multprec_Complex_Polynomials.Term;

    begin
      nt.dg := t.dg;
      Multprec_Complex_Numbers.Copy(t.cf,nt.cf);
      Set_Size(nt.cf,size);
      Multprec_Complex_Polynomials.Add(res,nt);
      cont := true;
    end Set_Size_Term;
    procedure Set_Size_Terms is
      new Multprec_Complex_Polynomials.Visiting_Iterator(Set_Size_Term);
 
  begin
    Set_Size_Terms(p);
    Multprec_Complex_Polynomials.Clear(p);
    p := res;
  end Set_Size;

  procedure Set_Size ( p : in out Multprec_Complex_Laurentials.Poly;
                       size : in natural32 ) is

    res : Multprec_Complex_Laurentials.Poly
        := Multprec_Complex_Laurentials.Null_Poly;

    procedure Set_Size_Term ( t : in Multprec_Complex_Laurentials.Term;
                              cont : out boolean ) is

      nt : Multprec_Complex_Laurentials.Term;

    begin
      nt.dg := t.dg;
      Multprec_Complex_Numbers.Copy(t.cf,nt.cf);
      Set_Size(nt.cf,size);
      Multprec_Complex_Laurentials.Add(res,nt);
      cont := true;
    end Set_Size_Term;
    procedure Set_Size_Terms is
      new Multprec_Complex_Laurentials.Visiting_Iterator(Set_Size_Term);
 
  begin
    Set_Size_Terms(p);
    Multprec_Complex_Laurentials.Clear(p);
    p := res;
  end Set_Size;

  procedure Set_Size ( p : in out Multprec_Floating_Polynomials.Poly;
                       size : in natural32 ) is

    res : Multprec_Floating_Polynomials.Poly
        := Multprec_Floating_Polynomials.Null_Poly;

    procedure Set_Size_Term ( t : in Multprec_Floating_Polynomials.Term;
                              cont : out boolean ) is

      nt : Multprec_Floating_Polynomials.Term;

    begin
      nt.dg := t.dg;
      Multprec_Floating_Numbers.Copy(t.cf,nt.cf);
      Multprec_Floating_Numbers.Set_Size(nt.cf,size);
      Multprec_Floating_Polynomials.Add(res,nt);
      cont := true;
    end Set_Size_Term;
    procedure Set_Size_Terms is
      new Multprec_Floating_Polynomials.Visiting_Iterator(Set_Size_Term);
 
  begin
    Set_Size_Terms(p);
    Multprec_Floating_Polynomials.Clear(p);
    p := res;
  end Set_Size;

  procedure Set_Size ( p : in out Multprec_Complex_Poly_Systems.Poly_Sys;
                       size : in natural32 ) is
  begin
    for i in p'range loop
      Set_Size(p(i),size);
    end loop;
  end Set_Size;

  procedure Set_Size ( p : in out Multprec_Complex_Laur_Systems.Laur_Sys;
                       size : in natural32 ) is
  begin
    for i in p'range loop
      Set_Size(p(i),size);
    end loop;
  end Set_Size;

  procedure Set_Size ( p : in out Multprec_Floating_Poly_Systems.Poly_Sys;
                       size : in natural32 ) is
  begin
    for i in p'range loop
      Set_Size(p(i),size);
    end loop;
  end Set_Size;

end Standard_to_Multprec_Convertors;
