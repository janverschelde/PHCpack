with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Random_Numbers;
with TripDobl_Complex_Numbers;
with TripDobl_Random_Numbers;
with QuadDobl_Complex_Numbers;
with QuadDobl_Random_Numbers;
with PentDobl_Complex_Numbers;
with PentDobl_Random_Numbers;
with OctoDobl_Complex_Numbers;
with OctoDobl_Random_Numbers;
with DecaDobl_Complex_Numbers;
with DecaDobl_Random_Numbers;
with HexaDobl_Complex_Numbers;
with HexaDobl_Random_Numbers;
with Degrees_in_Sets_of_Unknowns;

package body Multi_Projective_Transformations is

  function Multiset_Degrees
             ( p : in Standard_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in DoblDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in TripDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in QuadDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in PentDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in OctoDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in DecaDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Multiset_Degrees
             ( p : in HexaDobl_Complex_Polynomials.Poly;
               m : in natural32; z : in Partition )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in z'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,z(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Make_Homogeneous
             ( t : Standard_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : DoblDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return DoblDobl_Complex_Polynomials.Term is

    res : DoblDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : TripDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return TripDobl_Complex_Polynomials.Term is

    res : TripDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : QuadDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : PentDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return PentDobl_Complex_Polynomials.Term is

    res : PentDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : OctoDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return OctoDobl_Complex_Polynomials.Term is

    res : OctoDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : DecaDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return DecaDobl_Complex_Polynomials.Term is

    res : DecaDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( t : HexaDobl_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; z : Partition )
             return HexaDobl_Complex_Polynomials.Term is

    res : HexaDobl_Complex_Polynomials.Term;
    itm : constant integer32 := integer32(m);
    lst : constant integer32 := t.dg'last;
    deg : integer32;
 
  begin
    res.cf := t.cf;
    res.dg := new Standard_Natural_Vectors.Vector(1..lst+itm);
    for i in t.dg'range loop
      res.dg(i) := t.dg(i);
    end loop;
    for i in 1..itm loop
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,z(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in Standard_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return Standard_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in Standard_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : Standard_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      Standard_Complex_Polynomials.Add(res,rt);
      Standard_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      Standard_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in DoblDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return DoblDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DoblDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : DoblDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      DoblDobl_Complex_Polynomials.Add(res,rt);
      DoblDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      DoblDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in TripDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return TripDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : TripDobl_Complex_Polynomials.Poly
        := TripDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in TripDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : TripDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      TripDobl_Complex_Polynomials.Add(res,rt);
      TripDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      TripDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in QuadDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return QuadDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in QuadDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : QuadDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      QuadDobl_Complex_Polynomials.Add(res,rt);
      QuadDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      QuadDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in PentDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return PentDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : PentDobl_Complex_Polynomials.Poly
        := PentDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in PentDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : PentDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      PentDobl_Complex_Polynomials.Add(res,rt);
      PentDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      PentDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in OctoDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return OctoDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : OctoDobl_Complex_Polynomials.Poly
        := OctoDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in OctoDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : OctoDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      OctoDobl_Complex_Polynomials.Add(res,rt);
      OctoDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      OctoDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in DecaDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return DecaDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : DecaDobl_Complex_Polynomials.Poly
        := DecaDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in DecaDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : DecaDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      DecaDobl_Complex_Polynomials.Add(res,rt);
      DecaDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      DecaDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in HexaDobl_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return HexaDobl_Complex_Polynomials.Poly is

    deg : constant Standard_Integer_Vectors.Vector(1..integer32(m))
        := Multiset_Degrees(p,m,z);
    res : HexaDobl_Complex_Polynomials.Poly
        := HexaDobl_Complex_Polynomials.Null_Poly;

    procedure Visit_Term ( t : in HexaDobl_Complex_Polynomials.Term;
                           continue : out boolean ) is

      rt : HexaDobl_Complex_Polynomials.Term
         := Make_Homogeneous(t,deg,m,z);

    begin
      HexaDobl_Complex_Polynomials.Add(res,rt);
      HexaDobl_Complex_Polynomials.Clear(rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new
      HexaDobl_Complex_Polynomials.Visiting_Iterator(Visit_Term);
 
  begin
    Visit_Terms(p);
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in TripDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in PentDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    res : PentDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in OctoDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in DecaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in HexaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Standard_Random_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;

  begin
    res.cf := Standard_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end Standard_Random_Linear_Term;

  function DoblDobl_Random_Linear_Term
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Term is

    res : DoblDobl_Complex_Polynomials.Term;

  begin
    res.cf := DoblDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end DoblDobl_Random_Linear_Term;

  function TripDobl_Random_Linear_Term
             ( n,i : natural32 )
             return TripDobl_Complex_Polynomials.Term is

    res : TripDobl_Complex_Polynomials.Term;

  begin
    res.cf := TripDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end TripDobl_Random_Linear_Term;

  function QuadDobl_Random_Linear_Term
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;

  begin
    res.cf := QuadDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end QuadDobl_Random_Linear_Term;

  function PentDobl_Random_Linear_Term
             ( n,i : natural32 )
             return PentDobl_Complex_Polynomials.Term is

    res : PentDobl_Complex_Polynomials.Term;

  begin
    res.cf := PentDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end PentDobl_Random_Linear_Term;

  function OctoDobl_Random_Linear_Term
             ( n,i : natural32 )
             return OctoDobl_Complex_Polynomials.Term is

    res : OctoDobl_Complex_Polynomials.Term;

  begin
    res.cf := OctoDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end OctoDobl_Random_Linear_Term;

  function DecaDobl_Random_Linear_Term
             ( n,i : natural32 )
             return DecaDobl_Complex_Polynomials.Term is

    res : DecaDobl_Complex_Polynomials.Term;

  begin
    res.cf := DecaDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end DecaDobl_Random_Linear_Term;

  function HexaDobl_Random_Linear_Term
             ( n,i : natural32 )
             return HexaDobl_Complex_Polynomials.Term is

    res : HexaDobl_Complex_Polynomials.Term;

  begin
    res.cf := HexaDobl_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end HexaDobl_Random_Linear_Term;

  function Standard_Start_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term is

    res : Standard_Complex_Polynomials.Term;

  begin
    res.cf := Standard_Complex_Numbers.Create(1.0);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end Standard_Start_Linear_Term;

  function DoblDobl_Start_Linear_Term
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Term is

    res : DoblDobl_Complex_Polynomials.Term;
    one : constant double_double := create(1.0);

  begin
    res.cf := DoblDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end DoblDobl_Start_Linear_Term;

  function TripDobl_Start_Linear_Term
             ( n,i : natural32 )
             return TripDobl_Complex_Polynomials.Term is

    res : TripDobl_Complex_Polynomials.Term;
    one : constant triple_double := create(1.0);

  begin
    res.cf := TripDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end TripDobl_Start_Linear_Term;

  function QuadDobl_Start_Linear_Term
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Term is

    res : QuadDobl_Complex_Polynomials.Term;
    one : constant quad_double := create(1.0);

  begin
    res.cf := QuadDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end QuadDobl_Start_Linear_Term;

  function PentDobl_Start_Linear_Term
             ( n,i : natural32 )
             return PentDobl_Complex_Polynomials.Term is

    res : PentDobl_Complex_Polynomials.Term;
    one : constant penta_double := create(1.0);

  begin
    res.cf := PentDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end PentDobl_Start_Linear_Term;

  function OctoDobl_Start_Linear_Term
             ( n,i : natural32 )
             return OctoDobl_Complex_Polynomials.Term is

    res : OctoDobl_Complex_Polynomials.Term;
    one : constant octo_double := create(1.0);

  begin
    res.cf := OctoDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end OctoDobl_Start_Linear_Term;

  function DecaDobl_Start_Linear_Term
             ( n,i : natural32 )
             return DecaDobl_Complex_Polynomials.Term is

    res : DecaDobl_Complex_Polynomials.Term;
    one : constant deca_double := create(1.0);

  begin
    res.cf := DecaDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end DecaDobl_Start_Linear_Term;

  function HexaDobl_Start_Linear_Term
             ( n,i : natural32 )
             return HexaDobl_Complex_Polynomials.Term is

    res : HexaDobl_Complex_Polynomials.Term;
    one : constant hexa_double := create(1.0);

  begin
    res.cf := HexaDobl_Complex_Numbers.Create(one);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end HexaDobl_Start_Linear_Term;

  function Standard_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return Standard_Complex_Polynomials.Poly is

    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : Standard_Complex_Polynomials.Term
            := Standard_Random_Linear_Term(n,i);
        begin
          Standard_Complex_Polynomials.Add(res,t);
          Standard_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end Standard_Random_Linear_Polynomial;

  function DoblDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return DoblDobl_Complex_Polynomials.Poly is

    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : DoblDobl_Complex_Polynomials.Term
            := DoblDobl_Random_Linear_Term(n,i);
        begin
          DoblDobl_Complex_Polynomials.Add(res,t);
          DoblDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end DoblDobl_Random_Linear_Polynomial;

  function TripDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return TripDobl_Complex_Polynomials.Poly is

    res : TripDobl_Complex_Polynomials.Poly
        := TripDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : TripDobl_Complex_Polynomials.Term
            := TripDobl_Random_Linear_Term(n,i);
        begin
          TripDobl_Complex_Polynomials.Add(res,t);
          TripDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end TripDobl_Random_Linear_Polynomial;

  function QuadDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return QuadDobl_Complex_Polynomials.Poly is

    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : QuadDobl_Complex_Polynomials.Term
            := QuadDobl_Random_Linear_Term(n,i);
        begin
          QuadDobl_Complex_Polynomials.Add(res,t);
          QuadDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end QuadDobl_Random_Linear_Polynomial;

  function PentDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return PentDobl_Complex_Polynomials.Poly is

    res : PentDobl_Complex_Polynomials.Poly
        := PentDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : PentDobl_Complex_Polynomials.Term
            := PentDobl_Random_Linear_Term(n,i);
        begin
          PentDobl_Complex_Polynomials.Add(res,t);
          PentDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end PentDobl_Random_Linear_Polynomial;

  function OctoDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return OctoDobl_Complex_Polynomials.Poly is

    res : OctoDobl_Complex_Polynomials.Poly
        := OctoDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : OctoDobl_Complex_Polynomials.Term
            := OctoDobl_Random_Linear_Term(n,i);
        begin
          OctoDobl_Complex_Polynomials.Add(res,t);
          OctoDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end OctoDobl_Random_Linear_Polynomial;

  function DecaDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return DecaDobl_Complex_Polynomials.Poly is

    res : DecaDobl_Complex_Polynomials.Poly
        := DecaDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : DecaDobl_Complex_Polynomials.Term
            := DecaDobl_Random_Linear_Term(n,i);
        begin
          DecaDobl_Complex_Polynomials.Add(res,t);
          DecaDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end DecaDobl_Random_Linear_Polynomial;

  function HexaDobl_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return HexaDobl_Complex_Polynomials.Poly is

    res : HexaDobl_Complex_Polynomials.Poly
        := HexaDobl_Complex_Polynomials.Null_Poly;

  begin
    for i in 1..Sets_of_Unknowns.Dimension(s) loop
      if Sets_of_Unknowns.Is_In(s,i) then
        declare
          t : HexaDobl_Complex_Polynomials.Term
            := HexaDobl_Random_Linear_Term(n,i);
        begin
          HexaDobl_Complex_Polynomials.Add(res,t);
          HexaDobl_Complex_Polynomials.Clear(t);
        end;
      end if;
    end loop;
    return res;
  end HexaDobl_Random_Linear_Polynomial;

  function Standard_Start_Linear_Polynomial
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Poly is

    trm : Standard_Complex_Polynomials.Term
        := Standard_Start_Linear_Term(n,i);
    res : Standard_Complex_Polynomials.Poly
        := Standard_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    Standard_Complex_Polynomials.Sub(res,trm);
    Standard_Complex_Polynomials.Clear(trm);
    return res;
  end Standard_Start_Linear_Polynomial;

  function DoblDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Poly is

    trm : DoblDobl_Complex_Polynomials.Term
        := DoblDobl_Start_Linear_Term(n,i);
    res : DoblDobl_Complex_Polynomials.Poly
        := DoblDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    DoblDobl_Complex_Polynomials.Sub(res,trm);
    DoblDobl_Complex_Polynomials.Clear(trm);
    return res;
  end DoblDobl_Start_Linear_Polynomial;

  function TripDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return TripDobl_Complex_Polynomials.Poly is

    trm : TripDobl_Complex_Polynomials.Term
        := TripDobl_Start_Linear_Term(n,i);
    res : TripDobl_Complex_Polynomials.Poly
        := TripDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    TripDobl_Complex_Polynomials.Sub(res,trm);
    TripDobl_Complex_Polynomials.Clear(trm);
    return res;
  end TripDobl_Start_Linear_Polynomial;

  function QuadDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Poly is

    trm : QuadDobl_Complex_Polynomials.Term
        := QuadDobl_Start_Linear_Term(n,i);
    res : QuadDobl_Complex_Polynomials.Poly
        := QuadDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    QuadDobl_Complex_Polynomials.Sub(res,trm);
    QuadDobl_Complex_Polynomials.Clear(trm);
    return res;
  end QuadDobl_Start_Linear_Polynomial;

  function PentDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return PentDobl_Complex_Polynomials.Poly is

    trm : PentDobl_Complex_Polynomials.Term
        := PentDobl_Start_Linear_Term(n,i);
    res : PentDobl_Complex_Polynomials.Poly
        := PentDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    PentDobl_Complex_Polynomials.Sub(res,trm);
    PentDobl_Complex_Polynomials.Clear(trm);
    return res;
  end PentDobl_Start_Linear_Polynomial;

  function OctoDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return OctoDobl_Complex_Polynomials.Poly is

    trm : OctoDobl_Complex_Polynomials.Term
        := OctoDobl_Start_Linear_Term(n,i);
    res : OctoDobl_Complex_Polynomials.Poly
        := OctoDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    OctoDobl_Complex_Polynomials.Sub(res,trm);
    OctoDobl_Complex_Polynomials.Clear(trm);
    return res;
  end OctoDobl_Start_Linear_Polynomial;

  function DecaDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return DecaDobl_Complex_Polynomials.Poly is

    trm : DecaDobl_Complex_Polynomials.Term
        := DecaDobl_Start_Linear_Term(n,i);
    res : DecaDobl_Complex_Polynomials.Poly
        := DecaDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    DecaDobl_Complex_Polynomials.Sub(res,trm);
    DecaDobl_Complex_Polynomials.Clear(trm);
    return res;
  end DecaDobl_Start_Linear_Polynomial;

  function HexaDobl_Start_Linear_Polynomial
             ( n,i : natural32 )
             return HexaDobl_Complex_Polynomials.Poly is

    trm : HexaDobl_Complex_Polynomials.Term
        := HexaDobl_Start_Linear_Term(n,i);
    res : HexaDobl_Complex_Polynomials.Poly
        := HexaDobl_Complex_Polynomials.Create(trm);

  begin
    trm.dg(integer32(i)) := 0;
    HexaDobl_Complex_Polynomials.Sub(res,trm);
    HexaDobl_Complex_Polynomials.Clear(trm);
    return res;
  end HexaDobl_Start_Linear_Polynomial;

  function Standard_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : Standard_Complex_Polynomials.Term;
    ztm : Standard_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := Standard_Random_Linear_Polynomial(dim,z(i));
      cst.cf := Standard_Random_Numbers.Random1;
      ztm.cf := Standard_Random_Numbers.Random1;
      Standard_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      Standard_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    Standard_Complex_Polynomials.Clear(cst);
    Standard_Complex_Polynomials.Clear(ztm);
    return res;
  end Standard_Random_Linear_Polynomials;

  function DoblDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : DoblDobl_Complex_Polynomials.Term;
    ztm : DoblDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := DoblDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := DoblDobl_Random_Numbers.Random1;
      ztm.cf := DoblDobl_Random_Numbers.Random1;
      DoblDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      DoblDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    DoblDobl_Complex_Polynomials.Clear(cst);
    DoblDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end DoblDobl_Random_Linear_Polynomials;

  function TripDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : TripDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : TripDobl_Complex_Polynomials.Term;
    ztm : TripDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := TripDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := TripDobl_Random_Numbers.Random1;
      ztm.cf := TripDobl_Random_Numbers.Random1;
      TripDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      TripDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    TripDobl_Complex_Polynomials.Clear(cst);
    TripDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end TripDobl_Random_Linear_Polynomials;

  function QuadDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : QuadDobl_Complex_Polynomials.Term;
    ztm : QuadDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := QuadDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := QuadDobl_Random_Numbers.Random1;
      ztm.cf := QuadDobl_Random_Numbers.Random1;
      QuadDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      QuadDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    QuadDobl_Complex_Polynomials.Clear(cst);
    QuadDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end QuadDobl_Random_Linear_Polynomials;

  function PentDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : PentDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : PentDobl_Complex_Polynomials.Term;
    ztm : PentDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := PentDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := PentDobl_Random_Numbers.Random1;
      ztm.cf := PentDobl_Random_Numbers.Random1;
      PentDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      PentDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    PentDobl_Complex_Polynomials.Clear(cst);
    PentDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end PentDobl_Random_Linear_Polynomials;

  function OctoDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : OctoDobl_Complex_Polynomials.Term;
    ztm : OctoDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := OctoDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := OctoDobl_Random_Numbers.Random1;
      ztm.cf := OctoDobl_Random_Numbers.Random1;
      OctoDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      OctoDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    OctoDobl_Complex_Polynomials.Clear(cst);
    OctoDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end OctoDobl_Random_Linear_Polynomials;

  function DecaDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : DecaDobl_Complex_Polynomials.Term;
    ztm : DecaDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := DecaDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := DecaDobl_Random_Numbers.Random1;
      ztm.cf := DecaDobl_Random_Numbers.Random1;
      DecaDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      DecaDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    DecaDobl_Complex_Polynomials.Clear(cst);
    DecaDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end DecaDobl_Random_Linear_Polynomials;

  function HexaDobl_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    cst : HexaDobl_Complex_Polynomials.Term;
    ztm : HexaDobl_Complex_Polynomials.Term;

  begin
    cst.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    ztm.dg := new Standard_Natural_Vectors.Vector'(1..integer32(dim) => 0);
    for i in 1..m loop
      res(integer32(i)) := HexaDobl_Random_Linear_Polynomial(dim,z(i));
      cst.cf := HexaDobl_Random_Numbers.Random1;
      ztm.cf := HexaDobl_Random_Numbers.Random1;
      HexaDobl_Complex_Polynomials.Add(res(integer32(i)),cst);
      ztm.dg(integer32(n+i)) := 1;
      HexaDobl_Complex_Polynomials.Add(res(integer32(i)),ztm);
      ztm.dg(integer32(n+i)) := 0;
    end loop;
    HexaDobl_Complex_Polynomials.Clear(cst);
    HexaDobl_Complex_Polynomials.Clear(ztm);
    return res;
  end HexaDobl_Random_Linear_Polynomials;

  function Standard_Start_Linear_Polynomials
             ( n,m : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := Standard_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end Standard_Start_Linear_Polynomials;

  function DoblDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := DoblDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end DoblDobl_Start_Linear_Polynomials;

  function TripDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : TripDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := TripDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end TripDobl_Start_Linear_Polynomials;

  function QuadDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := QuadDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end QuadDobl_Start_Linear_Polynomials;

  function PentDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : PentDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := PentDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end PentDobl_Start_Linear_Polynomials;

  function OctoDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := OctoDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end OctoDobl_Start_Linear_Polynomials;

  function DecaDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := DecaDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end DecaDobl_Start_Linear_Polynomials;

  function HexaDobl_Start_Linear_Polynomials
             ( n,m : natural32 )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant natural32 := n+m;
    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := HexaDobl_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end HexaDobl_Start_Linear_Polynomials;

  function Add_Ones ( s : Standard_Complex_Solutions.Solution;
                      m : natural32 )
                    return Standard_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : Standard_Complex_Solutions.Solution(dim+integer32(m));

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := Standard_Complex_Numbers.Create(1.0);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : DoblDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return DoblDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : DoblDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant double_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := DoblDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : TripDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return TripDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : TripDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant triple_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := TripDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : QuadDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return QuadDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : QuadDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant quad_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := QuadDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : PentDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return PentDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : PentDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant penta_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := PentDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : OctoDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return OctoDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : OctoDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant octo_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := OctoDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : DecaDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return DecaDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : DecaDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant deca_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := DecaDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( s : HexaDobl_Complex_Solutions.Solution;
                      m : natural32 )
                    return HexaDobl_Complex_Solutions.Solution is

    dim : constant integer32 := s.n;
    res : HexaDobl_Complex_Solutions.Solution(dim+integer32(m));
    one : constant hexa_double := create(1.0);

  begin
    res.v(1..dim) := s.v(1..dim);
    for k in 1..integer32(m) loop
      res.v(dim+k) := HexaDobl_Complex_Numbers.Create(one);
    end loop;
    res.t := s.t;
    res.m := s.m;
    res.err := s.err;
    res.rco := s.rco;
    res.res := s.res;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : Standard_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return Standard_Complex_Solutions.Solution_List is

    res,res_last : Standard_Complex_Solutions.Solution_List;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : DoblDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return DoblDobl_Complex_Solutions.Solution_List is

    res,res_last : DoblDobl_Complex_Solutions.Solution_List;
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      DoblDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : TripDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return TripDobl_Complex_Solutions.Solution_List is

    res,res_last : TripDobl_Complex_Solutions.Solution_List;
    tmp : TripDobl_Complex_Solutions.Solution_List := sols;
    ls : TripDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not TripDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := TripDobl_Complex_Solutions.Head_Of(tmp);
      TripDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := TripDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : QuadDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return QuadDobl_Complex_Solutions.Solution_List is

    res,res_last : QuadDobl_Complex_Solutions.Solution_List;
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : PentDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return PentDobl_Complex_Solutions.Solution_List is

    res,res_last : PentDobl_Complex_Solutions.Solution_List;
    tmp : PentDobl_Complex_Solutions.Solution_List := sols;
    ls : PentDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not PentDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := PentDobl_Complex_Solutions.Head_Of(tmp);
      PentDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := PentDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : OctoDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return OctoDobl_Complex_Solutions.Solution_List is

    res,res_last : OctoDobl_Complex_Solutions.Solution_List;
    tmp : OctoDobl_Complex_Solutions.Solution_List := sols;
    ls : OctoDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not OctoDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := OctoDobl_Complex_Solutions.Head_Of(tmp);
      OctoDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := OctoDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : DecaDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return DecaDobl_Complex_Solutions.Solution_List is

    res,res_last : DecaDobl_Complex_Solutions.Solution_List;
    tmp : DecaDobl_Complex_Solutions.Solution_List := sols;
    ls : DecaDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DecaDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DecaDobl_Complex_Solutions.Head_Of(tmp);
      DecaDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := DecaDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  function Add_Ones ( sols : HexaDobl_Complex_Solutions.Solution_List;
                      m : natural32 )
                    return HexaDobl_Complex_Solutions.Solution_List is

    res,res_last : HexaDobl_Complex_Solutions.Solution_List;
    tmp : HexaDobl_Complex_Solutions.Solution_List := sols;
    ls : HexaDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not HexaDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := HexaDobl_Complex_Solutions.Head_Of(tmp);
      HexaDobl_Complex_Solutions.Append(res,res_last,Add_Ones(ls.all,m));
      tmp := HexaDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Add_Ones;

  procedure Add_Ones ( sols : in out Standard_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant Standard_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        Standard_Complex_Solutions.Clear(ls);
        ls := new Standard_Complex_Solutions.Solution'(sol);
        Standard_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant DoblDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        DoblDobl_Complex_Solutions.Clear(ls);
        ls := new DoblDobl_Complex_Solutions.Solution'(sol);
        DoblDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out TripDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : TripDobl_Complex_Solutions.Solution_List := sols;
    ls : TripDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not TripDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := TripDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant TripDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        TripDobl_Complex_Solutions.Clear(ls);
        ls := new TripDobl_Complex_Solutions.Solution'(sol);
        TripDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := TripDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant QuadDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        QuadDobl_Complex_Solutions.Clear(ls);
        ls := new QuadDobl_Complex_Solutions.Solution'(sol);
        QuadDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out PentDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : PentDobl_Complex_Solutions.Solution_List := sols;
    ls : PentDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not PentDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := PentDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant PentDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        PentDobl_Complex_Solutions.Clear(ls);
        ls := new PentDobl_Complex_Solutions.Solution'(sol);
        PentDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := PentDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out OctoDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : OctoDobl_Complex_Solutions.Solution_List := sols;
    ls : OctoDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not OctoDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := OctoDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant OctoDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        OctoDobl_Complex_Solutions.Clear(ls);
        ls := new OctoDobl_Complex_Solutions.Solution'(sol);
        OctoDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := OctoDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out DecaDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : DecaDobl_Complex_Solutions.Solution_List := sols;
    ls : DecaDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DecaDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DecaDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant DecaDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        DecaDobl_Complex_Solutions.Clear(ls);
        ls := new DecaDobl_Complex_Solutions.Solution'(sol);
        DecaDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := DecaDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  procedure Add_Ones ( sols : in out HexaDobl_Complex_Solutions.Solution_List;
                       m : in natural32 ) is

    tmp : HexaDobl_Complex_Solutions.Solution_List := sols;
    ls : HexaDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not HexaDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := HexaDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant HexaDobl_Complex_Solutions.Solution(dim+integer32(m))
            := Add_Ones(ls.all,m);
      begin
        HexaDobl_Complex_Solutions.Clear(ls);
        ls := new HexaDobl_Complex_Solutions.Solution'(sol);
        HexaDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := HexaDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Add_Ones;

  function Make_Affine ( sol : Standard_Complex_Solutions.Solution;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return Standard_Complex_Solutions.Solution is

    res : Standard_Complex_Solutions.Solution(sol.n-integer32(m));
    idx : integer32;

    use Standard_Complex_Numbers;

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    for k in 1..m loop
      idx := idz'last + integer32(k); -- homogeneous variable for k-th set
      for i in idz'range loop
        if idz(i) = k then
          res.v(i) := sol.v(i)/sol.v(idx);
        end if;
      end loop;
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine ( sol : DoblDobl_Complex_Solutions.Solution;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return DoblDobl_Complex_Solutions.Solution is

    res : DoblDobl_Complex_Solutions.Solution(sol.n-integer32(m));
    idx : integer32;

    use DoblDobl_Complex_Numbers;

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    for k in 1..m loop
      idx := idz'last + integer32(k); -- homogeneous variable for k-th set
      for i in idz'range loop
        if idz(i) = k then
          res.v(i) := sol.v(i)/sol.v(idx);
        end if;
      end loop;
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine ( sol : QuadDobl_Complex_Solutions.Solution;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return QuadDobl_Complex_Solutions.Solution is

    res : QuadDobl_Complex_Solutions.Solution(sol.n-integer32(m));
    idx : integer32;

    use QuadDobl_Complex_Numbers;

  begin
    res.t := sol.t;
    res.m := sol.m;
    res.err := sol.err;
    res.rco := sol.rco;
    res.res := sol.res;
    for k in 1..m loop
      idx := idz'last + integer32(k); -- homogeneous variable for k-th set
      for i in idz'range loop
        if idz(i) = k then
          res.v(i) := sol.v(i)/sol.v(idx);
        end if;
      end loop;
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine ( sols : Standard_Complex_Solutions.Solution_List;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return Standard_Complex_Solutions.Solution_List is

    res,res_last : Standard_Complex_Solutions.Solution_List;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Complex_Solutions.Append(res,res_last,Make_Affine(ls.all,m,idz));
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine ( sols : DoblDobl_Complex_Solutions.Solution_List;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return DoblDobl_Complex_Solutions.Solution_List is

    res,res_last : DoblDobl_Complex_Solutions.Solution_List;
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      DoblDobl_Complex_Solutions.Append(res,res_last,Make_Affine(ls.all,m,idz));
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Make_Affine;

  function Make_Affine ( sols : QuadDobl_Complex_Solutions.Solution_List;
                         m : natural32;
                         idz : Standard_Natural_Vectors.Vector )
                       return QuadDobl_Complex_Solutions.Solution_List is

    res,res_last : QuadDobl_Complex_Solutions.Solution_List;
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Complex_Solutions.Append(res,res_last,Make_Affine(ls.all,m,idz));
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    return res;
  end Make_Affine;

  procedure Make_Affine
              ( sols : in out Standard_Complex_Solutions.Solution_List;
                m : in natural32;
                idz : in Standard_Natural_Vectors.Vector ) is

    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;

  begin
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant Standard_Complex_Solutions.Solution(dim-integer32(m))
            := Make_Affine(ls.all,m,idz);
      begin
        Standard_Complex_Solutions.Clear(ls);
        ls := new Standard_Complex_Solutions.Solution'(sol);
        Standard_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Make_Affine;

  procedure Make_Affine
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List;
                m : in natural32;
                idz : in Standard_Natural_Vectors.Vector ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant DoblDobl_Complex_Solutions.Solution(dim-integer32(m))
            := Make_Affine(ls.all,m,idz);
      begin
        DoblDobl_Complex_Solutions.Clear(ls);
        ls := new DoblDobl_Complex_Solutions.Solution'(sol);
        DoblDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Make_Affine;
 
  procedure Make_Affine
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List;
                m : in natural32;
                idz : in Standard_Natural_Vectors.Vector ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;

  begin
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      declare
        dim : constant integer32 := ls.n;
        sol : constant QuadDobl_Complex_Solutions.Solution(dim-integer32(m))
            := Make_Affine(ls.all,m,idz);
      begin
        QuadDobl_Complex_Solutions.Clear(ls);
        ls := new QuadDobl_Complex_Solutions.Solution'(sol);
        QuadDobl_Complex_Solutions.Set_Head(tmp,ls);
      end;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Make_Affine;

  function Multi_Projective_Transformation
             ( p : Standard_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    dim : constant integer32 := p'last + integer32(m);
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant Standard_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);

  begin
    res(mhp'range) := mhp;
    if start
     then lhp := Standard_Start_Linear_Polynomials(nbr,m);
     else lhp := Standard_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

  function Multi_Projective_Transformation
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is
 
    dim : constant integer32 := p'last + integer32(m);
    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);

  begin
    res(mhp'range) := mhp;
    if start
     then lhp := DoblDobl_Start_Linear_Polynomials(nbr,m);
     else lhp := DoblDobl_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

  function Multi_Projective_Transformation
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is
 
    dim : constant integer32 := p'last + integer32(m);
    res : TripDobl_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant TripDobl_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : TripDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);

  begin
    res(mhp'range) := mhp;
    if start
     then lhp := TripDobl_Start_Linear_Polynomials(nbr,m);
     else lhp := TripDobl_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

  function Multi_Projective_Transformation
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant integer32 := p'last + integer32(m);
    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);


  begin
    res(mhp'range) := mhp;
    if start
     then lhp := QuadDobl_Start_Linear_Polynomials(nbr,m);
     else lhp := QuadDobl_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

  function Multi_Projective_Transformation
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant integer32 := p'last + integer32(m);
    res : PentDobl_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant PentDobl_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : PentDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);


  begin
    res(mhp'range) := mhp;
    if start
     then lhp := PentDobl_Start_Linear_Polynomials(nbr,m);
     else lhp := PentDobl_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

  function Multi_Projective_Transformation
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant integer32 := p'last + integer32(m);
    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant OctoDobl_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : OctoDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);


  begin
    res(mhp'range) := mhp;
    if start
     then lhp := OctoDobl_Start_Linear_Polynomials(nbr,m);
     else lhp := OctoDobl_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

  function Multi_Projective_Transformation
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant integer32 := p'last + integer32(m);
    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant DecaDobl_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : DecaDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);


  begin
    res(mhp'range) := mhp;
    if start
     then lhp := DecaDobl_Start_Linear_Polynomials(nbr,m);
     else lhp := DecaDobl_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

  function Multi_Projective_Transformation
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys; 
               m : natural32; z : Partition; start : boolean := false )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    dim : constant integer32 := p'last + integer32(m);
    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(1..dim);
    mhp : constant HexaDobl_Complex_Poly_Systems.Poly_Sys
        := Make_Homogeneous(p,m,z);
    lhp : HexaDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    nbr : constant natural32 := natural32(p'last);


  begin
    res(mhp'range) := mhp;
    if start
     then lhp := HexaDobl_Start_Linear_Polynomials(nbr,m);
     else lhp := HexaDobl_Random_Linear_Polynomials(nbr,m,z);
    end if;
    for i in lhp'range loop
      res(p'last+i) := lhp(i);
    end loop;
    return res;
  end Multi_Projective_Transformation;

end Multi_Projective_Transformations;
