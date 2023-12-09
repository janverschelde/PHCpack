with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Penta_Double_Numbers;               use Penta_Double_Numbers;
with Octo_Double_Numbers;                use Octo_Double_Numbers;
with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with TripDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with PentDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with TripDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with PentDobl_Random_Numbers;
with OctoDobl_Random_Numbers;
with DecaDobl_Random_Numbers;
with HexaDobl_Random_Numbers;
with Standard_Natural_Vectors;

package body Homogenization is

  function Homogeneous_Part
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;
    d : constant integer32 := Degree(p);

    procedure Homogeneous_Term ( t : in Term; continue : out boolean ) is
    begin
      if integer32(Sum(t.dg)) = d 
       then continue := true; Add(res,t);
       else continue := false;
      end if;
    end Homogeneous_Term;
    procedure Homogeneous_Terms is new Visiting_Iterator(Homogeneous_Term);

  begin
    Homogeneous_Terms(p);
    return res;
  end Homogeneous_Part;

  function Homogeneous_Part
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Poly_Systems;

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Homogeneous_Part(p(i));
    end loop;
    return res;
  end Homogeneous_Part;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : double_float;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := Standard_Random_Numbers.Random;
    t.cf := Standard_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := Standard_Random_Numbers.Random;
      t.cf := Standard_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : double_double;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := DoblDobl_Random_Numbers.Random;
    t.cf := DoblDobl_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := DoblDobl_Random_Numbers.Random;
      t.cf := DoblDobl_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return TripDobl_Complex_Polynomials.Poly is

    use TripDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : triple_double;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := TripDobl_Random_Numbers.Random;
    t.cf := TripDobl_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := TripDobl_Random_Numbers.Random;
      t.cf := TripDobl_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : quad_double;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := QuadDobl_Random_Numbers.Random;
    t.cf := QuadDobl_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := QuadDobl_Random_Numbers.Random;
      t.cf := QuadDobl_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return PentDobl_Complex_Polynomials.Poly is

    use PentDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : penta_double;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := PentDobl_Random_Numbers.Random;
    t.cf := PentDobl_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := PentDobl_Random_Numbers.Random;
      t.cf := PentDobl_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return OctoDobl_Complex_Polynomials.Poly is

    use OctoDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : octo_double;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := OctoDobl_Random_Numbers.Random;
    t.cf := OctoDobl_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := OctoDobl_Random_Numbers.Random;
      t.cf := OctoDobl_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return DecaDobl_Complex_Polynomials.Poly is

    use DecaDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : deca_double;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := DecaDobl_Random_Numbers.Random;
    t.cf := DecaDobl_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := DecaDobl_Random_Numbers.Random;
      t.cf := DecaDobl_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Real_Random_Hyperplane
             ( n : natural32 )
             return HexaDobl_Complex_Polynomials.Poly is

    use HexaDobl_Complex_Polynomials;

    res : Poly;
    t : Term;
    ranflt : hexa_double;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    ranflt := HexaDobl_Random_Numbers.Random;
    t.cf := HexaDobl_Complex_Numbers.Create(ranflt);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      ranflt := HexaDobl_Random_Numbers.Random;
      t.cf := HexaDobl_Complex_Numbers.Create(ranflt);
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Real_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := Standard_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := Standard_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := DoblDobl_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := DoblDobl_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return TripDobl_Complex_Polynomials.Poly is

    use TripDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := TripDobl_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := TripDobl_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := QuadDobl_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := QuadDobl_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return PentDobl_Complex_Polynomials.Poly is

    use PentDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := PentDobl_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := PentDobl_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return OctoDobl_Complex_Polynomials.Poly is

    use OctoDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := OctoDobl_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := OctoDobl_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return DecaDobl_Complex_Polynomials.Poly is

    use DecaDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := DecaDobl_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := DecaDobl_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Complex_Random_Hyperplane
             ( n : natural32 )
             return HexaDobl_Complex_Polynomials.Poly is

    use HexaDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := HexaDobl_Random_Numbers.Random1;
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    for i in 1..integer32(n) loop
      t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
      t.dg(i) := 1;
      t.cf := HexaDobl_Random_Numbers.Random1;
      Add(res,t);
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    end loop;
    return res;
  end Complex_Random_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use Standard_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := Standard_Complex_Numbers.Create(-1.0);
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := Standard_Complex_Numbers.Create(1.0);
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return DoblDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use DoblDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := DoblDobl_Complex_Numbers.Create(integer32(-1));
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := DoblDobl_Complex_Numbers.Create(integer32(1));
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return TripDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use TripDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := TripDobl_Complex_Numbers.Create(integer32(-1));
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := TripDobl_Complex_Numbers.Create(integer32(1));
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return QuadDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use QuadDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := QuadDobl_Complex_Numbers.Create(integer32(-1));
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := QuadDobl_Complex_Numbers.Create(integer32(1));
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return PentDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use PentDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := PentDobl_Complex_Numbers.Create(integer32(-1));
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := PentDobl_Complex_Numbers.Create(integer32(1));
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return OctoDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use OctoDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := OctoDobl_Complex_Numbers.Create(integer32(-1));
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := OctoDobl_Complex_Numbers.Create(integer32(1));
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return DecaDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use DecaDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := DecaDobl_Complex_Numbers.Create(integer32(-1));
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := DecaDobl_Complex_Numbers.Create(integer32(1));
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  function Standard_Hyperplane
             ( n,i : natural32 )
             return HexaDobl_Complex_Polynomials.Poly is

  -- DESCRIPTION : Returns x_i - 1.

    use HexaDobl_Complex_Polynomials;

    res : Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.cf := HexaDobl_Complex_Numbers.Create(integer32(-1));
    res := Create(t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    t.dg(integer32(i)) := 1;
    t.cf := HexaDobl_Complex_Numbers.Create(integer32(1));
    Add(res,t);
    Standard_Natural_Vectors.Clear
      (Standard_Natural_Vectors.Link_to_Vector(t.dg));
    return res;
  end Standard_Hyperplane;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out Standard_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      Standard_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      DoblDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out TripDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      TripDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      QuadDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out PentDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      PentDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out OctoDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      OctoDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out DecaDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      DecaDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Real_Random_Hyperplanes
              ( s : in out HexaDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   the polynomial system s will be filled with polynomials in m unknowns,
  --   with real coefficients.

  begin
    for i in s'range loop
      HexaDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Real_Random_Hyperplane(m);
    end loop;
  end Construct_Real_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out Standard_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      Standard_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      DoblDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out TripDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      TripDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      QuadDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out PentDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      PentDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out OctoDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      OctoDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out DecaDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      DecaDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Complex_Random_Hyperplanes
              ( s : in out HexaDobl_Complex_Poly_Systems.Poly_Sys;
                m : natural32 ) is

  -- DESCRIPTION :
  --   The polynomial system s will be filled with polynomials in m unknowns,
  --   with complex coefficients.

  begin
    for i in s'range loop
      HexaDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Complex_Random_Hyperplane(m);
    end loop;
  end Construct_Complex_Random_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      Standard_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      DoblDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out TripDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      TripDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      QuadDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out PentDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      PentDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out OctoDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      OctoDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out DecaDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      DecaDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Construct_Standard_Hyperplanes
              ( s : in out HexaDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   for i in s'range : s(i) := x_i - 1.

    n : constant natural32 := natural32(s'length);

  begin
    for i in s'range loop
      HexaDobl_Complex_Polynomials.Clear(s(i));
      s(i) := Standard_Hyperplane(n,natural32(i));
    end loop;
  end Construct_Standard_Hyperplanes;

  procedure Enlarge_Before
              ( p : in out Standard_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use Standard_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_Before
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use DoblDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_Before
              ( p : in out TripDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use TripDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_Before
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use QuadDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_Before
              ( p : in out PentDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use PentDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_Before
              ( p : in out OctoDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use OctoDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_Before
              ( p : in out DecaDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use DecaDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_Before
              ( p : in out HexaDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be inserted to t.dg

    use HexaDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in 1..integer32(m) loop
        d(i) := 0;
      end loop;
      for i in t.dg'range loop
        d(i+integer32(m)) := t.dg(i);
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_Before;

  procedure Enlarge_After
              ( p : in out Standard_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use Standard_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  procedure Enlarge_After
              ( p : in out DoblDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use DoblDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  procedure Enlarge_After
              ( p : in out TripDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use TripDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  procedure Enlarge_After
              ( p : in out QuadDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use QuadDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  procedure Enlarge_After
              ( p : in out PentDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use PentDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  procedure Enlarge_After
              ( p : in out OctoDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use OctoDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  procedure Enlarge_After
              ( p : in out DecaDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use DecaDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  procedure Enlarge_After
              ( p : in out HexaDobl_Complex_Polynomials.Poly;
                m : in natural32 ) is

  -- DESCRIPTION :
  --   To each term t of p, m additional zero entries will be added to t.dg

    use HexaDobl_Complex_Polynomials;

    procedure Enlarge_Term ( t : in out Term; continue : out boolean ) is

      d : constant Degrees
        := new Standard_Natural_Vectors.Vector(1..(t.dg'last+integer32(m)));

    begin
      for i in t.dg'range loop
        d(i) := t.dg(i);
      end loop;
      for i in (t.dg'last+1)..integer32(m) loop
        d(i) := 0;
      end loop;
      Standard_Natural_Vectors.Clear
        (Standard_Natural_Vectors.Link_to_Vector(t.dg));
      t.dg := d;
      continue := true;
    end Enlarge_Term;
    procedure Enlarge_Terms is new Changing_Iterator(Enlarge_Term);

  begin
    Enlarge_Terms(p);
  end Enlarge_After;

  function Add_Equations
             ( s1 : Standard_Complex_Poly_Systems.Poly_Sys;
               s2 : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equations
             ( s1 : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s2 : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equations
             ( s1 : TripDobl_Complex_Poly_Systems.Poly_Sys;
               s2 : TripDobl_Complex_Poly_Systems.Poly_Sys )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    use TripDobl_Complex_Polynomials;
    use TripDobl_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equations
             ( s1 : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s2 : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equations
             ( s1 : PentDobl_Complex_Poly_Systems.Poly_Sys;
               s2 : PentDobl_Complex_Poly_Systems.Poly_Sys )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    use PentDobl_Complex_Polynomials;
    use PentDobl_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equations
             ( s1 : OctoDobl_Complex_Poly_Systems.Poly_Sys;
               s2 : OctoDobl_Complex_Poly_Systems.Poly_Sys )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    use OctoDobl_Complex_Polynomials;
    use OctoDobl_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equations
             ( s1 : DecaDobl_Complex_Poly_Systems.Poly_Sys;
               s2 : DecaDobl_Complex_Poly_Systems.Poly_Sys )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    use DecaDobl_Complex_Polynomials;
    use DecaDobl_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equations
             ( s1 : HexaDobl_Complex_Poly_Systems.Poly_Sys;
               s2 : HexaDobl_Complex_Poly_Systems.Poly_Sys )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    use HexaDobl_Complex_Polynomials;
    use HexaDobl_Complex_Poly_Systems;

    n1 : constant natural32 := natural32(s1'last - s1'first + 1);
    n2 : constant natural32 := natural32(s2'last - s2'first + 1);
    res : Poly_Sys(1..integer32(n1+n2));
    m : natural32;

  begin
    for i in 1..integer32(n1) loop
      Copy(s1(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    for i in 1..integer32(n2) loop
      Copy(s2(i),res(integer32(n1)+i));
      m := Number_Of_Unknowns(res(integer32(n1)+i));
      if m < (n1+n2)
       then m := n1 + n2 - m; Enlarge_Before(res(integer32(n1)+i),m);
      end if;
    end loop;
    return res;
  end Add_Equations;

  function Add_Equation
             ( s : Standard_Complex_Poly_Systems.Poly_Sys;
               p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'last - s'first + 1);
    m : natural32;
    res : Poly_Sys(1..integer32(n+1));

  begin
    for i in 1..integer32(n) loop
      Copy(s(i),res(i));
      m := Number_Of_Unknowns(res(i));
      if m < n+1
       then m := n + 1 - m; Enlarge_After(res(i),m);
      end if;
    end loop;
    m := Number_Of_Unknowns(p);
    if m < (n+1)
     then m := n + 1 - m; Enlarge_Before(res(integer32(n+1)),m);
    end if;
    return res;
  end Add_Equation;

  function Add_Random_Hyperplanes
             ( s : Standard_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Random_Hyperplanes
             ( s : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Random_Hyperplanes
             ( s : TripDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    use TripDobl_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Random_Hyperplanes
             ( s : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Random_Hyperplanes
             ( s : PentDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    use PentDobl_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Random_Hyperplanes
             ( s : OctoDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    use OctoDobl_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Random_Hyperplanes
             ( s : DecaDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    use DecaDobl_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Random_Hyperplanes
             ( s : HexaDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32; re : boolean )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    use HexaDobl_Complex_Poly_Systems;

    s2 : Poly_Sys(1..integer32(m));
    n : constant natural32 := natural32(s'last - s'first + 1);
    res : Poly_Sys(1..integer32(m+n));

  begin
    if re
     then Construct_Real_Random_Hyperplanes(s2,m+n);
     else Construct_Complex_Random_Hyperplanes(s2,m+n);
    end if;
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Random_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : Standard_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : TripDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    use TripDobl_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : PentDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    use PentDobl_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : OctoDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    use OctoDobl_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : DecaDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    use DecaDobl_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

  function Add_Standard_Hyperplanes
             ( s : HexaDobl_Complex_Poly_Systems.Poly_Sys;
               m : natural32 )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    use HexaDobl_Complex_Poly_Systems;

    n : constant natural32 := natural32(s'length);
    res : Poly_Sys(1..integer32(n+m));
    s2 : Poly_Sys(1..integer32(m));

  begin
    Construct_Standard_Hyperplanes(s2);
    res := Add_Equations(s,s2);
    Clear(s2);
    return res;
  end Add_Standard_Hyperplanes;

end Homogenization;
