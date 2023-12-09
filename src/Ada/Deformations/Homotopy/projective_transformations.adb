with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
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
with Standard_Natural_Vectors;

package body Projective_Transformations is

  function Head_Degree
             ( p : Standard_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use Standard_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Head_Degree
             ( p : DoblDobl_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use DoblDobl_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Head_Degree
             ( p : TripDobl_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use TripDobl_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Head_Degree
             ( p : QuadDobl_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use QuadDobl_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Head_Degree
             ( p : PentDobl_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use PentDobl_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Head_Degree
             ( p : OctoDobl_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use OctoDobl_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Head_Degree
             ( p : DecaDobl_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use DecaDobl_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Head_Degree
             ( p : HexaDobl_Complex_Polynomials.Poly ) return natural32 is

  -- DESCRIPTION :
  --   Returns the degree of the head term.
  --   This is a patch, as the Degree() function causes problems...

    use HexaDobl_Complex_Polynomials;

    res : natural32 := 0;
    trm : constant Term := Head(p);

  begin
    for i in trm.dg'range loop
      res := res + trm.dg(i);
    end loop;
    return res;
  end Head_Degree;

  function Projective_Transformation
             ( p : Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : DoblDobl_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : TripDobl_Complex_Polynomials.Poly )
             return TripDobl_Complex_Polynomials.Poly is

    use TripDobl_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : PentDobl_Complex_Polynomials.Poly )
             return PentDobl_Complex_Polynomials.Poly is

    use PentDobl_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : OctoDobl_Complex_Polynomials.Poly )
             return OctoDobl_Complex_Polynomials.Poly is

    use OctoDobl_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : DecaDobl_Complex_Polynomials.Poly )
             return DecaDobl_Complex_Polynomials.Poly is

    use DecaDobl_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : HexaDobl_Complex_Polynomials.Poly )
             return HexaDobl_Complex_Polynomials.Poly is

    use HexaDobl_Complex_Polynomials;
  
    deg : constant natural32 := Head_Degree(p);
    htd : Degrees
        := new Standard_Natural_Vectors.Vector
                 (1..integer32(Number_of_Unknowns(p))+1);
    res : Poly := Null_Poly;

    procedure Embed_Term ( t : in Term; continue : out boolean ) is

      ht : Term;
      sum : natural32 := 0;

    begin
      ht.cf := t.cf;
      for i in t.dg'range loop
        sum := sum + t.dg(i);
        htd(i) := t.dg(i);
      end loop;
      htd(htd'last) := deg-sum;
      ht.dg := htd;
      Add(res,ht);
      continue := true;
    end Embed_Term;
    procedure Embed_Terms is new Visiting_Iterator(Embed_Term);

  begin
    Embed_Terms(p);
    Clear(htd);
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out TripDobl_Complex_Polynomials.Poly ) is

    use TripDobl_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out PentDobl_Complex_Polynomials.Poly ) is

    use PentDobl_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out OctoDobl_Complex_Polynomials.Poly ) is

    use OctoDobl_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out DecaDobl_Complex_Polynomials.Poly ) is

    use DecaDobl_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out HexaDobl_Complex_Polynomials.Poly ) is

    use HexaDobl_Complex_Polynomials;
  
    res : constant Poly := Projective_Transformation(p);

  begin
    Clear(p); p := res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : TripDobl_Complex_Poly_Systems.Poly_Sys )
             return TripDobl_Complex_Poly_Systems.Poly_Sys is

    res : TripDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : PentDobl_Complex_Poly_Systems.Poly_Sys )
             return PentDobl_Complex_Poly_Systems.Poly_Sys is

    res : PentDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : OctoDobl_Complex_Poly_Systems.Poly_Sys )
             return OctoDobl_Complex_Poly_Systems.Poly_Sys is

    res : OctoDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : DecaDobl_Complex_Poly_Systems.Poly_Sys )
             return DecaDobl_Complex_Poly_Systems.Poly_Sys is

    res : DecaDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( p : HexaDobl_Complex_Poly_Systems.Poly_Sys )
             return HexaDobl_Complex_Poly_Systems.Poly_Sys is

    res : HexaDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := Projective_Transformation(p(k));
    end loop;
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out TripDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out PentDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out OctoDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out DecaDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  procedure Projective_Transformation
              ( p : in out HexaDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for k in p'range loop
      Projective_Transformation(p(k));
    end loop;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : Standard_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(1.0);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : DoblDobl_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);
    one : constant double_double := create(1.0);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(one);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : TripDobl_Complex_Solutions.Solution )
             return TripDobl_Complex_Solutions.Solution is

    use TripDobl_Complex_Numbers;
    use TripDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);
    one : constant triple_double := create(1.0);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(one);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : QuadDobl_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);
    one : constant quad_double := create(1.0);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(one);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : PentDobl_Complex_Solutions.Solution )
             return PentDobl_Complex_Solutions.Solution is

    use PentDobl_Complex_Numbers;
    use PentDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);
    one : constant penta_double := create(1.0);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(one);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : OctoDobl_Complex_Solutions.Solution )
             return OctoDobl_Complex_Solutions.Solution is

    use OctoDobl_Complex_Numbers;
    use OctoDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);
    one : constant octo_double := create(1.0);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(one);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : DecaDobl_Complex_Solutions.Solution )
             return DecaDobl_Complex_Solutions.Solution is

    use DecaDobl_Complex_Numbers;
    use DecaDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);
    one : constant deca_double := create(1.0);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(one);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( s : HexaDobl_Complex_Solutions.Solution )
             return HexaDobl_Complex_Solutions.Solution is

    use HexaDobl_Complex_Numbers;
    use HexaDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n+1);
    one : constant hexa_double := create(1.0);

  begin
    r.v(1..n) := s.v(1..n);
    r.v(n+1) := Create(one);
    r.t := s.t;
    r.m := s.m;
    r.err := s.err;
    r.rco := s.rco;
    r.res := s.res;
    return r;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : Standard_Complex_Solutions.Solution_List )
             return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : DoblDobl_Complex_Solutions.Solution_List )
             return DoblDobl_Complex_Solutions.Solution_List is

    use DoblDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : TripDobl_Complex_Solutions.Solution_List )
             return TripDobl_Complex_Solutions.Solution_List is

    use TripDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : QuadDobl_Complex_Solutions.Solution_List )
             return QuadDobl_Complex_Solutions.Solution_List is

    use QuadDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : PentDobl_Complex_Solutions.Solution_List )
             return PentDobl_Complex_Solutions.Solution_List is

    use PentDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : OctoDobl_Complex_Solutions.Solution_List )
             return OctoDobl_Complex_Solutions.Solution_List is

    use OctoDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : DecaDobl_Complex_Solutions.Solution_List )
             return DecaDobl_Complex_Solutions.Solution_List is

    use DecaDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  function Projective_Transformation
             ( sols : HexaDobl_Complex_Solutions.Solution_List )
             return HexaDobl_Complex_Solutions.Solution_List is

    use HexaDobl_Complex_Solutions;

    res,res_last : Solution_List;
    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Append(res,res_last,Projective_Transformation(ls.all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out TripDobl_Complex_Solutions.Solution_List ) is

    use TripDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out PentDobl_Complex_Solutions.Solution_List ) is

    use PentDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out OctoDobl_Complex_Solutions.Solution_List ) is

    use OctoDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out DecaDobl_Complex_Solutions.Solution_List ) is

    use DecaDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  procedure Projective_Transformation 
              ( sols : in out HexaDobl_Complex_Solutions.Solution_List ) is

    use HexaDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        temp : Solution_List := sols;
        n : constant integer32 := Head_Of(sols).n;
        ls : Link_To_Solution;
        s : Solution(n);
        s2 : Solution(n+1);
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s := ls.all;
          s2 := Projective_Transformation(s);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Projective_Transformation;

  function Affine_Transformation
             ( s : Standard_Complex_Solutions.Solution )
             return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n-1);
    absvn : constant double_float := AbsVal(s.v(n));

  begin
    for i in 1..(n-1) loop
      if absvn + 1.0 = 1.0 
       then r.v(i) := Create(1.0E+10);
       else r.v(i) := s.v(i) / s.v(n);
      end if;
     end loop;
     r.t := s.t;
     r.m := s.m;
     r.err := s.err;
     r.rco := s.rco;
     r.res := s.res;
     return r;
  exception
    when constraint_error =>
       r.v(1..(n-1)) := (1..(n-1) => Create(1.0E+10));
       return r;
  end Affine_Transformation;

  function Affine_Transformation
             ( s : DoblDobl_Complex_Solutions.Solution )
             return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n-1);
    absvn : constant double_double := AbsVal(s.v(n));
    largecst : constant double_double := create(1.0E+20);
    one : constant double_double := create(1.0);

  begin
    for i in 1..(n-1) loop
      if absvn + one = one
       then r.v(i) := Create(largecst);
       else r.v(i) := s.v(i) / s.v(n);
      end if;
     end loop;
     r.t := s.t;
     r.m := s.m;
     r.err := s.err;
     r.rco := s.rco;
     r.res := s.res;
     return r;
  exception
    when constraint_error =>
       r.v(1..(n-1)) := (1..(n-1) => Create(largecst));
       return r;
  end Affine_Transformation;

  function Affine_Transformation
             ( s : QuadDobl_Complex_Solutions.Solution )
             return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    n : constant integer32 := s.n;
    r : Solution(n-1);
    absvn : constant quad_double := AbsVal(s.v(n));
    largecst : constant quad_double := create(1.0E+20);
    one : constant quad_double := create(1.0);

  begin
    for i in 1..(n-1) loop
      if absvn + one = one
       then r.v(i) := Create(largecst);
       else r.v(i) := s.v(i) / s.v(n);
      end if;
     end loop;
     r.t := s.t;
     r.m := s.m;
     r.err := s.err;
     r.rco := s.rco;
     r.res := s.res;
     return r;
  exception
    when constraint_error =>
       r.v(1..(n-1)) := (1..(n-1) => Create(largecst));
       return r;
  end Affine_Transformation;

  procedure Affine_Transformation
              ( sols : in out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        n : constant integer32 := Head_Of(sols).n;
        s1 : Solution(n);
        s2 : Solution(n-1);
        temp : Solution_List := sols;
        ls : Link_To_Solution;
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s1 := ls.all;
          s2 := Affine_Transformation(s1);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Affine_Transformation;

  procedure Affine_Transformation
              ( sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        n : constant integer32 := Head_Of(sols).n;
        s1 : Solution(n);
        s2 : Solution(n-1);
        temp : Solution_List := sols;
        ls : Link_To_Solution;
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s1 := ls.all;
          s2 := Affine_Transformation(s1);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Affine_Transformation;

  procedure Affine_Transformation
              ( sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

  begin
    if Is_Null(sols) then
      null;
    else
      declare
        n : constant integer32 := Head_Of(sols).n;
        s1 : Solution(n);
        s2 : Solution(n-1);
        temp : Solution_List := sols;
        ls : Link_To_Solution;
      begin
        while not Is_Null(temp) loop
          ls := Head_Of(temp);
          s1 := ls.all;
          s2 := Affine_Transformation(s1);
          Clear(ls);
          ls := new Solution'(s2);
          Set_Head(temp,ls);
          temp := Tail_Of(temp);
        end loop;
      end;
    end if;
  end Affine_Transformation;

end Projective_Transformations;
