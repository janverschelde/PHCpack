with Communications_with_User;           use Communications_with_User;
with Numbers_io;                         use Numbers_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Symbol_Table;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;   use QuadDobl_Complex_Laur_Systems_io;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets;

package body Square_and_Embed_Systems is

  function Restrict ( t : Standard_Complex_Polynomials.Term;
                      m,k : integer32 )
                    return Standard_Complex_Polynomials.Term is

    use Standard_Complex_Polynomials;

    rt : Term;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
    for i in m+1..rt.dg'last-k loop
      rt.dg(i) := 0;
    end loop;
    return rt;
  end Restrict;

  function Restrict ( t : Standard_Complex_Laurentials.Term;
                      m,k : integer32 )
                    return Standard_Complex_Laurentials.Term is

    use Standard_Complex_Laurentials;

    rt : Term;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
    for i in m+1..rt.dg'last-k loop
      rt.dg(i) := 0;
    end loop;
    return rt;
  end Restrict;

  function Restrict ( t : DoblDobl_Complex_Polynomials.Term;
                      m,k : integer32 )
                    return DoblDobl_Complex_Polynomials.Term is

    use DoblDobl_Complex_Polynomials;

    rt : Term;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
    for i in m+1..rt.dg'last-k loop
      rt.dg(i) := 0;
    end loop;
    return rt;
  end Restrict;

  function Restrict ( t : DoblDobl_Complex_Laurentials.Term;
                      m,k : integer32 )
                    return DoblDobl_Complex_Laurentials.Term is

    use DoblDobl_Complex_Laurentials;

    rt : Term;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
    for i in m+1..rt.dg'last-k loop
      rt.dg(i) := 0;
    end loop;
    return rt;
  end Restrict;

  function Restrict ( t : QuadDobl_Complex_Polynomials.Term;
                      m,k : integer32 )
                    return QuadDobl_Complex_Polynomials.Term is

    use QuadDobl_Complex_Polynomials;

    rt : Term;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Natural_Vectors.Vector'(t.dg.all);
    for i in m+1..rt.dg'last-k loop
      rt.dg(i) := 0;
    end loop;
    return rt;
  end Restrict;

  function Restrict ( t : QuadDobl_Complex_Laurentials.Term;
                      m,k : integer32 )
                    return QuadDobl_Complex_Laurentials.Term is

    use QuadDobl_Complex_Laurentials;

    rt : Term;

  begin
    rt.cf := t.cf;
    rt.dg := new Standard_Integer_Vectors.Vector'(t.dg.all);
    for i in m+1..rt.dg'last-k loop
      rt.dg(i) := 0;
    end loop;
    return rt;
  end Restrict;

  function Restrict ( p : Standard_Complex_Polynomials.Poly;
                      m,k : integer32 )
                    return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Restrict_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Restrict(t,m,k);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Restrict_Term;
    procedure Restrict_Terms is new Visiting_Iterator(Restrict_Term);      

  begin
    Restrict_Terms(p);
    return res;
  end Restrict;

  function Restrict ( p : Standard_Complex_Laurentials.Poly;
                      m,k : integer32 )
                    return Standard_Complex_Laurentials.Poly is

    use Standard_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Restrict_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Restrict(t,m,k);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Restrict_Term;
    procedure Restrict_Terms is new Visiting_Iterator(Restrict_Term);      

  begin
    Restrict_Terms(p);
    return res;
  end Restrict;

  function Restrict ( p : DoblDobl_Complex_Polynomials.Poly;
                      m,k : integer32 )
                    return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Restrict_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Restrict(t,m,k);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Restrict_Term;
    procedure Restrict_Terms is new Visiting_Iterator(Restrict_Term);      

  begin
    Restrict_Terms(p);
    return res;
  end Restrict;

  function Restrict ( p : DoblDobl_Complex_Laurentials.Poly;
                      m,k : integer32 )
                    return DoblDobl_Complex_Laurentials.Poly is

    use DoblDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Restrict_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Restrict(t,m,k);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Restrict_Term;
    procedure Restrict_Terms is new Visiting_Iterator(Restrict_Term);      

  begin
    Restrict_Terms(p);
    return res;
  end Restrict;

  function Restrict ( p : QuadDobl_Complex_Polynomials.Poly;
                      m,k : integer32 )
                    return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

    res : Poly := Null_Poly;

    procedure Restrict_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Restrict(t,m,k);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Restrict_Term;
    procedure Restrict_Terms is new Visiting_Iterator(Restrict_Term);      

  begin
    Restrict_Terms(p);
    return res;
  end Restrict;

  function Restrict ( p : QuadDobl_Complex_Laurentials.Poly;
                      m,k : integer32 )
                    return QuadDobl_Complex_Laurentials.Poly is

    use QuadDobl_Complex_Laurentials;

    res : Poly := Null_Poly;

    procedure Restrict_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Restrict(t,m,k);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Restrict_Term;
    procedure Restrict_Terms is new Visiting_Iterator(Restrict_Term);      

  begin
    Restrict_Terms(p);
    return res;
  end Restrict;

  function Maximum ( a,b : natural32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the maximum of a and b.

  begin
    if a > b
     then return a;
     else return b;
    end if;
  end Maximum;

  procedure Interactive_Embed_Square_System 
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    k,m : natural32 := 0;
    ans : character;

  begin
    put("Give the expected top dimension : "); Read_Natural(k);
    Witness_Sets_io.Add_Embed_Symbols(k); topdim := k;
    declare
      ep : Poly_Sys(p'first..p'last+integer32(k));
    begin
      ep := Slice_and_Embed(p,k);
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-integer32(k)+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),integer32(k));
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
      end if;
      put_line(file,ep);
      embsys := new Poly_Sys'(ep);
    end;
  end Interactive_Embed_Square_System;

  procedure Interactive_Embed_Square_System 
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                embsys : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                topdim : out natural32 ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;

    k,m : natural32 := 0;
    ans : character;

  begin
    put("Give the expected top dimension : "); Read_Natural(k);
    Witness_Sets_io.Add_Embed_Symbols(k); topdim := k;
    declare
      ep : Laur_Sys(p'first..p'last+integer32(k));
    begin
      ep := Slice_and_Embed(p,k);
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-integer32(k)+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),integer32(k));
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
      end if;
      put_line(file,ep);
      embsys := new Laur_Sys'(ep);
    end;
  end Interactive_Embed_Square_System;

  procedure Interactive_Embed_Square_System 
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                embsys : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    k,m : natural32 := 0;
    ans : character;

  begin
    put("Give the expected top dimension : "); Read_Natural(k);
    Witness_Sets_io.Add_Embed_Symbols(k); topdim := k;
    declare
      ep : Poly_Sys(p'first..p'last+integer32(k));
    begin
      ep := Slice_and_Embed(p,k);
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-integer32(k)+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),integer32(k));
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
      end if;
      put_line(file,ep);
      embsys := new Poly_Sys'(ep);
    end;
  end Interactive_Embed_Square_System;

  procedure Interactive_Embed_Square_System 
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                embsys : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                topdim : out natural32 ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;

    k,m : natural32 := 0;
    ans : character;

  begin
    put("Give the expected top dimension : "); Read_Natural(k);
    Witness_Sets_io.Add_Embed_Symbols(k); topdim := k;
    declare
      ep : Laur_Sys(p'first..p'last+integer32(k));
    begin
      ep := Slice_and_Embed(p,k);
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-integer32(k)+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),integer32(k));
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
      end if;
      put_line(file,ep);
      embsys := new Laur_Sys'(ep);
    end;
  end Interactive_Embed_Square_System;

  procedure Interactive_Embed_Square_System 
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                embsys : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    k,m : natural32 := 0;
    ans : character;

  begin
    put("Give the expected top dimension : "); Read_Natural(k);
    Witness_Sets_io.Add_Embed_Symbols(k); topdim := k;
    declare
      ep : Poly_Sys(p'first..p'last+integer32(k));
    begin
      ep := Slice_and_Embed(p,k);
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-integer32(k)+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),integer32(k));
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
      end if;
      put_line(file,ep);
      embsys := new Poly_Sys'(ep);
    end;
  end Interactive_Embed_Square_System;

  procedure Interactive_Embed_Square_System 
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                embsys : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                topdim : out natural32 ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;

    k,m : natural32 := 0;
    ans : character;

  begin
    put("Give the expected top dimension : "); Read_Natural(k);
    Witness_Sets_io.Add_Embed_Symbols(k); topdim := k;
    declare
      ep : Laur_Sys(p'first..p'last+integer32(k));
    begin
      ep := Slice_and_Embed(p,k);
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-integer32(k)+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),integer32(k));
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
      end if;
      put_line(file,ep);
      embsys := new Laur_Sys'(ep);
    end;
  end Interactive_Embed_Square_System;

  procedure Embed_Square_System 
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

  begin
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Poly_Sys(p'first..p'last+integer32(topdim));
    begin
      ep := Slice_and_Embed(p,topdim);
      embsys := new Poly_Sys'(ep);
    end;
  end Embed_Square_System;

  procedure Embed_Square_System 
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                topdim : in natural32;
                embsys : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use Standard_Complex_Laur_Systems;

  begin
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Laur_Sys(p'first..p'last+integer32(topdim));
    begin
      ep := Slice_and_Embed(p,topdim);
      embsys := new Laur_Sys'(ep);
    end;
  end Embed_Square_System;

  procedure Embed_Square_System 
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                embsys : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

  begin
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Poly_Sys(p'first..p'last+integer32(topdim));
    begin
      ep := Slice_and_Embed(p,topdim);
      embsys := new Poly_Sys'(ep);
    end;
  end Embed_Square_System;

  procedure Embed_Square_System 
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                topdim : in natural32;
                embsys : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use DoblDobl_Complex_Laur_Systems;

  begin
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Laur_Sys(p'first..p'last+integer32(topdim));
    begin
      ep := Slice_and_Embed(p,topdim);
      embsys := new Laur_Sys'(ep);
    end;
  end Embed_Square_System;

  procedure Embed_Square_System 
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                embsys : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

  begin
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Poly_Sys(p'first..p'last+integer32(topdim));
    begin
      ep := Slice_and_Embed(p,topdim);
      embsys := new Poly_Sys'(ep);
    end;
  end Embed_Square_System;

  procedure Embed_Square_System 
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                topdim : in natural32;
                embsys : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use QuadDobl_Complex_Laur_Systems;

  begin
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Laur_Sys(p'first..p'last+integer32(topdim));
    begin
      ep := Slice_and_Embed(p,topdim);
      embsys := new Laur_Sys'(ep);
    end;
  end Embed_Square_System;

  function Full_Embed_Nonsquare_System
              ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                nq,nv,k : natural32 )
              return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Poly_Systems;
 
    embedded : Poly_Sys(p'first..p'last+integer32(k));
    d : constant integer32 := integer32(nv - nq);

  begin
    if nv < nq then
      declare
        sp : constant Poly_Sys := Square(p);
      begin
        embedded := Slice_and_Embed(sp,k);
      end;
    else
      if integer32(k) <= d then
        embedded := Embed_with_Dummies(p,k);
      else
        declare
          aux : Poly_Sys(p'first..p'last+d);
        begin
          aux := Embed_with_Dummies(p,natural32(d));
          embedded := Slice_and_Embed(aux,k-natural32(d));
        end;
      end if;
    end if;
    return embedded;
  end Full_Embed_Nonsquare_System;

  function Full_Embed_Nonsquare_System
              ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                nq,nv,k : natural32 )
              return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laur_Systems;
 
    embedded : Laur_Sys(p'first..p'last+integer32(k));
    d : constant integer32 := integer32(nv - nq);

  begin
    if nv < nq then
      declare
        sp : constant Laur_Sys := Square(p);
      begin
        embedded := Slice_and_Embed(sp,k);
      end;
    else
      if integer32(k) <= d then
        embedded := Embed_with_Dummies(p,k);
      else
        declare
          aux : Laur_Sys(p'first..p'last+d);
        begin
          aux := Embed_with_Dummies(p,natural32(d));
          embedded := Slice_and_Embed(aux,k-natural32(d));
        end;
      end if;
    end if;
    return embedded;
  end Full_Embed_Nonsquare_System;

  function Full_Embed_Nonsquare_System
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                nq,nv,k : natural32 )
              return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;
 
    embedded : Poly_Sys(p'first..p'last+integer32(k));
    d : constant integer32 := integer32(nv - nq);

  begin
    if nv < nq then
      declare
        sp : constant Poly_Sys := Square(p);
      begin
        embedded := Slice_and_Embed(sp,k);
      end;
    else
      if integer32(k) <= d then
        embedded := Embed_with_Dummies(p,k);
      else
        declare
          aux : Poly_Sys(p'first..p'last+d);
        begin
          aux := Embed_with_Dummies(p,natural32(d));
          embedded := Slice_and_Embed(aux,k-natural32(d));
        end;
      end if;
    end if;
    return embedded;
  end Full_Embed_Nonsquare_System;

  function Full_Embed_Nonsquare_System
              ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv,k : natural32 )
              return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laur_Systems;
 
    embedded : Laur_Sys(p'first..p'last+integer32(k));
    d : constant integer32 := integer32(nv - nq);

  begin
    if nv < nq then
      declare
        sp : constant Laur_Sys := Square(p);
      begin
        embedded := Slice_and_Embed(sp,k);
      end;
    else
      if integer32(k) <= d then
        embedded := Embed_with_Dummies(p,k);
      else
        declare
          aux : Laur_Sys(p'first..p'last+d);
        begin
          aux := Embed_with_Dummies(p,natural32(d));
          embedded := Slice_and_Embed(aux,k-natural32(d));
        end;
      end if;
    end if;
    return embedded;
  end Full_Embed_Nonsquare_System;

  function Full_Embed_Nonsquare_System
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                nq,nv,k : natural32 )
              return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;
 
    embedded : Poly_Sys(p'first..p'last+integer32(k));
    d : constant integer32 := integer32(nv - nq);

  begin
    if nv < nq then
      declare
        sp : constant Poly_Sys := Square(p);
      begin
        embedded := Slice_and_Embed(sp,k);
      end;
    else
      if integer32(k) <= d then
        embedded := Embed_with_Dummies(p,k);
      else
        declare
          aux : Poly_Sys(p'first..p'last+d);
        begin
          aux := Embed_with_Dummies(p,natural32(d));
          embedded := Slice_and_Embed(aux,k-natural32(d));
        end;
      end if;
    end if;
    return embedded;
  end Full_Embed_Nonsquare_System;

  function Full_Embed_Nonsquare_System
              ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nq,nv,k : natural32 )
              return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laur_Systems;
 
    embedded : Laur_Sys(p'first..p'last+integer32(k));
    d : constant integer32 := integer32(nv - nq);

  begin
    if nv < nq then
      declare
        sp : constant Laur_Sys := Square(p);
      begin
        embedded := Slice_and_Embed(sp,k);
      end;
    else
      if integer32(k) <= d then
        embedded := Embed_with_Dummies(p,k);
      else
        declare
          aux : Laur_Sys(p'first..p'last+d);
        begin
          aux := Embed_with_Dummies(p,natural32(d));
          embedded := Slice_and_Embed(aux,k-natural32(d));
        end;
      end if;
    end if;
    return embedded;
  end Full_Embed_Nonsquare_System;

  procedure Interactive_Embed_Nonsquare_System
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk : in natural32;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 ) is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Poly_Sys(1..integer32(max)) := Square(p);
    m,ns : natural32;
    k,a : integer32 := 0;
    ans : character;

  begin
    if nbequ > nbunk
     then Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
    end if;
    put("Give the expected top dimension : "); Read_Integer(k);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(natural32(k));
    topdim := natural32(k);
    declare
      ep : Poly_Sys(sp'first..sp'last+k)
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,natural32(k));
    begin
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-k+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),k);
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
        a := integer32(nbunk - nbequ) - k;
        if a > 0 then
          for i in ep'last-2*k-a+1..ep'last-2*k loop
            declare
              rp : constant Poly := Restrict(ep(i),integer32(m),k);
            begin
              Clear(ep(i));
              ep(i) := rp;
            end;
          end loop;
        end if; 
      end if;
      put_line(file,ep);
      embsys := new Poly_Sys'(ep);
    end;
  end Interactive_Embed_Nonsquare_System;

  procedure Interactive_Embed_Nonsquare_System
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                nbequ,nbunk : in natural32;
                embsys : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                topdim : out natural32 ) is

    use Standard_Complex_Laurentials;
    use Standard_Complex_Laur_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Laur_Sys(1..integer32(max)) := Square(p);
    m,ns : natural32;
    k,a : integer32 := 0;
    ans : character;

  begin
    if nbequ > nbunk
     then Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
    end if;
    put("Give the expected top dimension : "); Read_Integer(k);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(natural32(k));
    topdim := natural32(k);
    declare
      ep : Laur_Sys(sp'first..sp'last+k)
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,natural32(k));
    begin
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-k+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),k);
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
        a := integer32(nbunk - nbequ) - k;
        if a > 0 then
          for i in ep'last-2*k-a+1..ep'last-2*k loop
            declare
              rp : constant Poly := Restrict(ep(i),integer32(m),k);
            begin
              Clear(ep(i));
              ep(i) := rp;
            end;
          end loop;
        end if; 
      end if;
      put_line(file,ep);
      embsys := new Laur_Sys'(ep);
    end;
  end Interactive_Embed_Nonsquare_System;

  procedure Interactive_Embed_Nonsquare_System
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk : in natural32;
                embsys : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 ) is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Poly_Sys(1..integer32(max)) := Square(p);
    m,ns : natural32;
    k,a : integer32 := 0;
    ans : character;

  begin
    if nbequ > nbunk
     then Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
    end if;
    put("Give the expected top dimension : "); Read_Integer(k);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(natural32(k));
    topdim := natural32(k);
    declare
      ep : Poly_Sys(sp'first..sp'last+k)
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,natural32(k));
    begin
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-k+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),k);
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
        a := integer32(nbunk - nbequ) - k;
        if a > 0 then
          for i in ep'last-2*k-a+1..ep'last-2*k loop
            declare
              rp : constant Poly := Restrict(ep(i),integer32(m),k);
            begin
              Clear(ep(i));
              ep(i) := rp;
            end;
          end loop;
        end if; 
      end if;
      put_line(file,ep);
      embsys := new Poly_Sys'(ep);
    end;
  end Interactive_Embed_Nonsquare_System;

  procedure Interactive_Embed_Nonsquare_System
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nbequ,nbunk : in natural32;
                embsys : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                topdim : out natural32 ) is

    use DoblDobl_Complex_Laurentials;
    use DoblDobl_Complex_Laur_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Laur_Sys(1..integer32(max)) := Square(p);
    m,ns : natural32;
    k,a : integer32 := 0;
    ans : character;

  begin
    if nbequ > nbunk
     then Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
    end if;
    put("Give the expected top dimension : "); Read_Integer(k);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(natural32(k));
    topdim := natural32(k);
    declare
      ep : Laur_Sys(sp'first..sp'last+k)
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,natural32(k));
    begin
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-k+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),k);
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
        a := integer32(nbunk - nbequ) - k;
        if a > 0 then
          for i in ep'last-2*k-a+1..ep'last-2*k loop
            declare
              rp : constant Poly := Restrict(ep(i),integer32(m),k);
            begin
              Clear(ep(i));
              ep(i) := rp;
            end;
          end loop;
        end if; 
      end if;
      put_line(file,ep);
      embsys := new Laur_Sys'(ep);
    end;
  end Interactive_Embed_Nonsquare_System;

  procedure Interactive_Embed_Nonsquare_System
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk : in natural32;
                embsys : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                topdim : out natural32 ) is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Poly_Sys(1..integer32(max)) := Square(p);
    m,ns : natural32;
    k,a : integer32 := 0;
    ans : character;

  begin
    if nbequ > nbunk
     then Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
    end if;
    put("Give the expected top dimension : "); Read_Integer(k);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(natural32(k));
    topdim := natural32(k);
    declare
      ep : Poly_Sys(sp'first..sp'last+k)
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,natural32(k));
    begin
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-k+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),k);
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
        a := integer32(nbunk - nbequ) - k;
        if a > 0 then
          for i in ep'last-2*k-a+1..ep'last-2*k loop
            declare
              rp : constant Poly := Restrict(ep(i),integer32(m),k);
            begin
              Clear(ep(i));
              ep(i) := rp;
            end;
          end loop;
        end if; 
      end if;
      put_line(file,ep);
      embsys := new Poly_Sys'(ep);
    end;
  end Interactive_Embed_Nonsquare_System;

  procedure Interactive_Embed_Nonsquare_System
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nbequ,nbunk : in natural32;
                embsys : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                topdim : out natural32 ) is

    use QuadDobl_Complex_Laurentials;
    use QuadDobl_Complex_Laur_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Laur_Sys(1..integer32(max)) := Square(p);
    m,ns : natural32;
    k,a : integer32 := 0;
    ans : character;

  begin
    if nbequ > nbunk
     then Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
    end if;
    put("Give the expected top dimension : "); Read_Integer(k);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(natural32(k));
    topdim := natural32(k);
    declare
      ep : Laur_Sys(sp'first..sp'last+k)
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,natural32(k));
    begin
      put("Should the slices be restricted to a subspace ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put("Give the dimension of the subspace : "); Read_Natural(m);
        put("The first "); put(m,1);
        put_line(" variables span the subspace...");
        Witness_Sets_io.Determine_Order(ep);
        for i in ep'last-k+1..ep'last loop
          declare
            rp : constant Poly := Restrict(ep(i),integer32(m),k);
          begin
            Clear(ep(i));
            ep(i) := rp;
          end;
        end loop;
        a := integer32(nbunk - nbequ) - k;
        if a > 0 then
          for i in ep'last-2*k-a+1..ep'last-2*k loop
            declare
              rp : constant Poly := Restrict(ep(i),integer32(m),k);
            begin
              Clear(ep(i));
              ep(i) := rp;
            end;
          end loop;
        end if; 
      end if;
      put_line(file,ep);
      embsys := new Laur_Sys'(ep);
    end;
  end Interactive_Embed_Nonsquare_System;

  procedure Embed_Nonsquare_System
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk,topdim : in natural32;
                embsys : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Poly_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Poly_Sys(1..integer32(max)) := Square(p);
    ns : natural32;

  begin
    if nbequ > nbunk then
      Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
     -- put("added "); put(nbequ - nbunk,1); put_line(" slack variables");
    end if;
   -- put_line("The squared polynomial system :"); put_line(sp);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Poly_Sys(sp'first..sp'last+integer32(topdim))
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,topdim);
    begin
      embsys := new Poly_Sys'(ep);
    end;
  end Embed_Nonsquare_System;

  procedure Embed_Nonsquare_System
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                nbequ,nbunk,topdim : in natural32;
                embsys : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use Standard_Complex_Laur_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Laur_Sys(1..integer32(max)) := Square(p);
    ns : natural32;

  begin
    if nbequ > nbunk then
      Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
     -- put("added "); put(nbequ - nbunk,1); put_line(" slack variables");
    end if;
   -- put_line("The squared polynomial system :"); put_line(sp);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Laur_Sys(sp'first..sp'last+integer32(topdim))
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,topdim);
    begin
      embsys := new Laur_Sys'(ep);
    end;
  end Embed_Nonsquare_System;

  procedure Embed_Nonsquare_System
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk,topdim : in natural32;
                embsys : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Poly_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Poly_Sys(1..integer32(max)) := Square(p);
    ns : natural32;

  begin
    if nbequ > nbunk then
      Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
     -- put("added "); put(nbequ - nbunk,1); put_line(" slack variables");
    end if;
   -- put_line("The squared polynomial system :"); put_line(sp);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Poly_Sys(sp'first..sp'last+integer32(topdim))
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,topdim);
    begin
      embsys := new Poly_Sys'(ep);
    end;
  end Embed_Nonsquare_System;

  procedure Embed_Nonsquare_System
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                nbequ,nbunk,topdim : in natural32;
                embsys : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use DoblDobl_Complex_Laur_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Laur_Sys(1..integer32(max)) := Square(p);
    ns : natural32;

  begin
    if nbequ > nbunk then
      Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
     -- put("added "); put(nbequ - nbunk,1); put_line(" slack variables");
    end if;
   -- put_line("The squared polynomial system :"); put_line(sp);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Laur_Sys(sp'first..sp'last+integer32(topdim))
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,topdim);
    begin
      embsys := new Laur_Sys'(ep);
    end;
  end Embed_Nonsquare_System;

  procedure Embed_Nonsquare_System
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                nbequ,nbunk,topdim : in natural32;
                embsys : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Poly_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Poly_Sys(1..integer32(max)) := Square(p);
    ns : natural32;

  begin
    if nbequ > nbunk then
      Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
     -- put("added "); put(nbequ - nbunk,1); put_line(" slack variables");
    end if;
   -- put_line("The squared polynomial system :"); put_line(sp);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Poly_Sys(sp'first..sp'last+integer32(topdim))
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,topdim);
    begin
      embsys := new Poly_Sys'(ep);
    end;
  end Embed_Nonsquare_System;

  procedure Embed_Nonsquare_System
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                nbequ,nbunk,topdim : in natural32;
                embsys : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use QuadDobl_Complex_Laur_Systems;

    max : constant natural32 := Maximum(nbequ,nbunk);
    sp : constant Laur_Sys(1..integer32(max)) := Square(p);
    ns : natural32;

  begin
    if nbequ > nbunk then
      Witness_Sets_io.Add_Slack_Symbols(nbequ-nbunk);
     -- put("added "); put(nbequ - nbunk,1); put_line(" slack variables");
    end if;
   -- put_line("The squared polynomial system :"); put_line(sp);
    ns := Symbol_Table.Number;
    if ns < nbunk
     then Witness_Sets_io.Add_Extra_Symbols(nbunk-ns);
    end if;
    Witness_Sets_io.Add_Embed_Symbols(topdim);
    declare
      ep : Laur_Sys(sp'first..sp'last+integer32(topdim))
         := Full_Embed_Nonsquare_System(sp,nbequ,nbunk,topdim);
    begin
      embsys := new Laur_Sys'(ep);
    end;
  end Embed_Nonsquare_System;

  procedure Interactive_Square_and_Embed
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                ep : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
                k : out natural32 ) is

    use Standard_Complex_Polynomials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    put("The number of equations : "); put(nq,1); new_line;
    put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Interactive_Embed_Nonsquare_System(file,p,nq,nv,ep,k);
     else Interactive_Embed_Square_System(file,p,ep,k);
    end if;
  end Interactive_Square_and_Embed;

  procedure Interactive_Square_and_Embed
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                ep : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
                k : out natural32 ) is

    use Standard_Complex_Laurentials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    put("The number of equations : "); put(nq,1); new_line;
    put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Interactive_Embed_Nonsquare_System(file,p,nq,nv,ep,k);
     else Interactive_Embed_Square_System(file,p,ep,k);
    end if;
  end Interactive_Square_and_Embed;

  procedure Interactive_Square_and_Embed
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                ep : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                k : out natural32 ) is

    use DoblDobl_Complex_Polynomials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    put("The number of equations : "); put(nq,1); new_line;
    put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Interactive_Embed_Nonsquare_System(file,p,nq,nv,ep,k);
     else Interactive_Embed_Square_System(file,p,ep,k);
    end if;
  end Interactive_Square_and_Embed;

  procedure Interactive_Square_and_Embed
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                ep : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                k : out natural32 ) is

    use DoblDobl_Complex_Laurentials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    put("The number of equations : "); put(nq,1); new_line;
    put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Interactive_Embed_Nonsquare_System(file,p,nq,nv,ep,k);
     else Interactive_Embed_Square_System(file,p,ep,k);
    end if;
  end Interactive_Square_and_Embed;

  procedure Interactive_Square_and_Embed
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                ep : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
                k : out natural32 ) is

    use QuadDobl_Complex_Polynomials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    put("The number of equations : "); put(nq,1); new_line;
    put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Interactive_Embed_Nonsquare_System(file,p,nq,nv,ep,k);
     else Interactive_Embed_Square_System(file,p,ep,k);
    end if;
  end Interactive_Square_and_Embed;

  procedure Interactive_Square_and_Embed
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                ep : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
                k : out natural32 ) is

    use QuadDobl_Complex_Laurentials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
    put("The number of equations : "); put(nq,1); new_line;
    put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Interactive_Embed_Nonsquare_System(file,p,nq,nv,ep,k);
     else Interactive_Embed_Square_System(file,p,ep,k);
    end if;
  end Interactive_Square_and_Embed;

  procedure Square_and_Embed
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                ep : out Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use Standard_Complex_Polynomials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
   -- put("The number of equations : "); put(nq,1); new_line;
   -- put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Embed_Nonsquare_System(p,nq,nv,topdim,ep);
     else Embed_Square_System(p,topdim,ep);
    end if;
  end Square_and_Embed;

  procedure Square_and_Embed
              ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                topdim : in natural32;
                ep : out Standard_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use Standard_Complex_Laurentials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
   -- put("The number of equations : "); put(nq,1); new_line;
   -- put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Embed_Nonsquare_System(p,nq,nv,topdim,ep);
     else Embed_Square_System(p,topdim,ep);
    end if;
  end Square_and_Embed;

  procedure Square_and_Embed
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                ep : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use DoblDobl_Complex_Polynomials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
   -- put("The number of equations : "); put(nq,1); new_line;
   -- put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Embed_Nonsquare_System(p,nq,nv,topdim,ep);
     else Embed_Square_System(p,topdim,ep);
    end if;
  end Square_and_Embed;

  procedure Square_and_Embed
              ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                topdim : in natural32;
                ep : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use DoblDobl_Complex_Laurentials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
   -- put("The number of equations : "); put(nq,1); new_line;
   -- put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Embed_Nonsquare_System(p,nq,nv,topdim,ep);
     else Embed_Square_System(p,topdim,ep);
    end if;
  end Square_and_Embed;

  procedure Square_and_Embed
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                topdim : in natural32;
                ep : out QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    use QuadDobl_Complex_Polynomials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
   -- put("The number of equations : "); put(nq,1); new_line;
   -- put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Embed_Nonsquare_System(p,nq,nv,topdim,ep);
     else Embed_Square_System(p,topdim,ep);
    end if;
  end Square_and_Embed;

  procedure Square_and_Embed
              ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                topdim : in natural32;
                ep : out QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys ) is

    use QuadDobl_Complex_Laurentials;

    nq : constant natural32 := natural32(p'last);
    nv : constant natural32 := Number_of_Unknowns(p(p'first));

  begin
   -- put("The number of equations : "); put(nq,1); new_line;
   -- put("The number of variables : "); put(nv,1); new_line;
    if nq /= nv
     then Embed_Nonsquare_System(p,nq,nv,topdim,ep);
     else Embed_Square_System(p,topdim,ep);
    end if;
  end Square_and_Embed;

  function Remove_Last_Variables
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               n : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Embedding(p(i),n);
    end loop;
    return res;
  end Remove_Last_Variables;

  function Remove_Last_Variables
             ( p : Standard_Complex_Laur_Systems.Laur_Sys;
               n : natural32 )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    res : Standard_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Embedding(p(i),n);
    end loop;
    return res;
  end Remove_Last_Variables;

  function Remove_Last_Variables
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               n : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    res : DoblDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Embedding(p(i),n);
    end loop;
    return res;
  end Remove_Last_Variables;

  function Remove_Last_Variables
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               n : natural32 )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Embedding(p(i),n);
    end loop;
    return res;
  end Remove_Last_Variables;

  function Remove_Last_Variables
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               n : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    res : QuadDobl_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Embedding(p(i),n);
    end loop;
    return res;
  end Remove_Last_Variables;

  function Remove_Last_Variables
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               n : natural32 )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Remove_Embedding(p(i),n);
    end loop;
    return res;
  end Remove_Last_Variables;

  procedure Remove_Last_Variables
              ( p : in out Standard_Complex_Poly_Systems.Poly_Sys;
                n : in natural32 ) is

    use Standard_Complex_Polynomials;
    res : Poly;

  begin
    for i in p'range loop
      res := Remove_Embedding(p(i),n);
      Copy(res,p(i)); Clear(res);
    end loop;
  end Remove_Last_Variables;

  procedure Remove_Last_Variables
              ( p : in out Standard_Complex_Laur_Systems.Laur_Sys;
                n : in natural32 ) is

    use Standard_Complex_Laurentials;
    res : Poly;

  begin
    for i in p'range loop
      res := Remove_Embedding(p(i),n);
      Copy(res,p(i)); Clear(res);
    end loop;
  end Remove_Last_Variables;

  procedure Remove_Last_Variables
              ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys;
                n : in natural32 ) is

    use DoblDobl_Complex_Polynomials;
    res : Poly;

  begin
    for i in p'range loop
      res := Remove_Embedding(p(i),n);
      Copy(res,p(i)); Clear(res);
    end loop;
  end Remove_Last_Variables;

  procedure Remove_Last_Variables
              ( p : in out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                n : in natural32 ) is

    use DoblDobl_Complex_Laurentials;
    res : Poly;

  begin
    for i in p'range loop
      res := Remove_Embedding(p(i),n);
      Copy(res,p(i)); Clear(res);
    end loop;
  end Remove_Last_Variables;

  procedure Remove_Last_Variables
              ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys;
                n : in natural32 ) is

    use QuadDobl_Complex_Polynomials;
    res : Poly;

  begin
    for i in p'range loop
      res := Remove_Embedding(p(i),n);
      Copy(res,p(i)); Clear(res);
    end loop;
  end Remove_Last_Variables;

  procedure Remove_Last_Variables
              ( p : in out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                n : in natural32 ) is

    use QuadDobl_Complex_Laurentials;
    res : Poly;

  begin
    for i in p'range loop
      res := Remove_Embedding(p(i),n);
      Copy(res,p(i)); Clear(res);
    end loop;
  end Remove_Last_Variables;

  function Remove_Embedding
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               dim,ns : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Poly_Systems;
    res : Poly_Sys(p'range);

  begin
    if dim > 0 then
      declare
        wrk : constant Poly_Sys := Remove_Embedding1(p,dim);
        nz : constant natural32 := Number_of_Zero_Equations(wrk);
        res1 : Poly_Sys(1..wrk'last-integer32(nz))
             := wrk(1..wrk'last-integer32(nz));
      begin
        if ns > 0 then
          Remove_Last_Variables(res1,ns);
        end if;
        return res1;
      end;
    elsif ns > 0 then
      res := Remove_Last_Variables(p,ns);
    end if;
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( p : Standard_Complex_Laur_Systems.Laur_Sys;
               dim,ns : natural32 )
             return Standard_Complex_Laur_Systems.Laur_Sys is

    use Standard_Complex_Laur_Systems;
    res : Laur_Sys(p'range);

  begin
    if dim > 0 then
      declare
        wrk : constant Laur_Sys := Remove_Embedding1(p,dim);
        nz : constant natural32 := Number_of_Zero_Equations(wrk);
        res1 : Laur_Sys(1..wrk'last-integer32(nz))
             := wrk(1..wrk'last-integer32(nz));
      begin
        if ns > 0 then
          Remove_Last_Variables(res1,ns);
        end if;
        return res1;
      end;
    elsif ns > 0 then
      res := Remove_Last_Variables(p,ns);
    end if;
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
               dim,ns : natural32 )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;
    res : Poly_Sys(p'range);

  begin
    if dim > 0 then
      declare
        wrk : constant Poly_Sys := Remove_Embedding1(p,dim);
        nz : constant natural32 := Number_of_Zero_Equations(wrk);
        res1 : Poly_Sys(1..wrk'last-integer32(nz))
             := wrk(1..wrk'last-integer32(nz));
      begin
        if ns > 0 then
          Remove_Last_Variables(res1,ns);
        end if;
        return res1;
      end;
    elsif ns > 0 then
      res := Remove_Last_Variables(p,ns);
    end if;
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
               dim,ns : natural32 )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

    use DoblDobl_Complex_Laur_Systems;
    res : Laur_Sys(p'range);

  begin
    if dim > 0 then
      declare
        wrk : constant Laur_Sys := Remove_Embedding1(p,dim);
        nz : constant natural32 := Number_of_Zero_Equations(wrk);
        res1 : Laur_Sys(1..wrk'last-integer32(nz))
             := wrk(1..wrk'last-integer32(nz));
      begin
        if ns > 0 then
          Remove_Last_Variables(res1,ns);
        end if;
        return res1;
      end;
    elsif ns > 0 then
      res := Remove_Last_Variables(p,ns);
    end if;
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
               dim,ns : natural32 )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;
    res : Poly_Sys(p'range);

  begin
    if dim > 0 then
      declare
        wrk : constant Poly_Sys := Remove_Embedding1(p,dim);
        nz : constant natural32 := Number_of_Zero_Equations(wrk);
        res1 : Poly_Sys(1..wrk'last-integer32(nz))
             := wrk(1..wrk'last-integer32(nz));
      begin
        if ns > 0 then
          Remove_Last_Variables(res1,ns);
        end if;
        return res1;
      end;
    elsif ns > 0 then
      res := Remove_Last_Variables(p,ns);
    end if;
    return res;
  end Remove_Embedding;

  function Remove_Embedding
             ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
               dim,ns : natural32 )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys is

    use QuadDobl_Complex_Laur_Systems;
    res : Laur_Sys(p'range);

  begin
    if dim > 0 then
      declare
        wrk : constant Laur_Sys := Remove_Embedding1(p,dim);
        nz : constant natural32 := Number_of_Zero_Equations(wrk);
        res1 : Laur_Sys(1..wrk'last-integer32(nz))
             := wrk(1..wrk'last-integer32(nz));
      begin
        if ns > 0 then
          Remove_Last_Variables(res1,ns);
        end if;
        return res1;
      end;
    elsif ns > 0 then
      res := Remove_Last_Variables(p,ns);
    end if;
    return res;
  end Remove_Embedding;

end Square_and_Embed_Systems;
