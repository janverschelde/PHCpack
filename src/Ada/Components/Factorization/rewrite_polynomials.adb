with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;

package body Rewrite_Polynomials is

  procedure Binary ( k,d : in natural32; deco : out Link_to_Vector ) is

    n : natural32 := d;
    rest : natural32;

  begin
    if d = 1 then
      put("At step "); put(k,1); put(" : "); put(d,1);
      deco := new Vector(0..integer32(k));
      deco(integer32(k)) := 1;
    else
      if n mod 2 = 1 then
        rest := 1;
        n := n-1;
      else
        rest := 0;
      end if;
      Binary(k+1,n/2,deco);
      put(rest,1);
      deco(integer32(k)) := rest;
    end if;
  end Binary;

  function Binary_Length ( d : natural32 ) return natural32 is

    n : natural32 := d;
    k : natural32 := 0;
 
  begin
    while n > 1 loop
      n := n/2;
      k := k+1;
    end loop;
    return k;
  end Binary_Length;

  function Recursive_Binary ( k,d : natural32 ) return Vector is

    rest : natural32;

  begin
    if d = 1 then
      declare
        deco : Vector(0..integer32(k));
      begin
        deco(integer32(k)) := 1;
        return deco;
      end;
    else
      rest := d mod 2;
      declare
        deco1 : constant Vector := Recursive_Binary(k+1,d/2);
        deco2 : Vector(deco1'range) := deco1;
      begin
        deco2(integer32(k)) := rest;
        return deco2;
      end;
    end if;
  end Recursive_Binary;

  function Binary ( d : natural32 ) return Vector is
  begin
    if d = 0 then
      declare
        deco : Vector(0..0) := (0..0 => 0);
      begin
        return deco;
      end;
    else
      return Recursive_Binary(0,d);
    end if;
  end Binary;

  function Multi_Degrees
             ( p : Poly ) return Standard_Natural_Vectors.Vector is

    nvr : constant integer32 := integer32(Number_of_Unknowns(p));
    res : Standard_Natural_Vectors.Vector(1..nvr);

    procedure Scan_Term ( t : in Term; cont : out boolean ) is
    begin
      for i in t.dg'range loop
        if t.dg(i) > res(i)
         then res(i) := t.dg(i);
        end if;
      end loop;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    res := (res'range => 0);
    Scan_Terms(p);
    return res;
  end Multi_Degrees;

  procedure Number_of_Variables
               ( deg : in Standard_Natural_Vectors.Vector;
                 nvr : out Standard_Natural_Vectors.Vector;
                 tnv : out natural32 ) is
  begin
    tnv := 0;
    for i in deg'range loop
      if deg(i) > 0
       then nvr(i) := Binary_Length(deg(i)) + 1;
            tnv := tnv + nvr(i);
       else nvr(i) := 0;
      end if;
    end loop;
  end Number_of_Variables;

  function Rewrite_Univariate_Term ( n : natural32; t : Term ) return Term is

    res : Term;
    bindg : constant Vector := Binary(t.dg(1));

  begin
    res.cf := t.cf;
    res.dg := new Vector'(1..integer32(n) => 0);
    for i in bindg'range loop              -- note : bindg'first = 0
      res.dg(i+1) := bindg(i);             -- therefore we shift
    end loop;
    return res;
  end Rewrite_Univariate_Term;

  function Rewrite_Multivariate_Term
             ( n : natural32; t : Term;
               nvr : Standard_Natural_Vectors.Vector ) return Term is

    res : Term;
    bindg : constant Vector := Binary(t.dg(1));
    ind : natural32 := 0;

  begin
    res.cf := t.cf;
    res.dg := new Vector'(1..integer32(n) => 0);
    for i in t.dg'range loop
      declare
        bindgi : constant Vector := Binary(t.dg(i));
      begin
        for j in bindgi'range loop               -- note : bindg'first = 0
          res.dg(integer32(ind)+j+1) := bindgi(j);   -- therefore we shift
        end loop;
      end;
      ind := ind + nvr(i);
    end loop;
    return res;
  end Rewrite_Multivariate_Term;

  function Rewrite_Univariate_Poly ( n : natural32; p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Rewrite_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Rewrite_Univariate_Term(n,t);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Rewrite_Term;
    procedure Rewrite_Terms is new Visiting_Iterator(Rewrite_Term);

  begin
    Rewrite_Terms(p);
    return res;
  end Rewrite_Univariate_Poly;

  function Rewrite_Multivariate_Poly
             ( n : natural32; p : Poly;
               nvr : Standard_Natural_Vectors.Vector ) return Poly is

    res : Poly := Null_Poly;

    procedure Rewrite_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Rewrite_Multivariate_Term(n,t,nvr);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Rewrite_Term;
    procedure Rewrite_Terms is new Visiting_Iterator(Rewrite_Term);

  begin
    Rewrite_Terms(p);
    return res;
  end Rewrite_Multivariate_Poly;

  procedure Telescope ( sys : in out Poly_Sys; n : in natural32 ) is

    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Vector'(1..integer32(n) => 0);
    for i in 1..integer32(n)-1 loop
      t.dg(i) := 2;
      sys(i) := Create(t);
      t.dg(i) := 0;
      t.dg(i+1) := 1;
      Sub(sys(i),t);
    end loop;
    Clear(t);
  end Telescope;

  procedure Telescope ( sys : in out Poly_Sys; n : in natural32;
                        nvr : in Standard_Natural_Vectors.Vector ) is

    t : Term;
    indequ,indvar : natural32 := 0;

  begin
    t.cf := Create(1.0);
    t.dg := new Vector'(1..integer32(n) => 0);
    for i in nvr'range loop
      for j in 1..nvr(i)-1 loop
        indequ := indequ + 1;
        indvar := indvar + 1;
        t.dg(integer32(indvar)) := 2;
        sys(integer32(indequ)) := Create(t);
        t.dg(integer32(indvar)) := 0;
        t.dg(integer32(indvar)+1) := 1;
        Sub(sys(integer32(indequ)),t);
        t.dg(integer32(indvar)+1) := 0;
      end loop;
      indvar := indvar + 1;
    end loop;
    Clear(t);
  end Telescope;

  function Rewrite_Univariate_Polynomial ( p : Poly ) return Poly_Sys is

    len : constant natural32 := Binary_Length(natural32(Degree(p)));
    res : Poly_Sys(1..integer32(len)+1);

  begin
    Telescope(res,len+1);
    res(integer32(len)+1) := Rewrite_Univariate_Poly(len+1,p);
    return res;
  end Rewrite_Univariate_Polynomial;

  function Rewrite_Multivariate_Polynomial ( p : Poly ) return Poly_Sys is

    deg : constant Standard_Natural_Vectors.Vector := Multi_Degrees(p);
    nvr : Standard_Natural_Vectors.Vector(deg'range);
    tnv : natural32;

  begin
    Number_of_Variables(deg,nvr,tnv);
    declare
      res : Poly_Sys(1..integer32(tnv));
    begin
      Telescope(res,tnv,nvr);
      res(integer32(tnv)-deg'last+1) := Rewrite_Multivariate_Poly(tnv,p,nvr);
      return res;
    end;
  end Rewrite_Multivariate_Polynomial;

  procedure Enlarge_Symbol_Table ( n : in natural32; sb1 : in Symbol ) is

    use Symbol_Table;
   -- sb1 : Symbol := Symbol_Table.Get(1);
    ind : natural32 := natural32(sb1'first);
    sb2 : Symbol;

  begin
    while sb1(integer(ind)) /= ' ' loop
      ind := ind+1;
    end loop;
    Symbol_Table.Enlarge(n);
    for i in 1..n loop
      sb2 := sb1;
      declare
        order : constant String := Convert(integer64(i));
      begin
        for j in order'range loop
          sb2(integer(ind)+j-1) := order(j);
        end loop;
      end;
      Symbol_Table.Add(sb2);
    end loop;
  end Enlarge_Symbol_Table;

  procedure Define_Symbol_Table
              ( n : in natural32; nvr : in Standard_Natural_Vectors.Vector ) is

    sb : Symbol;

  begin
    Symbol_Table.Init(n);
    sb := (sb'range => ' ');
    sb(1) := 'x';
    for i in nvr'range loop
      declare
        nbi : constant string := Convert(i);
      begin
        for j in nbi'range loop
          sb(1+j) := nbi(j);
        end loop;
        for j in 1..nvr(i) loop
          declare
            nbj : constant string := Convert(integer64(j));
          begin
            for k in nbj'range loop
              sb(1+nbi'last+k) := nbj(k);
            end loop;
            Symbol_Table.Add(sb);
          end;
        end loop;
      end;
    end loop;  
  end Define_Symbol_Table;

end Rewrite_Polynomials;
