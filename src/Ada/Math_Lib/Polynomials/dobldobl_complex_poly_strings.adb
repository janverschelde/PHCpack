with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Strings;
with DoblDobl_Complex_Numbers_cv;        use DoblDobl_Complex_Numbers_cv;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Term_Lists;
with Multprec_Complex_Poly_Strings;
with DoblDobl_Polynomial_Convertors;     use DoblDobl_Polynomial_Convertors;

package body DoblDobl_Complex_Poly_Strings is

-- NOTE : The implementation is a wrapper to Multprec_Complex_Poly_Strings.

  size : constant natural32 := 5;

  function Multprec_Terms_to_DoblDobl_Complex
             ( p : Multprec_Complex_Term_Lists.Term_List )
             return DoblDobl_Complex_Term_Lists.Term_List is

    res,res_last : DoblDobl_Complex_Term_Lists.Term_List;
    tmp : Multprec_Complex_Term_Lists.Term_List := p;
    mpt : Multprec_Complex_Polynomials.Term;

  begin
    while not Multprec_Complex_Term_Lists.Is_Null(tmp) loop
      mpt := Multprec_Complex_Term_Lists.Head_Of(tmp);
      declare
        ddt : DoblDobl_Complex_Polynomials.Term;
      begin
        ddt.cf := Multprec_to_DoblDobl_Complex(mpt.cf);
        ddt.dg := DoblDobl_Complex_Polynomials.Degrees(mpt.dg);
        -- the Append makes a copy anyway
        DoblDobl_Complex_Term_Lists.Append(res,res_last,ddt);
      end;
      tmp := Multprec_Complex_Term_Lists.Tail_Of(tmp);
    end loop;
    return res;
  end Multprec_Terms_to_DoblDobl_Complex;

  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p : in out Poly ) is

    q : Multprec_Complex_Polynomials.Poly;

  begin
    Multprec_Complex_Poly_Strings.Parse(s,k,n,size,q);
    p := Multprec_Polynomial_to_DoblDobl_Complex(q);
    Multprec_Complex_Polynomials.Clear(q);
  end Parse;

  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p,p_last : in out Term_List ) is

    q,q_last : Multprec_Complex_Term_Lists.Term_List;
    first,second : DoblDobl_Complex_Term_Lists.Term_List;

  begin
    Multprec_Complex_Poly_Strings.Parse(s,k,n,size,q,q_last);
    p := Multprec_Terms_to_DoblDobl_Complex(q);
    Multprec_Complex_Term_Lists.Clear(q);
    if DoblDobl_Complex_Term_Lists.Is_Null(p) then
      p_last := p;
    else
      first := p;
      second := Tail_Of(first);
      while not Is_Null(second) loop
        first := second;
        second := Tail_Of(second);
      end loop;
      p_last := first;
    end if;
  end Parse;

  function Parse ( n : natural32; s : string ) return Poly is

    res : Poly;
    p : integer := s'first;

  begin
    Parse(s,p,n,res);
    return res;
  end Parse;

  function Parse ( n : natural32; s : string ) return Term_List is

    res,res_last : Term_List;
    p : integer := s'first;

  begin
    Parse(s,p,n,res,res_last);
    return res;
  end Parse;

  function Parse ( n,m : natural32; s : string ) return Poly_Sys is

    res : Poly_Sys(1..integer32(n));
    ind : constant Standard_Natural_Vectors.Vector(1..integer32(n))
        := Standard_Complex_Poly_Strings.Delimiters(n,s);

  begin
   -- This function may be called when the symbol table has not yet
   -- been initialized and then it should not crash.
    if Symbol_Table.Number < m then
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    res(1) := Parse(m,s(s'first..integer(ind(1))));
    for i in 2..integer32(n) loop
      res(i) := Parse(m,s(integer(ind(i-1)+1)..integer(ind(i))));
    end loop;
    return res;
  end Parse;

  function Parse ( n,m : natural32; s : string ) return Array_of_Term_Lists is

    res : Array_of_Term_Lists(1..integer32(n));
    ind : constant Standard_Natural_Vectors.Vector(1..integer32(n))
        := Standard_Complex_Poly_Strings.Delimiters(n,s);

  begin
    if Symbol_Table.Number < m then -- same comment as other Parse
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    res(1) := Parse(m,s(s'first..integer(ind(1))));
    for i in 2..integer32(n) loop
      res(i) := Parse(m,s(integer(ind(i-1)+1)..integer(ind(i))));
    end loop;
    return res;
  end Parse;

  function Parse ( m : natural32; s : Array_of_Strings ) return Poly_Sys is

    res : Poly_Sys(integer32(s'first)..integer32(s'last));
 
  begin
    if Symbol_Table.Number < m then -- same comment as other Parse
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    for i in s'range loop
      declare
      begin
        res(integer32(i)) := Parse(m,s(i).all);
      exception
        when others => put("something is wrong with string ");
                       put(natural32(i),1);
                       new_line; put_line(s(i).all); raise;
      end;
    end loop;
    return res;
  end Parse;

  function Parse ( m : natural32; s : Array_of_Strings ) 
                 return Array_of_Term_Lists is

    res : Array_of_Term_Lists(integer32(s'first)..integer32(s'last));
 
  begin
    for i in s'range loop
      declare
      begin
        res(integer32(i)) := Parse(m,s(i).all);
      exception
        when others => put("something is wrong with string ");
                       put(natural32(i),1);
                       new_line; put_line(s(i).all); raise;
      end;
    end loop;
    return res;
  end Parse;

  function Size_Limit ( p : Poly ) return natural32 is

    nbtrm : constant natural64 := natural64(Number_of_Terms(p));
    nbvar : constant natural64 := natural64(Number_of_Unknowns(p));
    symsz : constant natural64 := 5;
    cffsz : constant natural64 := 80;
    bound : constant natural64 := 2**31 - 1;
    res : constant natural64 := nbtrm*nbvar*symsz*cffsz;

  begin
    if res > bound
     then return natural32(bound);
     else return natural32(res);
    end if;
  end Size_Limit;

  function Write ( p : Poly ) return string is

    q : Multprec_Complex_Polynomials.Poly
      := DoblDobl_Complex_to_Multprec_Polynomial(p);
    res : constant string := Multprec_Complex_Poly_Strings.Write(q);

  begin
    Multprec_Complex_Polynomials.Clear(q);
    return res;
  end Write;

  function Write ( p : Poly; s : Array_of_Symbols ) return string is

    q : Multprec_Complex_Polynomials.Poly
      := DoblDobl_Complex_to_Multprec_Polynomial(p);
    res : constant string := Multprec_Complex_Poly_Strings.Write(q,s);

  begin
    Multprec_Complex_Polynomials.Clear(q);
    return res;
  end Write;

  function Write ( p : Poly_Sys ) return string is
  begin
    if p'first = p'last
     then return Write(p(p'first));
     else return Write(p(p'first)) & Write(p(p'first+1..p'last));
    end if;
  end Write;

  function Write ( p : Poly_Sys; s : Array_of_Symbols ) return string is
  begin
    if p'first = p'last
     then return Write(p(p'first),s);
     else return Write(p(p'first),s) & Write(p(p'first+1..p'last),s);
    end if;
  end Write;

  function Write ( p : Poly_Sys ) return Array_of_Strings is

    res : Array_of_Strings(integer(p'first)..integer(p'last));

  begin
    for i in res'range loop
      declare
        s : constant string := Write(p(integer32(i)));
      begin
        res(i) := new string'(s);
      end;
    end loop;
    return res;
  end Write;

  function Write ( p : Poly_Sys; s : Array_of_Symbols )
                 return Array_of_Strings is

    res : Array_of_Strings(integer(p'first)..integer(p'last));

  begin
    for i in res'range loop
      declare
        pstr : constant string := Write(p(integer32(i)),s);
      begin
        res(i) := new string'(pstr);
      end;
    end loop;
    return res;
  end Write;

end DoblDobl_Complex_Poly_Strings;
