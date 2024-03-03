with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors;
with Symbol_Table;
with Standard_Complex_Poly_Strings;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laur_Strings;
with TripDobl_Polynomial_Convertors;     use TripDobl_Polynomial_Convertors;

package body TripDobl_Complex_Laur_Strings is

-- NOTE : The implementation is a wrapper to Multprec_Complex_Laur_Strings.

  size : constant natural32 := 8;

  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p : in out Poly ) is

    q : Multprec_Complex_Laurentials.Poly;

  begin
    Multprec_Complex_Laur_Strings.Parse(s,k,n,size,q);
    p := Multprec_Laurential_to_TripDobl_Complex(q);
    Multprec_Complex_Laurentials.Clear(q);
  end Parse;

  function Parse ( n : natural32; s : string ) return Poly is

    res : Poly;
    p : integer := s'first;

  begin
    Parse(s,p,n,res);
    return res;
  end Parse;

  function Parse ( n,m : natural32; s : string ) return Laur_Sys is

    res : Laur_Sys(1..integer32(n));
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

  function Parse ( m : natural32; s : Array_of_Strings ) return Laur_Sys is

    res : Laur_Sys(integer32(s'first)..integer32(s'last));
 
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

    q : Multprec_Complex_Laurentials.Poly
      := TripDobl_Complex_to_Multprec_Laurential(p);
    res : constant string := Multprec_Complex_Laur_Strings.Write(q);

  begin
    Multprec_Complex_Laurentials.Clear(q);
    return res;
  end Write;

  function Write ( p : Laur_Sys ) return string is
  begin
    if p'first = p'last
     then return Write(p(p'first));
     else return Write(p(p'first)) & Write(p(p'first+1..p'last));
    end if;
  end Write;

  function Write ( p : Laur_Sys ) return Array_of_Strings is

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

end TripDobl_Complex_Laur_Strings;
