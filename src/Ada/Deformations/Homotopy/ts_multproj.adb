with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Random_Numbers;
with Numbers_io;
with Characters_and_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;

procedure ts_multproj is

-- DESCRIPTION :
--   Development of multiprojective coordinate transformations.

  function Make_Partition
             ( n,m : natural32; p : Standard_Natural_Vectors.Vector )
             return partition is

  -- DESCRIPTION :
  --   Returns the m-partition p, represented as a partition of
  --   a set of unknowns, a subset of the total set of n variables.

    res : partition(1..m);

  begin
    for i in res'range loop
      res(i) := Sets_of_Unknowns.Create(n); -- initialize each set
    end loop;
    for i in p'range loop
      Sets_of_Unknowns.Add(res(p(i)),natural32(i));
    end loop;
    return res;
  end Make_Partition;

  function Multiset_Degrees
             ( p : in Standard_Complex_Polynomials.Poly;
               m : in natural32; mp : in Partition )
             return Standard_Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the degrees of the polynomial p in the m sets
  --   of the partition mp, as vector of range 1..m.

    res : Standard_Integer_Vectors.Vector(1..integer32(m));

  begin
    for i in mp'range loop
      res(integer32(i)) := Degrees_in_Sets_of_Unknowns.Degree(p,mp(i));
    end loop;
    return res;
  end Multiset_Degrees;

  function Make_Homogeneous
             ( t : Standard_Complex_Polynomials.Term; 
               d : Standard_Integer_Vectors.Vector;
               m : natural32; p : Partition )
             return Standard_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Returns the term with m variables added,
  --   with degrees to match the given multiset degrees in d.

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
      deg := Degrees_in_Sets_of_Unknowns.Degree(t,p(natural32(i)));
      res.dg(lst+i) := natural32(d(i) - deg);
    end loop;
    return res;
  end Make_Homogeneous;

  function Make_Homogeneous
             ( p : in Standard_Complex_Polynomials.Poly; 
               m : in natural32; z : in Partition )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the polynomial p with m variables added,
  --   with degrees to match the given multiset degrees in d.

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
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys; 
               m : in natural32; z : in Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns the polynomial system p with m variables added,
  --   with degrees to match the multiset degrees of each polynomial.

    res : Standard_Complex_Poly_Systems.Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Make_Homogeneous(p(i),m,z);
    end loop;
    return res;
  end Make_Homogeneous;

  function Standard_Random_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Returns a term in the i-th variable, with random coefficient,
  --   as a term in n variables.

    res : Standard_Complex_Polynomials.Term;

  begin
    res.cf := Standard_Random_Numbers.Random1;
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end Standard_Random_Linear_Term;

  function Standard_Start_Linear_Term
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Term is

  -- DESCRIPTION :
  --   Returns the i-th variable as a term in n variables,
  --   with coefficient equal to one.

    res : Standard_Complex_Polynomials.Term;

  begin
    res.cf := Standard_Complex_Numbers.Create(1.0);
    res.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    res.dg(integer32(i)) := 1;
    return res;
  end Standard_Start_Linear_Term;

  function Standard_Random_Linear_Polynomial
             ( n : natural32; s : Sets_of_Unknowns.Set )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns a linear polynomial in the variables in s,
  --   with nonzero constant coefficient as a polynomial in n variables.

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

  function Standard_Start_Linear_Polynomial
             ( n,i : natural32 )
             return Standard_Complex_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns the start polynomial Zi - 1.

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

  function Standard_Random_Linear_Polynomials
             ( n,m : natural32; z : Partition )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns m random linear polynomials in n+m variables in the sets
  --   of the partition z, with the constant added and the random term
  --   for the extra homogeneous m-th variable.

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

  function Standard_Start_Linear_Polynomials
             ( n,m : natural32 )
             return Standard_Complex_Poly_Systems.Poly_Sys is

  -- DESCRIPTION :
  --   Returns m start polynomials in n+m variables Zi - 1,
  --   for i in range 1..m.

    dim : constant natural32 := n+m;
    res : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m));
    
  begin
    for i in 1..m loop
      res(integer32(i)) := Standard_Start_Linear_Polynomial(dim,n+i);
    end loop;
    return res;
  end Standard_Start_Linear_Polynomials;

  function Interactive_Partition
             ( m : natural32 ) return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a partition of the symbols into m sets,
  --   prompting the user for each symbol in the table for its set.
  --   On return is a vector of range 1..n, where n is the number of symbols.

    nbr : constant integer32 := integer32(Symbol_Table.Number);
    res : Standard_Natural_Vectors.Vector(1..nbr) := (1..nbr => 0);
    idx : natural32 := 0;

  begin
    for i in 1..nbr loop
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(i));
      begin
        loop
          put("-> which set will "); Symbol_Table_io.put(sb);
          put(" be in? ");
          Numbers_io.Read_Natural(idx);
          exit when ((idx > 0) and (idx <= m));
          put("-> index must be larger than 0 and less than "); 
          put(m+1,1); put_line(".  Please try again.");
        end loop;
        res(i) := idx;
      end;
    end loop;
    return res;
  end Interactive_Partition;

  procedure Test ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                   m : in natural32 ) is

    mpart : constant Standard_Natural_Vectors.Vector
          := Interactive_Partition(m);
    nbr : constant natural32 := Symbol_Table.Number;
    spart : constant Partition := Make_Partition(nbr,m,mpart);
    deg : Standard_Integer_Vectors.Vector(1..integer32(m));
    mhp : Standard_Complex_Poly_Systems.Poly_Sys(p'range);
    lhp : constant Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := Standard_Random_Linear_Polynomials(nbr,m,spart);
    shp : constant Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(m))
        := Standard_Start_Linear_Polynomials(nbr,m);

  begin
    put("The partition :"); put(mpart); new_line;
    put("Its symbolic form : "); put(spart); new_line;
    for i in p'range loop
      deg := Multiset_Degrees(p(i),m,spart);
      put("degrees of polynomial "); put(i,1);
      put(" :"); put(deg); new_line;
    end loop;
    Symbol_Table.Enlarge(m);
    for i in 1..m loop
      declare
        sv : constant string := "Z" & Characters_and_Numbers.nConvert(i);
      begin
        Symbol_Table.Add_String(sv);
      end;
    end loop;
    mhp := Make_Homogeneous(p,m,spart);
    put("The "); put(m,1); put_line("-homogeneous form :"); put(mhp);
    put_line("The linear equations :"); put_line(lhp);
    put_line("The start equations :"); put(shp);
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for a polynomial system
  --   and then interactively defines a partition for the set of unknowns.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    m : natural32 := 0;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    put("-> read "); put(lp'last,1); put_line(" polynomials.");
    new_line;
    put("Give the number of sets in the partition : "); get(m);
    Test(lp.all,m);
  end Main;

begin
  Main;
end ts_multproj;
