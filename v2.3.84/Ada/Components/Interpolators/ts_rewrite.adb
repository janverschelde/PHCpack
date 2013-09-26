with text_io;                            use text_io;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;

procedure ts_rewrite is

-- DESCRIPTION :
--   This is a test program to rewrite polynomials of high degree into
--   polynomials of modest degree, at the expense of extra variables.

  procedure Binary ( k,d : in natural32; deco : out Link_to_Vector ) is

  -- DESCRIPTION :
  --   Writes the binary decomposition of d on standard output.
  --   The step counter is k.  The binary decomposition is in deco,
  --   starting with the least significant bit.
  --   This routine was useful for developing, but is not as user
  --   friendly as the function Binary.

  -- REQUIRED : d > 0;

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

  -- DESCRIPTION :
  --   Returns the length of the binary representation of d,
  --   this is the number of bits plus one.
    
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

  -- DESCRIPTION :
  --   This is a recursive auxiliary procedure to compute the
  --   binary decomposition of the number d.  The k counts the
  --   number of divisions by 2.  The routine should be called
  --   with k = 0.

  -- REQUIRED : d > 0.

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

  -- DESCRIPTION :
  --   Returns the binary decomposition of the natural number d.
  --   The result contains the coefficients of the binary respresentation
  --   of d, starting with the least significant bit.

  begin
    if d = 0 then
      declare
        deco : constant Vector(0..0) := (0..0 => 0);
      begin
        return deco;
      end;
    else
      return Recursive_Binary(0,d);
    end if;
  end Binary;

  procedure Test_Binary ( d : in natural32; deco : in Vector ) is

  -- DESCRIPTION : 
  --   Tests whether the binary decomposition of d matches the value.

    val : natural32 := deco(0);
    acc : natural32 := 1;
    len : constant natural32 := Binary_Length(d);

  begin
    put("Length : "); put(len,1);
    if len = natural32(deco'last)
     then put(" okay ");
     else put(" BUG! ");
    end if;
    for i in 1..deco'last loop
      acc := acc*2;
      if deco(i) = 1
       then val := val + acc;
      end if;
    end loop;
    put(d,1); 
    if d = val
     then put(" = "); put(val,1); put_line("  okay");
     else put(" <> "); put(val,1); put_line("  bug!!!");
    end if;
  end Test_Binary;

  procedure Binary_Decomposition is

  -- DESCRIPTION :
  --   Interactive test on computing binary decompositions of
  --   natural numbers.

    d : natural32 := 0;
    deco : Link_to_Vector;
    ans : character;

  begin
    new_line;
    loop
      put("Give a degree : "); get(d);
      if d > 0 then
        put("The binary decomposition of "); put(d,1);
        put_line(" :");
        Binary(0,d,deco); new_line;
        put("The decomposition vector : "); put(deco.all); new_line;
        Test_Binary(d,deco.all);
      end if;
      declare
        bindeco : constant Vector := Binary(d);
      begin
        put("The decomposition vector : "); put(bindeco); new_line;
        Test_Binary(d,bindeco);
      end;
      put("Do you want more tests ? (y/n) "); get(ans);
      exit when (ans /= 'y');
      Clear(deco);
    end loop;
  end Binary_Decomposition;

  function Rewrite_Term ( n : natural32; t : Term ) return Term is

  -- DESCRIPTION :
  --   Writes the term t as a term in n variables.

  -- REQUIRED : n >= length of binary decomposition of t.dg. 
  --            t is term in one variable.

    res : Term;
    bindg : constant Vector := Binary(t.dg(1));

  begin
   -- put(" bindg : "); put(bindg); new_line;
    res.cf := t.cf;
    res.dg := new Vector'(1..integer32(n) => 0);
    for i in bindg'range loop              -- note : bindg'first = 0
      res.dg(i+1) := bindg(i);             -- therefore we shift
    end loop;
    return res;
  end Rewrite_Term;

  function Rewrite_Poly ( n : natural32; p : Poly ) return Poly is

  -- DESCRIPTION :
  --   Rewrites all terms of the polynomial, n is the number of
  --   variables of the new polynomial on return.

    res : Poly := Null_Poly;

    procedure Rewrite_Term ( t : in Term; continue : out boolean ) is

      rt : Term := Rewrite_Term(n,t);

    begin
      Add(res,rt);
      Clear(rt);
      continue := true;
    end Rewrite_Term;
    procedure Rewrite_Terms is new Visiting_Iterator(Rewrite_Term);

  begin
    Rewrite_Terms(p);
    return res;
  end Rewrite_Poly;

  procedure Lexicon ( sys : in out Poly_Sys; n : in natural32 ) is

  -- DESCRIPTION :
  --   Fills in the equations x(i+1) - x(i)^2 into the system sys,
  --   where n is the number of variables used.  Thus i < n.

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
  end Lexicon;

  function Rewrite_Univariate_Polynomial ( p : Poly ) return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the polynomial system that is equivalent to p.

  -- REQUIRED : p must be a univariate polynomial.
 
    len : constant natural32 := Binary_Length(natural32(Degree(p)));
    res : Poly_Sys(1..integer32(len)+1);

  begin
    Lexicon(res,len+1);
    res(integer32(len)+1) := Rewrite_Poly(len+1,p);
    return res;
  end Rewrite_Univariate_Polynomial;

  procedure Enlarge_Symbol_Table ( n : natural32 ) is

  -- DESCRIPTION :
  --   Enlarges the symbol table with n symbols of the form
  --   z1,..,zn, where z is the first symbol in the table.

    use Symbol_Table;
    sb1 : constant Symbol := Symbol_Table.Get(1);
    ind : integer := sb1'first;
    sb2 : Symbol;

  begin
    while sb1(ind) /= ' ' loop
      ind := ind+1;
    end loop;
    Symbol_Table.Enlarge(n);
    for i in 1..n loop
      sb2 := sb1;
      declare
        order : constant String := Convert(integer32(i));
      begin
        for j in order'range loop
          sb2(ind+j-1) := order(j);
        end loop;
      end;
      Symbol_Table.Add(sb2);
    end loop;
  end Enlarge_Symbol_Table;

  procedure Rewrite_Polynomials is

    n,len : natural32 := 0;
    p : Poly;

  begin
    new_line;
    put_line("Give the number of unknowns, followed by polynomial : ");
    get(n,p);
    put("Your polynomial : "); put(p); new_line;
    put("The symbol used : "); put(Symbol_Table.get(1)); new_line;
    len := Binary_Length(natural32(Degree(p)));
    Enlarge_Symbol_Table(len+1);
    put_line("The equivalent polynomial system : ");
    put(Rewrite_Univariate_Polynomial(p));
  end Rewrite_Polynomials;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Rewriting high degree polynomials into lower degree systems.");
    new_line;
    put_line("MENU for interactive testing : ");
    put_line("  1. test binary decomposition of natural numbers.");
    put_line("  2. rewrite polynomial in one variable.");
    put("Type 1 or 2 to select : "); get(ans);
    case ans is
      when '1' => Binary_Decomposition;
      when '2' => Rewrite_Polynomials;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_rewrite;
