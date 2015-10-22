with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Symbol_Table;                       use Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Matrix_Indeterminates;
with Brackets,Brackets_io;               use Brackets,Brackets_io;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Standard_Bracket_Polynomials_io;    use Standard_Bracket_Polynomials_io;
with Straightening_Syzygies;             use Straightening_Syzygies;
with Bracket_Expansions;                 use Bracket_Expansions;

procedure ts_expand is

-- DESCRIPTION :
--   Test on the implementation of bracket expansion.
---  This procedure pre-dates the package SAGBI_Homotopies.

  function Number_of_Brackets ( n,d : natural32 ) return natural32 is

  -- DESCIPTION :
  --   Returns the number of brackets of d entries chosen from n numbers.

    a,b : natural32;

  begin
    a := 1;
    for i in d+1..n loop
      a := a*i;
    end loop;
    b := 1;
    for i in 1..n-d loop
      b := b*i;
    end loop;
    return a/b;
  end Number_of_Brackets;

  procedure Write_Laplace_Expansion
              ( n,d : in natural32; p : in Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Writes the Laplace expansion in expanded form, i.e.: with xij's.

    procedure Write_Term ( t : in Bracket_Term; continue : out boolean ) is

      first : boolean;

      procedure Write_Bracket ( b : in Bracket; cont : out boolean ) is
      begin
        if first then
          if REAL_PART(t.coeff) > 0.0
           then put("+");
           else put("-");
          end if;
          put(b); put("*(");
          first := false;
        else
          put(Expand(n,d,b));
          put(")");
          new_line;
        end if;
        cont := true;
      end Write_Bracket;
      procedure Write_Brackets is new Enumerate_Brackets(Write_Bracket);

    begin
      first := true;
      Write_Brackets(t.monom);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Enumerate_Terms(Write_Term);

  begin
    Write_Terms(p);
  end Write_Laplace_Expansion;

  function Localize ( p : Poly; n,d : natural32 ) return Poly is

  -- DESCRIPTION :
  --   The last (n-d)*d variables are replaced by ones or zeros.

    res,tmp : Poly;
    last : integer32 := integer32(d*n);

  begin
    Copy(p,res);
    for i in 1..d loop
      for j in 1..d loop
        if i = j
         then tmp := Eval(res,Create(1.0),last);
         else tmp := Eval(res,Create(0.0),last);
        end if;
        Copy(tmp,res); Clear(tmp);
        last := last-1;
      end loop;
    end loop;
    return res;
  end Localize;

  procedure Write_Localized_Laplace_Expansion
                ( p : in Bracket_Polynomial; n,d : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the Laplace expansion in expanded form, i.e.: with xij's
  --   and with the last d*d variables set to zeros or ones.

    procedure Write_Term ( t : in Bracket_Term; continue : out boolean ) is

      first : boolean;

      procedure Write_Bracket ( b : in Bracket; cont : out boolean ) is
      begin
        if first then
          if REAL_PART(t.coeff) > 0.0
           then put("+");
           else put("-");
          end if;
          put(b); put("*(");
          first := false;
        else
          put(Localize(Expand(n,d,b),n,d));
          put(")");
          new_line;
        end if;
        cont := true;
      end Write_Bracket;
      procedure Write_Brackets is new Enumerate_Brackets(Write_Bracket);

    begin
      first := true;
      Write_Brackets(t.monom);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Enumerate_Terms(Write_Term);

  begin
    Write_Terms(p);
  end Write_Localized_Laplace_Expansion;

  function Coordinatize ( b : Bracket ) return natural32 is

  -- DESCRIPTION :
  --   Returns the decimal expansion of the bracket.

    res : natural32 := 0;

  begin
    for i in b'range loop
      res := res*10 + b(i);
    end loop;
    return res;
  end Coordinatize;

  function Coordinatize ( p : Bracket_Polynomial ) return Bracket_Polynomial is

  -- DESCRIPTION :
  --   Replaces the first bracket in every monomial by the decimal expansion.

    res : Bracket_Polynomial;

    procedure Coordinatize_Term ( t : in Bracket_Term; cont1 : out boolean ) is

      first,second : boolean;
      bm : Bracket_Monomial;
      bt : Bracket_Term;

      procedure Coordinatize_Bracket ( b : in Bracket; cont2 : out boolean ) is
      begin
        if first then
          bt.coeff := Create(double_float(Coordinatize(b)));
          first := false;
          second := true;
        elsif second then
          bm := Create(b);
        else
          Multiply(bm,b);
        end if;
        cont2 := true;
      end Coordinatize_Bracket;
      procedure Coordinatize_Brackets is
        new Enumerate_Brackets(Coordinatize_Bracket);

    begin
      first := true; second := false;
      Coordinatize_Brackets(t.monom);
      bt.monom := bm;
      if REAL_PART(t.coeff) < 0.0
       then Min(res,bt);
       else Add(res,bt);
      end if;
      cont1 := true;
    end Coordinatize_Term;
    procedure Coordinatize_Terms is new Enumerate_Terms(Coordinatize_Term);

  begin
    Coordinatize_Terms(p);
    return res;
  end Coordinatize;

  function Weight ( e : Standard_Natural_Vectors.Vector; n,d : natural32 )
                  return natural32 is

  -- DESCRIPTION :
  --   Returns the weight of the exponent vector.

    res : natural32 := 0;
    jmp : natural32;
    ind : integer32;

  begin
    for j in 1..d loop
      jmp := d-j;
      for i in 1..n-d loop
        ind := integer32((i-1)*d + j);
        if e(ind) > 0
         then res := res + e(ind)*(i-1)*jmp;
        end if;
      end loop;
    end loop;
    return res;
  end Weight;

  procedure Divide ( p : in out Poly; w : in natural32 ) is

  -- DESCRIPTION :
  --   Divides the polynomial by t^w.

    procedure Divide_Term ( t : in out Term; continue : out boolean ) is
    begin
      t.dg(t.dg'last) := t.dg(t.dg'last) - w;
      continue := true;
    end Divide_Term;
    procedure Divide_Terms is new Changing_Iterator(Divide_Term);

  begin
    Divide_Terms(p);
  end Divide;

  function Lift ( p : Poly; n,d : natural32 ) return Poly is

  -- DESCRIPTION :
  --   Returns the lifted polynomial, which should be an expanded bracket.

    res : Poly := Null_Poly;
    minwei : natural32 := 10000;

    procedure Lift_Term ( t : in Term; continue : out boolean ) is

      tt : Term;
      wei : natural32;

    begin
      tt.cf := t.cf;
      tt.dg := new Standard_Natural_Vectors.Vector(1..integer32(n*d)+1);
      for i in t.dg'range loop
        tt.dg(i) := t.dg(i);
      end loop;
      tt.dg(t.dg'last+1..tt.dg'last) := (t.dg'last+1..tt.dg'last => 0);
      wei := Weight(tt.dg.all,n,d);
      tt.dg(tt.dg'last) := wei;
      Add(res,tt);
      Clear(tt.dg);
      if wei < minwei
       then minwei := wei;
      end if;
      continue := true;
    end Lift_Term;
    procedure Lift_Terms is new Visiting_Iterator(Lift_Term);

  begin
    Lift_Terms(p);
    if minwei /= 0
     then Divide(res,minwei);
    end if;
    return res;
  end Lift;

  procedure Write_Lifted_Localized_Laplace_Expansion
               ( p : in Bracket_Polynomial; n,d : in natural32;
                 L : out Poly ) is

  -- DESCRIPTION :
  --   Writes the Laplace expansion in expanded form, i.e.: with xij's
  --   and with the last d*d variables set to zeros or ones.
  --   Returns the lifted coordinatized polynomial.

    res : Poly := Null_Poly;

    procedure Write_Term ( t : in Bracket_Term; continue : out boolean ) is

      first : boolean;
      cf : integer32;

      procedure Write_Bracket ( b : in Bracket; cont : out boolean ) is

        lp : Poly;

      begin
        if first then
          cf := integer32(Coordinatize(b));
          if REAL_PART(t.coeff) > 0.0
           then put("+");
           else put("-"); cf := -cf;
          end if;
          put(b); put("*(");
          first := false;
        else
          lp := Lift(Localize(Expand(n,d,b),n,d),n,d);
          put(lp);
          put(")");
          new_line;
          Mul(lp,Create(double_float(cf)));
          Add(res,lp);
          Clear(lp);
        end if;
        cont := true;
      end Write_Bracket;
      procedure Write_Brackets is new Enumerate_Brackets(Write_Bracket);

    begin
      first := true;
      Write_Brackets(t.monom);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Enumerate_Terms(Write_Term);

  begin
    Write_Terms(p);
    L := res;
  end Write_Lifted_Localized_Laplace_Expansion;

  procedure Expand_Laplace ( n,d : natural32 ) is

  -- DESCRIPTION :
  --   Writes the expanded bracket polynomial.

    p : constant Bracket_Polynomial := Laplace_Expansion(n,n-d);
   -- q : Bracket_Polynomial; -- := Coordinatize(p);
   -- lq : Poly; --  := Expand(q);
   -- ld : Poly; --  := Localize(lq,n,d);
    lt,l0 : Poly;
    tsb : Symbol_Table.Symbol;

  begin
    put("The Laplace expansion of "); put(n,1); put("*"); put(n,1);
    put("-determinant as product of "); put(d,1); put("- and ");
    put(n-d,1); put_line("-blocks : ");
    put(p);
   -- return;
   -- put_line("The coordinatized Laplace expansion : ");
   -- put(q);
    put_line("Expanded in terms of xij's : ");
    Write_Laplace_Expansion(n,d,p);
   -- put_line("The coordinatized Laplace expansion in terms of xij's : ");
   -- put(lq); new_line;
    put_line("The localized version : ");
    Write_Localized_Laplace_Expansion(p,n,d);
   -- put(ld); new_line;
    Symbol_Table.Enlarge(1);
    tsb(1) := 't';
    for i in 2..tsb'last loop
      tsb(i) := ' ';
    end loop;
    Symbol_Table.Add(tsb);
   -- lt := Lift(ld,n,d);
    put_line("The lifted localized version :");
    Write_Lifted_Localized_Laplace_Expansion(p,n,d,lt);
    put_line("The coordinatized lifted localized polynomial : ");
    put(lt); new_line;
    put_line("The polynomial in the start system : ");
    l0 := Eval(lt,Create(0.0),integer32(n*d+1));
    put(l0); new_line;
  end Expand_Laplace;

  procedure Memory_Consumption ( n,d : natural32 ) is

    nb : natural32 := 0;
    bp : Bracket_Polynomial;

  begin
    put("Give number of expansions : "); get(nb);
    for i in 1..nb loop
      for j in 1..n loop
        bp := Laplace_Expansion(n,j);
        Clear(bp);
      end loop;
    end loop;
  end Memory_Consumption;

  procedure Expand_Brackets ( n,d : in natural32 ) is

  -- DESCRIPTION :
  --   Reads a bracket and writes the bracket expansion.

    ans : character;
    b : Bracket(1..integer32(d));

  begin
    loop
      put("Give "); put(d,1); put(" entries for the bracket : ");
      get(b);
      put("The expansion of the bracket "); put(b); put_line(" :");
      put(Expand(n,d,b)); new_line;
      put("Do you want to expand other brackets ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Expand_Brackets;

  procedure Localized_Expand_Brackets ( n,d : in natural32 ) is

  -- DESCRIPTION :
  --   Reads a bracket and writes the bracket expansion.

    ans : character;
    b : Bracket(1..integer32(d));

  begin
    loop
      put("Give "); put(d,1); put(" entries for the bracket : "); get(b);
      put("The expansion of the bracket "); put(b); put_line(" :");
      put(Localized_Expand(n,d,b)); new_line;
      put("Do you want to expand other brackets ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Localized_Expand_Brackets;

  procedure Main is

    n,d : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing bracket expansions");
    loop
      new_line;
      put("Give number of elements to choose from : "); get(n);
      put("Give the number of entries in bracket : "); get(d);
      Matrix_Indeterminates.Initialize_Symbols(n,d);
      put_line("Choose one of the following :");
      put_line("  1. Expand single brackets.");
      put_line("  2. Expand single brackets in local coordinates.");
      put_line("  3. Apply Laplace expansion.");
      put_line("  4. Test memory consumption.");
      put("Make your choice (1,2,3, or 4) : "); get(ans);
      put("(n,d) = ("); put(n,1); put(","); put(d,1); put(")");
      put("    #brackets : "); put(Number_of_Brackets(n,d),1);
      put("    #equations : "); put((n-d)*d,1); new_line;
      case ans is
        when '1' => Expand_Brackets(n,d);
        when '2' => Localized_Expand_Brackets(n,d);
        when '3' => Expand_Laplace(n,d);
        when '4' => Memory_Consumption(n,d);
        when others => put_line("option not available");
      end case;
      Matrix_Indeterminates.Clear_Symbols;
      put("Do you want more tests for other n and d ? (y/n) "); get(ans);
      exit when ans /= 'y';
    end loop;
  end Main;

begin
  Main;
end ts_expand;
