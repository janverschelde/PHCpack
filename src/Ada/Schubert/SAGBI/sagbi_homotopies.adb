with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Brackets;                           use Brackets;
with Bracket_Monomials;                  use Bracket_Monomials;
with Standard_Bracket_Polynomials;       use Standard_Bracket_Polynomials;
with Straightening_Syzygies;             use Straightening_Syzygies;
with Bracket_Expansions;                 use Bracket_Expansions;
with Evaluated_Minors;                   use Evaluated_Minors;

package body SAGBI_Homotopies is

  function Coordinatize_Hexadecimal ( b : Bracket ) return natural32 is

  -- DESCRIPTION :
  --   Returns the hexadecimal expansion of the entries in the bracket.

    res : natural32 := 0;

  begin
    for i in b'range loop
      res := res*16 + b(i);
    end loop;
    return res;
  end Coordinatize_Hexadecimal;

  function Unsigned ( i : integer32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the unsigned integer.

    n : natural32;

  begin
    if i < 0
     then n := natural32(-i);
     else n := natural32(i);
    end if;
    return n;
  end Unsigned;

  function Bracketize_Hexadecimal ( n,d : natural32 ) return Bracket is

  -- DESCRIPTION :
  --   Returns the d-bracket from the hexadecimal expansion n.

    res : Bracket(1..integer32(d));
    nn : natural32 := n;

  begin
    for i in reverse 1..d loop
      res(integer32(i)) := nn mod 16;
      nn := nn/16;
    end loop;
    return res;
  end Bracketize_Hexadecimal;

--  function Coordinatize 
--             ( p : Bracket_Polynomial ) return Bracket_Polynomial is

  -- DESCRIPTION :
  --   Replaces the first bracket in every monomial by the decimal expansion.

--    res : Bracket_Polynomial;
--
--    procedure Coordinatize_Term
--                  ( t : in Bracket_Term; cont1 : out boolean ) is
--
--      first,second : boolean;
--      bm : Bracket_Monomial;
--      bt : Bracket_Term;
--
--      procedure Coordinatize_Bracket
--                   ( b : in Bracket; cont2 : out boolean ) is
--      begin
--        if first then
--          bt.coeff := Create(double_float(Coordinatize_Hexadecimal(b)));
--          first := false;
--          second := true;
--        elsif second then
--          bm := Create(b);
--        else
--          Multiply(bm,b);
--        end if;
--        cont2 := true;
--      end Coordinatize_Bracket;
--      procedure Coordinatize_Brackets is
--        new Enumerate_Brackets(Coordinatize_Bracket);
--
--    begin
--      first := true; second := false;
--      Coordinatize_Brackets(t.monom);
--      bt.monom := bm;
--      if REAL_PART(t.coeff) < 0.0
--       then Min(res,bt);
--       else Add(res,bt);
--      end if;
--      cont1 := true;
--    end Coordinatize_Term;
--    procedure Coordinatize_Terms is new Enumerate_Terms(Coordinatize_Term);
--
--  begin
--    Coordinatize_Terms(p);
--    return res;
--  end Coordinatize;

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

  function Weight ( e : Standard_Natural_Vectors.Vector; n,d : natural32 )
                  return natural32 is

  -- DESCRIPTION :
  --   Returns the weight of the exponent vector for the localization that
  --   takes the d-by-d identitity matrix in the lower-right of the d-plane.
  --   The lifting recipe is xij*t^((i-1)*(d-j)).

    res : natural32 := 0;
    jmp : natural32;
    ind : integer32;

  begin
    for j in 1..integer32(d) loop
      jmp := d-natural32(j);
      for i in 1..integer32(n-d) loop
        ind := (i-1)*integer32(d) + j;
        if e(ind) > 0
         then res := res + natural32((i-1))*jmp;
        end if;
      end loop;
    end loop;
    return res;
  end Weight;

  function Weight ( locmap : Standard_Natural_Matrices.Matrix;
                    e : Standard_Natural_Vectors.Vector ) return natural32 is

  -- DESCRIPTION :
  --   Returns the weight of the exponent vector as xij*t^((i-1)*(d-j))
  --   for the localization pattern in locmap.

    res : natural32 := 0;
    d : constant natural32 := natural32(locmap'length(2));
    jmp : natural32;
    ind : integer32;

  begin
    ind := 0;
    for i in locmap'range(1) loop
      for j in locmap'range(2) loop
        jmp := d-natural32(j);
        if locmap(i,j) = 2 then
          ind := ind+1;
          if e(ind) > 0
           then res := res + (natural32(i)-1)*jmp;
          end if;
        end if;
      end loop;
    end loop;
    return res;
  end Weight;

  function Lift ( p : Poly; n,d : natural32 ) return Poly is

  -- DESCRIPTION :
  --   Returns the lifted polynomial, where the xij is lifted according
  --   to xij*t^((i-1)*(d-j)).  The lowest powers of t are divided out.
  --   The d-by-d identity matrix is the lower-right of the d-plane.

    res : Poly := Null_Poly;
    first : boolean := true;
    minwei : natural32;

    procedure Lift_Term ( t : in Term; continue : out boolean ) is

      tt : Term;
      wei : natural32;

    begin
      tt.cf := t.cf;
      tt.dg := new Standard_Natural_Vectors.Vector(1..t.dg'last+1);
      tt.dg(t.dg'range) := t.dg.all;
      wei := Weight(t.dg.all,n,d);
      tt.dg(tt.dg'last) := wei;
      Add(res,tt);
      Clear(tt.dg);
      if first then
        minwei := wei;
        first := false;
      elsif wei < minwei then
        minwei := wei;
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

  function Lift ( locmap : Standard_Natural_Matrices.Matrix; p : Poly )
                return Poly is

  -- DESCRIPTION :
  --   Lifts p as to xij*t^((i-1)*(d-j)) and divides by the lowest powers
  --   of t, respecting the localization pattern in locmap.

    res : Poly := Null_Poly;
    first : boolean := true;
    minwei : natural32;

    procedure Lift_Term ( t : in Term; continue : out boolean ) is

      tt : Term;
      wei : natural32;

    begin
      tt.cf := t.cf;
      tt.dg := new Standard_Natural_Vectors.Vector(1..t.dg'last+1);
      tt.dg(t.dg'range) := t.dg.all;
      wei := Weight(locmap,t.dg.all);
      tt.dg(tt.dg'last) := wei;
      Add(res,tt);
      Clear(tt.dg);
      if first then
        minwei := wei;
        first := false;
      elsif wei < minwei then
        minwei := wei;
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

-- TARGET ROUTINES :

  function Lifted_Localized_Laplace_Expansion
             ( n,d : natural32 ) return Poly is

    res : Poly := Null_Poly;
    p : constant Bracket_Polynomial := Laplace_Expansion(n,n-d);

    procedure Visit_Term ( t : in Bracket_Term; continue : out boolean ) is

      first : boolean := true;
      cf : integer32;

      procedure Visit_Bracket ( b : in Bracket; cont : out boolean ) is

        pb,lp : Poly;

      begin
        if first then
          cf := integer32(Coordinatize_Hexadecimal(b));
          first := false;
        else
          pb := Localized_Expand(n,d,b);
          lp := Lift(pb,n,d);  Clear(pb);
          Mul(lp,Create(double_float(cf)));
          Add(res,lp);
          Clear(lp);
        end if;
        cont := true;
      end Visit_Bracket;
      procedure Visit_Brackets is new Enumerate_Brackets(Visit_Bracket);

    begin
      Visit_Brackets(t.monom);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Enumerate_Terms(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Lifted_Localized_Laplace_Expansion;

  function Lifted_Localized_Laplace_Expansion
             ( locmap : Standard_Natural_Matrices.Matrix ) return Poly is

    res : Poly := Null_Poly;
    n : constant natural32 := natural32(locmap'length(1));
    d : constant natural32 := natural32(locmap'length(2));
    p : constant Bracket_Polynomial := Laplace_Expansion(n,n-d);

    procedure Visit_Term ( t : in Bracket_Term; continue : out boolean ) is

      first : boolean := true;
      cf : integer32;

      procedure Visit_Bracket ( b : in Bracket; cont : out boolean ) is

        pb,lp : Poly;

      begin
        if first then
          cf := integer32(Coordinatize_Hexadecimal(b));
          first := false;
        else
          pb := Expand(locmap,b);
          Reduce_Variables(locmap,pb);
          lp := Lift(locmap,pb);  Clear(pb);
          Mul(lp,Create(double_float(cf)));
          Add(res,lp);
          Clear(lp);
        end if;
        cont := true;
      end Visit_Bracket;
      procedure Visit_Brackets is new Enumerate_Brackets(Visit_Bracket);

    begin
      Visit_Brackets(t.monom);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Enumerate_Terms(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Lifted_Localized_Laplace_Expansion;

  function Intersection_Coefficients
              ( m : Standard_Floating_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);
    nmd : constant natural32 := natural32(m'last(2));
    ind : integer32;
    b : Bracket(1..integer32(nmd));

  begin
    for i in c'range loop
      ind := integer32(REAL_PART(c(i)));
      b := Bracketize_Hexadecimal(Unsigned(ind),nmd);
      if ind > 0
       then res(i) := Create(Determinant(m,b));
       else res(i) := Create(-Determinant(m,b));
      end if;
    end loop;
    return res;
  end Intersection_Coefficients;

  function Intersection_Coefficients
              ( m : Standard_Complex_Matrices.Matrix;
                c : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);
    nmd : constant natural32 := natural32(m'last(2));
    ind : integer32;
    b : Bracket(1..integer32(nmd));

  begin
    for i in c'range loop
      ind := integer32(REAL_PART(c(i)));
      b := Bracketize_Hexadecimal(Unsigned(ind),nmd);
      if ind > 0
       then res(i) := Determinant(m,b);
       else res(i) := -Determinant(m,b);
      end if;
    end loop;
    return res;
  end Intersection_Coefficients;

  function Intersection_Condition
             ( m : Standard_Floating_Matrices.Matrix; p : Poly ) return Poly is

    res : Poly := Null_Poly;
    nmd : constant natural32 := natural32(m'last(2));

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      c : constant integer32 := integer32(REAL_PART(t.cf));
      b : constant Bracket(1..integer32(nmd))
        := Bracketize_Hexadecimal(Unsigned(c),nmd);
      det : constant double_float := Determinant(m,b);
      rt : Term;

    begin
      if c > 0
       then rt.cf := Create(det);
       else rt.cf := Create(-det);
      end if;
      rt.dg := t.dg;
      Add(res,rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Intersection_Condition;

  function Intersection_Condition
             ( m : Standard_Complex_Matrices.Matrix; p : Poly ) return Poly is

    res : Poly := Null_Poly;
    nmd : constant natural32 := natural32(m'last(2));

    procedure Visit_Term ( t : in Term; continue : out boolean ) is

      c : constant integer32 := integer32(REAL_PART(t.cf));
      b : constant Bracket(1..integer32(nmd))
        := Bracketize_Hexadecimal(Unsigned(c),nmd);
      det : constant Complex_Number := Determinant(m,b);
      rt : Term;

    begin
      if c > 0
       then rt.cf :=  det;
       else rt.cf := -det;
      end if;
      rt.dg := t.dg;
      Add(res,rt);
      continue := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Intersection_Condition;

end SAGBI_Homotopies;
