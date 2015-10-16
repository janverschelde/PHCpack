with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with Brackets_io;                       use Brackets_io;
with Bracket_Monomials;                 use Bracket_Monomials;
with Symbolic_Minor_Equations;          use Symbolic_Minor_Equations;

package body Remember_Symbolic_Minors is

  function Number_of_Minors ( n,k : natural32 ) return natural32 is
  begin
    return Symbolic_Minor_Equations.Number_of_Maximal_Minors(n,k);
  end Number_of_Minors;

  function Create ( n,k : natural32;
                    x : Standard_Complex_Poly_Matrices.Matrix )
                  return Standard_Symbolic_Minors is

    nq : constant integer32 := integer32(Number_of_Maximal_Minors(n,k));
    bm : Bracket_Monomial := Maximal_Minors(n,k);
    res : Standard_Symbolic_Minors(nq);
    ind : integer32 := 0;

    procedure Expand_Minor ( b : in Bracket; continue : out boolean ) is
    begin
      ind := ind + 1;
      res.b(ind) := new Bracket'(b);
      res.p(ind) := Expanded_Minor(x,b);
      continue := true;
    end Expand_Minor;
    procedure Expand_Minors is new Enumerate_Brackets(Expand_Minor);

  begin
    Expand_Minors(bm);
    Clear(bm);
    return res;
  end Create;

  function Create ( n,k : natural32;
                    x : DoblDobl_Complex_Poly_Matrices.Matrix )
                  return DoblDobl_Symbolic_Minors is

    nq : constant integer32 := integer32(Number_of_Maximal_Minors(n,k));
    bm : Bracket_Monomial := Maximal_Minors(n,k);
    res : DoblDobl_Symbolic_Minors(nq);
    ind : integer32 := 0;

    procedure Expand_Minor ( b : in Bracket; continue : out boolean ) is
    begin
      ind := ind + 1;
      res.b(ind) := new Bracket'(b);
      res.p(ind) := Expanded_Minor(x,b);
      continue := true;
    end Expand_Minor;
    procedure Expand_Minors is new Enumerate_Brackets(Expand_Minor);

  begin
    Expand_Minors(bm);
    Clear(bm);
    return res;
  end Create;

  function Create ( n,k : natural32;
                    x : QuadDobl_Complex_Poly_Matrices.Matrix )
                  return QuadDobl_Symbolic_Minors is

    nq : constant integer32 := integer32(Number_of_Maximal_Minors(n,k));
    bm : Bracket_Monomial := Maximal_Minors(n,k);
    res : QuadDobl_Symbolic_Minors(nq);
    ind : integer32 := 0;

    procedure Expand_Minor ( b : in Bracket; continue : out boolean ) is
    begin
      ind := ind + 1;
      res.b(ind) := new Bracket'(b);
      res.p(ind) := Expanded_Minor(x,b);
      continue := true;
    end Expand_Minor;
    procedure Expand_Minors is new Enumerate_Brackets(Expand_Minor);

  begin
    Expand_Minors(bm);
    Clear(bm);
    return res;
  end Create;

  function Search ( t : Standard_Symbolic_Minors; b : Bracket )
                  return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Polynomials;

  begin
    for i in 1..t.m loop
      if Is_Equal(b,t.b(i).all)
       then return t.p(i);
      end if;
    end loop;
    return Null_Poly;
  end Search;

  function Search ( t : DoblDobl_Symbolic_Minors; b : Bracket )
                  return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Polynomials;

  begin
    for i in 1..t.m loop
      if Is_Equal(b,t.b(i).all)
       then return t.p(i);
      end if;
    end loop;
    return Null_Poly;
  end Search;

  function Search ( t : QuadDobl_Symbolic_Minors; b : Bracket )
                  return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Polynomials;

  begin
    for i in 1..t.m loop
      if Is_Equal(b,t.b(i).all)
       then return t.p(i);
      end if;
    end loop;
    return Null_Poly;
  end Search;

  procedure Write ( t : in Standard_Symbolic_Minors ) is
  begin
    for i in 1..t.m loop
      put(t.b(i).all); put(" : "); put(t.p(i)); new_line;
    end loop;
  end Write;

  procedure Write ( t : in DoblDobl_Symbolic_Minors ) is
  begin
    for i in 1..t.m loop
      put(t.b(i).all); put(" : "); put(t.p(i)); new_line;
    end loop;
  end Write;

  procedure Write ( t : in QuadDobl_Symbolic_Minors ) is
  begin
    for i in 1..t.m loop
      put(t.b(i).all); put(" : "); put(t.p(i)); new_line;
    end loop;
  end Write;

  procedure Query ( t : in Standard_Symbolic_Minors; k : in integer32 ) is

    use Standard_Complex_Polynomials;

    b : Bracket(1..k);
    p : Poly;
    ans : character;

  begin
    loop
      put("Query minor table ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a bracket : "); get(b);
      p := Search(t,b);
      if p = Null_Poly then
        put("There is no minor indexed by "); put(b); new_line;
      else
        put("Minor indexed by "); put(b); put(" is "); put(p); new_line;
      end if;
    end loop;
  end Query;

  procedure Query ( t : in DoblDobl_Symbolic_Minors; k : in integer32 ) is

    use DoblDobl_Complex_Polynomials;

    b : Bracket(1..k);
    p : Poly;
    ans : character;

  begin
    loop
      put("Query minor table ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a bracket : "); get(b);
      p := Search(t,b);
      if p = Null_Poly then
        put("There is no minor indexed by "); put(b); new_line;
      else
        put("Minor indexed by "); put(b); put(" is "); put(p); new_line;
      end if;
    end loop;
  end Query;

  procedure Query ( t : in QuadDobl_Symbolic_Minors; k : in integer32 ) is

    use QuadDobl_Complex_Polynomials;

    b : Bracket(1..k);
    p : Poly;
    ans : character;

  begin
    loop
      put("Query minor table ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a bracket : "); get(b);
      p := Search(t,b);
      if p = Null_Poly then
        put("There is no minor indexed by "); put(b); new_line;
      else
        put("Minor indexed by "); put(b); put(" is "); put(p); new_line;
      end if;
    end loop;
  end Query;

  procedure Clear ( t : in out Standard_Symbolic_Minors ) is
  begin 
    Clear(t.b);
    Standard_Complex_Poly_Systems.Clear(t.p);
  end Clear;

  procedure Clear ( t : in out DoblDobl_Symbolic_Minors ) is
  begin 
    Clear(t.b);
    DoblDobl_Complex_Poly_Systems.Clear(t.p);
  end Clear;

  procedure Clear ( t : in out QuadDobl_Symbolic_Minors ) is
  begin 
    Clear(t.b);
    QuadDobl_Complex_Poly_Systems.Clear(t.p);
  end Clear;

end Remember_Symbolic_Minors;
