with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Brackets_io;                       use Brackets_io;
with Bracket_Monomials;                 use Bracket_Monomials;
with Symbolic_Minor_Equations;          use Symbolic_Minor_Equations;
with Evaluated_Minors;                  use Evaluated_Minors;

package body Remember_Numeric_Minors is

  function Number_of_Minors ( n,k : natural32 ) return natural32 is
  begin
    return Symbolic_Minor_Equations.Number_of_Maximal_Minors(n,k);
  end Number_of_Minors;

  function Create ( n,k : natural32;
                    x : Standard_Complex_Matrices.Matrix )
                  return Standard_Numeric_Minors is

    nq : constant integer32 := integer32(Number_of_Maximal_Minors(n,k));
    bm : Bracket_Monomial := Maximal_Minors(n,k);
    res : Standard_Numeric_Minors(nq);
    ind : integer32 := 0;

    procedure Expand_Minor ( b : in Bracket; continue : out boolean ) is
    begin
      ind := ind + 1;
      res.b(ind) := new Bracket'(b);
      res.v(ind) := Determinant(x,b);
      continue := true;
    end Expand_Minor;
    procedure Expand_Minors is new Enumerate_Brackets(Expand_Minor);

  begin
    Expand_Minors(bm);
    Clear(bm);
    return res;
  end Create;

  function Create ( n,k : natural32;
                    x : DoblDobl_Complex_Matrices.Matrix )
                  return DoblDobl_Numeric_Minors is

    nq : constant integer32 := integer32(Number_of_Maximal_Minors(n,k));
    bm : Bracket_Monomial := Maximal_Minors(n,k);
    res : DoblDobl_Numeric_Minors(nq);
    ind : integer32 := 0;

    procedure Expand_Minor ( b : in Bracket; continue : out boolean ) is
    begin
      ind := ind + 1;
      res.b(ind) := new Bracket'(b);
      res.v(ind) := Determinant(x,b);
      continue := true;
    end Expand_Minor;
    procedure Expand_Minors is new Enumerate_Brackets(Expand_Minor);

  begin
    Expand_Minors(bm);
    Clear(bm);
    return res;
  end Create;

  function Create ( n,k : natural32;
                    x : QuadDobl_Complex_Matrices.Matrix )
                  return QuadDobl_Numeric_Minors is

    nq : constant integer32 := integer32(Number_of_Maximal_Minors(n,k));
    bm : Bracket_Monomial := Maximal_Minors(n,k);
    res : QuadDobl_Numeric_Minors(nq);
    ind : integer32 := 0;

    procedure Expand_Minor ( b : in Bracket; continue : out boolean ) is
    begin
      ind := ind + 1;
      res.b(ind) := new Bracket'(b);
      res.v(ind) := Determinant(x,b);
      continue := true;
    end Expand_Minor;
    procedure Expand_Minors is new Enumerate_Brackets(Expand_Minor);

  begin
    Expand_Minors(bm);
    Clear(bm);
    return res;
  end Create;

  function Search ( t : Standard_Numeric_Minors; b : Bracket )
                  return Standard_Complex_Numbers.Complex_Number is
  begin
    for i in 1..t.m loop
      if Is_Equal(b,t.b(i).all)
       then return t.v(i);
      end if;
    end loop;
    return Standard_Complex_Numbers.Create(0.0);
  end Search;

  function Search ( t : DoblDobl_Numeric_Minors; b : Bracket )
                  return DoblDobl_Complex_Numbers.Complex_Number is
  begin
    for i in 1..t.m loop
      if Is_Equal(b,t.b(i).all)
       then return t.v(i);
      end if;
    end loop;
    return DoblDobl_Complex_Numbers.Create(integer(0));
  end Search;

  function Search ( t : QuadDobl_Numeric_Minors; b : Bracket )
                  return QuadDobl_Complex_Numbers.Complex_Number is
  begin
    for i in 1..t.m loop
      if Is_Equal(b,t.b(i).all)
       then return t.v(i);
      end if;
    end loop;
    return QuadDobl_Complex_Numbers.Create(integer(0));
  end Search;

  procedure Write ( t : in Standard_Numeric_Minors ) is
  begin
    for i in 1..t.m loop
      put(t.b(i).all); put(" : "); put(t.v(i)); new_line;
    end loop;
  end Write;

  procedure Write ( t : in DoblDobl_Numeric_Minors ) is
  begin
    for i in 1..t.m loop
      put(t.b(i).all); put(" : "); put(t.v(i)); new_line;
    end loop;
  end Write;

  procedure Write ( t : in QuadDobl_Numeric_Minors ) is
  begin
    for i in 1..t.m loop
      put(t.b(i).all); put(" : "); put(t.v(i)); new_line;
    end loop;
  end Write;

  procedure Query ( t : in Standard_Numeric_Minors; k : in integer32 ) is

    b : Bracket(1..k);
    v : Standard_Complex_Numbers.Complex_Number;
    ans : character;

  begin
    loop
      put("Query minor table ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a bracket : "); get(b);
      v := Search(t,b);
      put("Value of minor "); put(b);
      put(" is "); put(v); new_line;
    end loop;
  end Query;

  procedure Query ( t : in DoblDobl_Numeric_Minors; k : in integer32 ) is

    b : Bracket(1..k);
    v : DoblDobl_Complex_Numbers.Complex_Number;
    ans : character;

  begin
    loop
      put("Query minor table ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a bracket : "); get(b);
      v := Search(t,b);
      put("Value of minor "); put(b);
      put(" is "); put(v); new_line;
    end loop;
  end Query;

  procedure Query ( t : in QuadDobl_Numeric_Minors; k : in integer32 ) is

    b : Bracket(1..k);
    v : QuadDobl_Complex_Numbers.Complex_Number;
    ans : character;

  begin
    loop
      put("Query minor table ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give a bracket : "); get(b);
      v := Search(t,b);
      put("Value of minor "); put(b);
      put(" is "); put(v); new_line;
    end loop;
  end Query;

  procedure Clear ( t : in out Standard_Numeric_Minors ) is
  begin 
    Clear(t.b);
  end Clear;

  procedure Clear ( t : in out DoblDobl_Numeric_Minors ) is
  begin 
    Clear(t.b);
  end Clear;

  procedure Clear ( t : in out QuadDobl_Numeric_Minors ) is
  begin 
    Clear(t.b);
  end Clear;

end Remember_Numeric_Minors;
