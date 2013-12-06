with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Brackets_io;                       use Brackets_io;
with Bracket_Monomials;                 use Bracket_Monomials;
with Symbolic_Minor_Equations;          use Symbolic_Minor_Equations;
with Evaluated_Minors;                  use Evaluated_Minors;

package body Remember_Numeric_Minors is

  function Number_of_Minors ( n,k : natural32 ) return natural32 is
  begin
    return Symbolic_Minor_Equations.Number_of_Maximal_Minors(n,k);
  end Number_of_Minors;

  function Create ( n,k : natural32; x : Matrix ) return Numeric_Minor_Table is

    nq : constant integer32 := integer32(Number_of_Maximal_Minors(n,k));
    bm : Bracket_Monomial := Maximal_Minors(n,k);
    res : Numeric_Minor_Table(nq);
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

  function Search ( t : Numeric_Minor_Table;
                    b : Bracket ) return Complex_Number is
  begin
    for i in 1..t.m loop
      if Is_Equal(b,t.b(i).all)
       then return t.v(i);
      end if;
    end loop;
    return Create(0.0);
  end Search;

  procedure Write ( t : in Numeric_Minor_Table ) is
  begin
    for i in 1..t.m loop
      put(t.b(i).all); put(" : "); put(t.v(i)); new_line;
    end loop;
  end Write;

  procedure Query ( t : in Numeric_Minor_Table; k : in integer32 ) is

    b : Bracket(1..k);
    v : Complex_Number;
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

  procedure Clear ( t : in out Numeric_Minor_Table ) is
  begin 
    Clear(t.b);
  end Clear;

end Remember_Numeric_Minors;
