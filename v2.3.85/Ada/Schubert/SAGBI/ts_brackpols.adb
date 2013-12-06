with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;         use Standard_Complex_Numbers;
with Brackets,Brackets_io;             use Brackets,Brackets_io;
with Bracket_Monomials;                use Bracket_Monomials;
with Bracket_Monomials_io;             use Bracket_Monomials_io;
with Bracket_Polynomials;              use Bracket_Polynomials;
with Bracket_Polynomials_io;           use Bracket_Polynomials_io;

procedure ts_brackpols is

  procedure Compare ( p : in Bracket_Polynomial; t : in Bracket_Term ) is

  -- DESCRIPTION :
  --   Compares every monomial in p with the term t.

    procedure Compare ( tt : in Bracket_Term; continue : out boolean ) is
    begin
      if tt < t
       then put(tt); put(" < "); put(t);
       elsif tt > t
           then put(tt); put(" > "); put(t);
           else put(tt); put(" = "); put(t);
      end if;
      new_line;
      if tt > t
       then put(tt); put(" > "); put(t);
       elsif tt < t
           then put(tt); put(" < "); put(t);
           else put(tt); put(" = "); put(t);
      end if;
      new_line;
      continue := true;
    end Compare;
    procedure Compare_Terms is new Enumerate_Terms(Compare);

  begin
    Compare_Terms(p);
  end Compare;

  procedure Interactive_Read ( d : in integer32; t : out Bracket_Term ) is

    m : natural32 := 0;
    b : Bracket(1..d);
    s,c : integer32;
    bm : Bracket_Monomial;

  begin
    c := 1;
    put("Give the number of brackets in the monomial : "); get(m);
    for i in 1..m loop
      put("Give "); put(d,1); put(" numbers for the ");
      put(i,1); put("th bracket : "); get(b,s);
      c := s*c;
      Multiply(bm,b);
    end loop;
    put_line("The bracket monomial : "); put(bm); new_line;
    if c > 0
     then t.coeff := Create(1.0);
     else t.coeff := -Create(1.0);
    end if;
    t.monom := bm;
  end Interactive_Read;

  procedure Interactive_Read
              ( d : in integer32; p : out Bracket_Polynomial ) is

    m : natural32 := 0;
    bp : Bracket_Polynomial;

  begin
    put_line("Reading a bracket polynomial.");
    put("Give the number of monomials : "); get(m);
    for i in 1..m loop
      declare
        bt : Bracket_Term;
      begin
        put("Reading bracket monomial "); put(i,1); new_line;
        Interactive_Read(d,bt);
        put_line("The bracket term : "); put(bt);
        Compare(bp,bt);
        Add(bp,bt);
        put_line("The bracket polynomial : "); put(bp);
      end;
    end loop;
    p := bp;
  end Interactive_Read;

  procedure Main is

    d : integer32 := 0;
    bp : Bracket_Polynomial;

  begin
    put("Give the number of entries in the brackets : "); get(d);
    Interactive_Read(d,bp);
    put_line("The bracket polynomial : "); put(bp);
  end Main;

begin
  Main;
end ts_brackpols;
