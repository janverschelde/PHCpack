with text_io;                          use text_io;
with Communications_with_User;         use Communications_with_User;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Brackets,Brackets_io;             use Brackets,Brackets_io;
with Bracket_Monomials;                use Bracket_Monomials;
with Bracket_Monomials_io;             use Bracket_Monomials_io;
with Standard_Bracket_Polynomials;
with Standard_Bracket_Polynomials_io;  use Standard_Bracket_Polynomials_io;
with DoblDobl_Bracket_Polynomials;
with DoblDobl_Bracket_Polynomials_io;  use DoblDobl_Bracket_Polynomials_io;
with QuadDobl_Bracket_Polynomials;
with QuadDobl_Bracket_Polynomials_io;  use QuadDobl_Bracket_Polynomials_io;

procedure ts_brackpols is

-- DESCRIPTION :
--   A simple test on the data structures to represent bracket polynomials
--   with coefficients in double, double double, or quad double precision.

  procedure Compare
              ( p : in Standard_Bracket_Polynomials.Bracket_Polynomial;
                t : in Standard_Bracket_Polynomials.Bracket_Term ) is

  -- DESCRIPTION :
  --   Compares every monomial in p with the term t.

    use Standard_Bracket_Polynomials;

    procedure Compare ( tt : in Bracket_Term; continue : out boolean ) is
    begin
      if tt < t then
        put(tt); put(" < "); put(t);
      elsif tt > t then
        put(tt); put(" > "); put(t);
      else
         put(tt); put(" = "); put(t);
      end if;
      new_line;
      if tt > t then
        put(tt); put(" > "); put(t);
      elsif tt < t then
        put(tt); put(" < "); put(t);
      else
        put(tt); put(" = "); put(t);
      end if;
      new_line;
      continue := true;
    end Compare;
    procedure Compare_Terms is new Enumerate_Terms(Compare);

  begin
    Compare_Terms(p);
  end Compare;

  procedure Compare
              ( p : in DoblDobl_Bracket_Polynomials.Bracket_Polynomial;
                t : in DoblDobl_Bracket_Polynomials.Bracket_Term ) is

  -- DESCRIPTION :
  --   Compares every monomial in p with the term t.

    use DoblDobl_Bracket_Polynomials;

    procedure Compare ( tt : in Bracket_Term; continue : out boolean ) is
    begin
      if tt < t then
        put(tt); put(" < "); put(t);
      elsif tt > t then
        put(tt); put(" > "); put(t);
      else
         put(tt); put(" = "); put(t);
      end if;
      new_line;
      if tt > t then
        put(tt); put(" > "); put(t);
      elsif tt < t then
        put(tt); put(" < "); put(t);
      else
        put(tt); put(" = "); put(t);
      end if;
      new_line;
      continue := true;
    end Compare;
    procedure Compare_Terms is new Enumerate_Terms(Compare);

  begin
    Compare_Terms(p);
  end Compare;

  procedure Compare
              ( p : in QuadDobl_Bracket_Polynomials.Bracket_Polynomial;
                t : in QuadDobl_Bracket_Polynomials.Bracket_Term ) is

  -- DESCRIPTION :
  --   Compares every monomial in p with the term t.

    use QuadDobl_Bracket_Polynomials;

    procedure Compare ( tt : in Bracket_Term; continue : out boolean ) is
    begin
      if tt < t then
        put(tt); put(" < "); put(t);
      elsif tt > t then
        put(tt); put(" > "); put(t);
      else
         put(tt); put(" = "); put(t);
      end if;
      new_line;
      if tt > t then
        put(tt); put(" > "); put(t);
      elsif tt < t then
        put(tt); put(" < "); put(t);
      else
        put(tt); put(" = "); put(t);
      end if;
      new_line;
      continue := true;
    end Compare;
    procedure Compare_Terms is new Enumerate_Terms(Compare);

  begin
    Compare_Terms(p);
  end Compare;

  procedure Interactive_Read
              ( d : in integer32;
                t : out Standard_Bracket_Polynomials.Bracket_Term ) is

  -- DESCRIPTION :
  --   Prompts the user for bracket monomials.

    use Standard_Complex_Numbers;
    use Standard_Bracket_Polynomials;

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
              ( d : in integer32;
                t : out DoblDobl_Bracket_Polynomials.Bracket_Term ) is

  -- DESCRIPTION :
  --   Prompts the user for bracket monomials.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Bracket_Polynomials;

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
     then t.coeff := Create(integer(1));
     else t.coeff := -Create(integer(1));
    end if;
    t.monom := bm;
  end Interactive_Read;

  procedure Interactive_Read
              ( d : in integer32;
                t : out QuadDobl_Bracket_Polynomials.Bracket_Term ) is

  -- DESCRIPTION :
  --   Prompts the user for bracket monomials.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Bracket_Polynomials;

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
     then t.coeff := Create(integer(1));
     else t.coeff := -Create(integer(1));
    end if;
    t.monom := bm;
  end Interactive_Read;

  procedure Interactive_Read
              ( d : in integer32;
                p : out Standard_Bracket_Polynomials.Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Reads a bracket polynomial, one monomial at a time.

    use Standard_Bracket_Polynomials;

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

  procedure Interactive_Read
              ( d : in integer32;
                p : out DoblDobl_Bracket_Polynomials.Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Reads a bracket polynomial, one monomial at a time.

    use DoblDobl_Bracket_Polynomials;

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

  procedure Interactive_Read
              ( d : in integer32;
                p : out QuadDobl_Bracket_Polynomials.Bracket_Polynomial ) is

  -- DESCRIPTION :
  --   Reads a bracket polynomial, one monomial at a time.

    use QuadDobl_Bracket_Polynomials;

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

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Test on bracket polynomials with coefficients
  --   in standard double precision.

    use Standard_Bracket_Polynomials;

    d : integer32 := 0;
    bp : Bracket_Polynomial;

  begin
    put("Give the number of entries in the brackets : "); get(d);
    Interactive_Read(d,bp);
    put_line("The bracket polynomial : "); put(bp);
  end Standard_Main;

  procedure DoblDobl_Main is

  -- DESCRIPTION :
  --   Test on bracket polynomials with coefficients
  --   in double double precision.

    use DoblDobl_Bracket_Polynomials;

    d : integer32 := 0;
    bp : Bracket_Polynomial;

  begin
    put("Give the number of entries in the brackets : "); get(d);
    Interactive_Read(d,bp);
    put_line("The bracket polynomial : "); put(bp);
  end DoblDobl_Main;

  procedure QuadDobl_Main is

  -- DESCRIPTION :
  --   Test on bracket polynomials with coefficients
  --   in quad double precision.

    use QuadDobl_Bracket_Polynomials;

    d : integer32 := 0;
    bp : Bracket_Polynomial;

  begin
    put("Give the number of entries in the brackets : "); get(d);
    Interactive_Read(d,bp);
    put_line("The bracket polynomial : "); put(bp);
  end QuadDobl_Main;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision and invokes
  --   the the appropriate test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Main;
      when '1' => DoblDobl_Main;
      when '2' => QuadDobl_Main;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_brackpols;
