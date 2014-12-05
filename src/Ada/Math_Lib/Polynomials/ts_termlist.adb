with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Strings;
with Standard_Random_Polynomials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Strings;
with DoblDobl_Random_Polynomials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Strings;
with QuadDobl_Random_Polynomials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Strings;
with Multprec_Random_Polynomials;
with Standard_Complex_Term_Lists;
with Standard_Complex_Term_Lists_io;     use Standard_Complex_Term_Lists_io;
with DoblDobl_Complex_Term_Lists;
with DoblDobl_Complex_Term_Lists_io;     use DoblDobl_Complex_Term_Lists_io;
with QuadDobl_Complex_Term_Lists;
with QuadDobl_Complex_Term_Lists_io;     use QuadDobl_Complex_Term_Lists_io;
with Multprec_Complex_Term_Lists;
with Multprec_Complex_Term_Lists_io;     use Multprec_Complex_Term_Lists_io;

procedure ts_termlist is

-- DESCRIPTION :
--   Test on the package list of terms.

  procedure Standard_Test_Parse ( n : in natural32; s : in string ) is

  -- DESCRIPTION :
  --   Tests the parsing of string into a list of terms,
  --   in n variable.

    use Standard_Complex_Term_Lists;
    t : Term_List := Standard_Complex_Poly_Strings.Parse(n,s);
    p : Standard_Complex_Polynomials.Poly;

  begin
   -- put_line("The string :"); put_line(s);
    put_line("The list of parsed terms :"); put(t);
    p := Create(t);
    Clear(t);
    put("The reconstructed polynomial :"); put_line(p);
  end Standard_Test_Parse;

  procedure Multprec_Test_Parse ( n : in natural32; s : in string ) is

  -- DESCRIPTION :
  --   Tests the parsing of string into a list of terms,
  --   in n variable.

    use Multprec_Complex_Term_Lists;
    size : constant natural32 := 5;
    t : Term_List := Multprec_Complex_Poly_Strings.Parse(n,size,s);
    p : Multprec_Complex_Polynomials.Poly;

  begin
   -- put_line("The string :"); put_line(s);
    put_line("The list of parsed terms :"); put(t);
    p := Create(t);
   -- Clear(t);
    put("The reconstructed polynomial :"); put_line(p);
  end Multprec_Test_Parse;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Tests the creation of a list of terms on a random polynomial
  --   with complex coefficients in standard double precision.

    n,d,m : natural32 := 0;
    p : Standard_Complex_Polynomials.Poly;
    t : Standard_Complex_Term_Lists.Term_List;
    s : Link_to_String;

  begin
    new_line;
    put_line("Generating a random polynomial with complex coefficients ...");
    put("Give the number of variables : "); get(n);
    put("Give the largest degree : "); get(d);
    put("Give the number of terms : "); get(m);
    p := Standard_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
    put("The random polynomial : "); put_line(p);
    t := Standard_Complex_Term_Lists.Create(p);
    put_line("The list of terms :"); put(t);
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(n)));
    s := new string'(Standard_Complex_Poly_Strings.Write(p));
    put_line("The string respresentation of the polynomial :");
    put_line(s.all);
    Standard_Complex_Term_Lists.Clear(t);
    put_line("The list of terms after clear :"); put(t);
    Standard_Test_Parse(n,s.all);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Tests the creation of a list of terms on a random polynomial
  --   with complex coefficients in double double precision.

    n,d,m : natural32 := 0;
    p : DoblDobl_Complex_Polynomials.Poly;
    t : DoblDobl_Complex_Term_Lists.Term_List;
    s : Link_to_String;

  begin
    new_line;
    put_line("Generating a random polynomial with complex coefficients ...");
    put("Give the number of variables : "); get(n);
    put("Give the largest degree : "); get(d);
    put("Give the number of terms : "); get(m);
    p := DoblDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
    put("The random polynomial : "); put_line(p);
    t := DoblDobl_Complex_Term_Lists.Create(p);
    put_line("The list of terms :"); put(t);
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(n)));
    s := new string'(DoblDobl_Complex_Poly_Strings.Write(p));
    put_line("The string respresentation of the polynomial :");
    put_line(s.all);
    DoblDobl_Complex_Term_Lists.Clear(t);
    put_line("The list of terms after clear :"); put(t);
   -- DoblDobl_Test_Parse(n,s.all);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Tests the creation of a list of terms on a random polynomial
  --   with complex coefficients in quad double precision.

    n,d,m : natural32 := 0;
    p : QuadDobl_Complex_Polynomials.Poly;
    t : QuadDobl_Complex_Term_Lists.Term_List;
    s : Link_to_String;

  begin
    new_line;
    put_line("Generating a random polynomial with complex coefficients ...");
    put("Give the number of variables : "); get(n);
    put("Give the largest degree : "); get(d);
    put("Give the number of terms : "); get(m);
    p := QuadDobl_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
    put("The random polynomial : "); put_line(p);
    t := QuadDobl_Complex_Term_Lists.Create(p);
    put_line("The list of terms :"); put(t);
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(n)));
    s := new string'(QuadDobl_Complex_Poly_Strings.Write(p));
    put_line("The string respresentation of the polynomial :");
    put_line(s.all);
    QuadDobl_Complex_Term_Lists.Clear(t);
    put_line("The list of terms after clear :"); put(t);
   -- QuadDobl_Test_Parse(n,s.all);
  end QuadDobl_Test;

  procedure Multprec_Test is

  -- DESCRIPTION :
  --   Tests the creation of a list of terms on a random polynomial
  --   with complex coefficients in arbitrary multiprecision.

    n,d,m : natural32 := 0;
    p : Multprec_Complex_Polynomials.Poly;
    t : Multprec_Complex_Term_Lists.Term_List;
    s : Link_to_String;

  begin
    new_line;
    put_line("Generating a random polynomial with complex coefficients ...");
    put("Give the number of variables : "); get(n);
    put("Give the largest degree : "); get(d);
    put("Give the number of terms : "); get(m);
    p := Multprec_Random_Polynomials.Random_Sparse_Poly(n,d,m,0);
    put("The random polynomial : "); put_line(p);
    t := Multprec_Complex_Term_Lists.Create(p);
    put_line("The list of terms :"); put(t);
    Symbol_Table.Init(Symbol_Table.Standard_Symbols(integer32(n)));
    s := new string'(Multprec_Complex_Poly_Strings.Write(p));
    put_line("The string respresentation of the polynomial :");
    put_line(s.all);
    Multprec_Complex_Term_Lists.Clear(t);
    put_line("The list of terms after clear :"); put(t);
    declare
      q : Multprec_Complex_Polynomials.Poly
        := Multprec_Complex_Poly_Strings.Parse(n,5,s.all);
    begin
      put("The random polynomial parsed from string :"); put_line(q);
    end;
    Multprec_Test_Parse(n,s.all);
  end Multprec_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to choose the precision of the tests.

    ans : character;

  begin
    new_line;
    put_line("MENU to test operations on lists of terms ...");
    put_line("  0. in standard double precision;");
    put_line("  1. in double double precision;");
    put_line("  2. in quad double precision;");
    put_line("  3. in arbitrary multiprecision;");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when '3' => Multprec_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_termlist;
