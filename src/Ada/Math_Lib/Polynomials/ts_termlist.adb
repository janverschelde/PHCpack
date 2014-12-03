with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Strings;
with Standard_Random_Polynomials;
with Standard_Complex_Term_Lists;
with Standard_Complex_Term_Lists_io;     use Standard_Complex_Term_Lists_io;

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

  procedure Main is

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
  end Main;

begin
  Main;
end ts_termlist;
