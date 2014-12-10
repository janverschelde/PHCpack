with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Multprec_Floating_Numbers;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_Strings;
with Standard_Random_Polynomials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Strings;
with DoblDobl_Random_Polynomials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Strings;
with QuadDobl_Random_Polynomials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
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

  procedure DoblDobl_Test_Parse ( n : in natural32; s : in string ) is

  -- DESCRIPTION :
  --   Tests the parsing of string into a list of terms,
  --   in n variable.

    use DoblDobl_Complex_Term_Lists;
    t : Term_List := DoblDobl_Complex_Poly_Strings.Parse(n,s);
    p : DoblDobl_Complex_Polynomials.Poly;

  begin
   -- put_line("The string :"); put_line(s);
    put_line("The list of parsed terms :"); put(t);
    p := Create(t);
    Clear(t);
    put("The reconstructed polynomial :"); put_line(p);
  end DoblDobl_Test_Parse;

  procedure QuadDobl_Test_Parse ( n : in natural32; s : in string ) is

  -- DESCRIPTION :
  --   Tests the parsing of string into a list of terms,
  --   in n variable.

    use QuadDobl_Complex_Term_Lists;
    t : Term_List := QuadDobl_Complex_Poly_Strings.Parse(n,s);
    p : QuadDobl_Complex_Polynomials.Poly;

  begin
   -- put_line("The string :"); put_line(s);
    put_line("The list of parsed terms :"); put(t);
    p := Create(t);
    Clear(t);
    put("The reconstructed polynomial :"); put_line(p);
  end QuadDobl_Test_Parse;

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
    Clear(t);
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
    DoblDobl_Test_Parse(n,s.all);
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
    QuadDobl_Test_Parse(n,s.all);
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

  procedure Standard_Test_User_Data
              ( nq,nv : in natural32; p : in Array_of_Strings ) is

  -- DESCRIPTION :
  --   Times the parsing of the array of strings in nv variables
  --   into lists of terms with standard double precision.
  --   The strings in p represent nq polynomials in nv variables.
  --   To make the comparison in efficiency, the user can also ask
  --   to parse the strings to polynomial data structures.

    use Standard_Complex_Term_Lists;
    use Standard_Complex_Poly_Systems;
    t : Array_of_Term_Lists(1..integer32(nq));
    f : Poly_Sys(t'range);
    timer : Timing_Widget;
    ans : character;
 
  begin
    Symbol_Table.Init(nv);
    new_line;
    put("Parse to a polynomial (if no, to list of terms) ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      f := Standard_Complex_Poly_Strings.Parse(nv,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the polynomials ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The polynomial system : "); put(f);
      end if;
    else
      tstart(timer);
      t := Standard_Complex_Poly_Strings.Parse(nv,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the lists of terms ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The lists of terms : "); put(t);
      end if;
    end if;
  end Standard_Test_User_Data;

  procedure DoblDobl_Test_User_Data
              ( nq,nv : in natural32; p : in Array_of_Strings ) is

  -- DESCRIPTION :
  --   Times the parsing of the array of strings in nv variables
  --   into lists of terms with double double precision.
  --   The strings in p represent nq polynomials in nv variables.

    use DoblDobl_Complex_Term_Lists;
    use DoblDobl_Complex_Poly_Systems;
    t : Array_of_Term_Lists(1..integer32(nq));
    f : Poly_Sys(t'range);
    timer : Timing_Widget;
    ans : character;
 
  begin
    Symbol_Table.Init(nv);
    new_line;
    put("Parse to a polynomial (if no, to list of terms) ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      f := DoblDobl_Complex_Poly_Strings.Parse(nv,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the polynomials ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The polynomial system : "); put(f);
      end if;
    else
      tstart(timer);
      t := DoblDobl_Complex_Poly_Strings.Parse(nv,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the lists of terms ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The lists of terms : "); put(t);
      end if;
    end if;
  end DoblDobl_Test_User_Data;

  procedure QuadDobl_Test_User_Data
              ( nq,nv : in natural32; p : in Array_of_Strings ) is

  -- DESCRIPTION :
  --   Times the parsing of the array of strings in nv variables
  --   into lists of terms with quad double precision.
  --   The strings in p represent nq polynomials in nv variables.

    use QuadDobl_Complex_Term_Lists;
    use QuadDobl_Complex_Poly_Systems;
    t : Array_of_Term_Lists(1..integer32(nq));
    f : Poly_Sys(t'range);
    timer : Timing_Widget;
    ans : character;
 
  begin
    Symbol_Table.Init(nv);
    new_line;
    put("Parse to a polynomial (if no, to list of terms) ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      f := QuadDobl_Complex_Poly_Strings.Parse(nv,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the polynomials ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The polynomial system : "); put(f);
      end if;
    else
      tstart(timer);
      t := QuadDobl_Complex_Poly_Strings.Parse(nv,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the lists of terms ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The lists of terms : "); put(t);
      end if;
    end if;
  end QuadDobl_Test_User_Data;

  procedure Multprec_Test_User_Data
              ( nq,nv : in natural32; p : in Array_of_Strings ) is

  -- DESCRIPTION :
  --   Times the parsing of the array of strings in nv variables
  --   into lists of terms with quad double precision.
  --   The strings in p represent nq polynomials in nv variables.

    use Multprec_Complex_Term_Lists;
    use Multprec_Complex_Poly_Systems;
    t : Array_of_Term_Lists(1..integer32(nq));
    f : Poly_Sys(t'range);
    timer : Timing_Widget;
    ans : character;
    deci,size : natural32 := 0;
 
  begin
    new_line;
    put("Give the number of decimal places in the precision : ");
    get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Symbol_Table.Init(nv);
    new_line;
    put("Parse to a polynomial (if no, to list of terms) ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      f := Multprec_Complex_Poly_Strings.Parse(nv,size,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the polynomials ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The polynomial system : "); put(f);
      end if;
    else
      tstart(timer);
      t := Multprec_Complex_Poly_Strings.Parse(nv,size,p);
      tstop(timer);
      new_line;
      print_times(standard_output,timer,"parsing the strings into term lists");
      new_line;
      put("Do you want to see the lists of terms ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        new_line;
        put_line("The lists of terms : "); put(t);
      end if;
    end if;
  end Multprec_Test_User_Data;

  procedure Test_User_Data is

  -- DESCRIPTION :
  --   Prompts the user to provide a file name and then reads
  --   an array of strings from that file for parsing.

    file : file_type;
    ls : Link_to_Array_of_Strings;
    timer : Timing_Widget;
    nq,nv : natural32;
    ans : character;

  begin
    new_line;
    put_line("Reading an array of strings from file ...");
    Read_Name_and_Open_File(file);
    tstart(timer);
    get(file,natural(nq),natural(nv),ls);
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"reading strings from file");
    new_line;
    put("Do you want to see the strings ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("The polynomials read from file :");
      for i in ls'range loop
        put_line(ls(i).all);
      end loop;
    end if;
    new_line;
    put_line("MENU for the precision level of parsing the coefficients :");
    put_line("  0. parse coefficients in standard double precision;");
    put_line("  1. parse coefficients in double double precision;");
    put_line("  2. parse coefficients in quad double precision;");
    put_line("  3. parse coefficients in arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Standard_Test_User_Data(nq,nv,ls.all);
      when '1' => DoblDobl_Test_User_Data(nq,nv,ls.all);
      when '2' => QuadDobl_Test_User_Data(nq,nv,ls.all);
      when '3' => Multprec_Test_User_Data(nq,nv,ls.all);
      when others => null;
    end case;
  end Test_User_Data;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to choose the precision of the tests.

    ans : character;

  begin
    new_line;
    put_line("MENU to test operations on lists of terms ...");
    put_line("  0. randomly generated data in standard double precision;");
    put_line("  1. randomly generated data in double double precision;");
    put_line("  2. randomly generated data in quad double precision;");
    put_line("  3. randomly generated data in arbitrary multiprecision;");
    put_line("  4. test on system given on file.");
    put("Type 0, 1, 2, 3, or 4 to select the test : ");
    Ask_Alternative(ans,"01234");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when '3' => Multprec_Test;
      when '4' => Test_USer_Data;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_termlist;
