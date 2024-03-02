with Communications_with_User;         use Communications_with_User;
with File_Scanning;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;      use Standard_Complex_Numbers_io;
with Standard_Parse_Numbers;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;      use Multprec_Complex_Numbers_io;
with Multprec_Parse_Numbers;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;  use Standard_Complex_Polynomials_io;
with Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;  use Standard_Complex_Laurentials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io; use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io; use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Poly_Strings;    use Standard_Complex_Poly_Strings;
with Standard_Complex_Laur_Strings;    use Standard_Complex_Laur_Strings;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;  use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Laurentials_io;  use DoblDobl_Complex_Laurentials_io;
with DoblDobl_Complex_Poly_Strings;    use DoblDobl_Complex_Poly_Strings;
with DoblDobl_Complex_Laur_Strings;    use DoblDobl_Complex_Laur_Strings;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;  use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laurentials_io;  use QuadDobl_Complex_Laurentials_io;
with QuadDobl_Complex_Poly_Strings;    use QuadDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Laur_Strings;    use QuadDobl_Complex_Laur_Strings;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;  use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Laurentials;
with Multprec_Complex_Laurentials_io;  use Multprec_Complex_Laurentials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io; use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Laur_Systems;
with Multprec_Complex_Laur_Systems_io; use Multprec_Complex_Laur_Systems_io;
with Multprec_Complex_Poly_Strings;    use Multprec_Complex_Poly_Strings;
with Multprec_Complex_Laur_Strings;    use Multprec_Complex_Laur_Strings;

package body Test_Parse_Polynomials is

  procedure Parse_Standard_Number is

    s : constant string := Read_String;
    c : Standard_Complex_Numbers.Complex_Number; 
    p : integer := s'first;

  begin
    Standard_Parse_Numbers.Parse(s & " ",p,c);
    put("The number : "); put(c); new_line;
  end Parse_Standard_Number;

  procedure Parse_Multiprecision_Number is

    size : natural32 := 0;

  begin
    put("Give the size of the precision : "); get(size);
    skip_line;
    declare
      s : constant string := Read_String;
      c : Multprec_Complex_Numbers.Complex_Number;
      p : integer := s'first;
    begin
      Multprec_Parse_Numbers.Parse(s & " ",size,p,c);
      put("The number : "); put(c); new_line;
    end;
  end Parse_Multiprecision_Number;

  procedure Parse_Standard_Polynomial is

    use Standard_Complex_Polynomials;

    n : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    p := Parse(n,s);
    put_line("The parsed polynomial : "); put(p); new_line;
    Symbol_Table.Init(n);
    put("Give a polynomial : "); get(p);
    put("-> Your polynomial : "); put(p);
    new_line;
    declare
      ps : constant string := Write(p);
      sp : Poly;
    begin
      put("The polynomial : "); put(ps); new_line;
      sp := Parse(n,ps);
      put("-> After parsing : "); put(sp); new_line;
    end;
  end Parse_Standard_Polynomial;

  procedure Parse_DoblDobl_Polynomial is

    use DoblDobl_Complex_Polynomials;

    n : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    p := Parse(n,s);
    put_line("The parsed polynomial : "); put(p); new_line;
    new_line;
    declare
      ps : constant string := Write(p);
      sp : Poly;
    begin
      put("The polynomial : "); put(ps); new_line;
      sp := Parse(n,ps);
      put("-> After parsing : "); put(sp); new_line;
    end;
  end Parse_DoblDobl_Polynomial;

  procedure Parse_QuadDobl_Polynomial is

    use QuadDobl_Complex_Polynomials;

    n : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    p := Parse(n,s);
    put_line("The parsed polynomial : "); put(p); new_line;
    new_line;
    declare
      ps : constant string := Write(p);
      sp : Poly;
    begin
      put("The polynomial : "); put(ps); new_line;
      sp := Parse(n,ps);
      put("-> After parsing : "); put(sp); new_line;
    end;
  end Parse_QuadDobl_Polynomial;

  procedure Parse_Multprec_Polynomial is

    use Multprec_Complex_Polynomials;

    n,size : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    put("Give the size of the numbers : "); get(size);
    Symbol_Table.Init(n);
    p := Parse(n,size,s);
    put_line("The parsed polynomial : "); put(p); new_line;
    new_line;
    declare
      ps : constant string := Write(p);
      sp : Poly;
    begin
      put("The polynomial : "); put(ps); new_line;
      sp := Parse(n,size,ps);
      put("-> After parsing : "); put(sp); new_line;
    end;
  end Parse_Multprec_Polynomial;

  procedure Parse_Standard_Laurent_Polynomial is

    use Standard_Complex_Laurentials;

    n : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    p := Parse(n,s);
    put_line("The parsed polynomial : "); put(p); new_line;
    Symbol_Table.Init(n);
    put("Give a polynomial : "); get(p);
    put("-> Your polynomial : "); put(p);
    new_line;
    declare
      ps : constant string := Write(p);
      sp : Poly;
    begin
      put("The polynomial : "); put(ps); new_line;
      sp := Parse(n,ps);
      put("-> After parsing : "); put(sp); new_line;
    end;
  end Parse_Standard_Laurent_Polynomial;

  procedure Parse_DoblDobl_Laurent_Polynomial is

    use DoblDobl_Complex_Laurentials;

    n : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    p := Parse(n,s);
    put_line("The parsed polynomial : "); put(p); new_line;
    new_line;
    declare
      ps : constant string := Write(p);
      sp : Poly;
    begin
      put("The polynomial : "); put(ps); new_line;
      sp := Parse(n,ps);
      put("-> After parsing : "); put(sp); new_line;
    end;
  end Parse_DoblDobl_Laurent_Polynomial;

  procedure Parse_QuadDobl_Laurent_Polynomial is

    use QuadDobl_Complex_Laurentials;

    n : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    p := Parse(n,s);
    put_line("The parsed polynomial : "); put(p); new_line;
    new_line;
    declare
      ps : constant string := Write(p);
      sp : Poly;
    begin
      put("The polynomial : "); put(ps); new_line;
      sp := Parse(n,ps);
      put("-> After parsing : "); put(sp); new_line;
    end;
  end Parse_QuadDobl_Laurent_Polynomial;

  procedure Parse_Multprec_Laurent_Polynomial is

    use Multprec_Complex_Laurentials;

    n,size : natural32 := 0;
    s : constant string := Read_String;
    p : Poly;

  begin
    put("-> Your string : "); put_line(s);
    put("Give the number of variables : "); get(n);
    put("Give the size of the numbers : "); get(size);
    Symbol_Table.Init(n);
    p := Parse(n,size,s);
    put_line("The parsed polynomial : "); put(p); new_line;
  end Parse_Multprec_Laurent_Polynomial;

  procedure Parse_String_from_System is

    use Standard_Complex_Polynomials,Standard_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;

  begin
    get(lp);
    declare
      sp : constant string := Write(lp.all);
      n : constant natural32 := natural32(lp'last);
      m : constant natural32 := Number_of_Unknowns(lp(lp'last));
      p : Poly_Sys(1..integer32(n));
    begin
      put_line("The polynomial system as string : ");
      put(sp); new_line;
      p := Parse(n,m,sp);
      put_line("The system after parsing : "); put(p);
    end;
  end Parse_String_from_System;

  function Scan_for_Delimiter 
             ( file : in file_type; d : in character ) return string is

    max : constant natural := 1024;
    cnt : integer := 0;
    res : string(1..max);
    c : character;

  begin
    loop
      get(file,c);
      cnt := cnt + 1;
      exit when (cnt > max);
      res(cnt) := c;
      exit when (c = d);
    end loop;
    return res(1..cnt);
  end Scan_for_Delimiter;

  procedure Parse_String_from_File is

    file : file_type;
    n,m : natural32 := 0;

  begin
    put_line("Reading a file name for a system...");
    Read_Name_and_Open_File(file);
    get(file,n);
    new_line;
    put("The number of polynomials : "); put(n,1); new_line;
    m := natural32(File_Scanning.Scan_Line_for_Number(file));
    if Symbol_Table.Empty then
      if m > 0
       then Symbol_Table.Init(m);
       else Symbol_Table.Init(n);
      end if;
    end if;
    if m = 0
     then m := n;
    end if;
    put("The number of variables : "); put(m,1); new_line;
    for i in 1..n loop
      declare
        s : constant string := Scan_for_Delimiter(file,';');
      begin
        put("polynomial "); put(i,1); put(" : "); put(s); new_line;
      end;
    end loop;
  end Parse_String_from_File;

  procedure Read_Strings
              ( n,m : out natural32; p : out Link_to_Array_of_Strings ) is

    file : file_type;

  begin
    new_line;
    put_line("Reading a file name for a polynomial system...");
    Read_Name_and_Open_File(file);
    get(file,natural(n),natural(m),p);
    new_line;
    put("number of equations : "); put(n,1); new_line;
    put("number of variables : "); put(m,1); new_line;
    for i in 1..n loop
      put("p["); put(i,1); put("] : ");
      put(p(integer(i)).all); new_line;
    end loop;
  end Read_Strings;

  procedure Parse_Standard_Complex_Strings_from_File is

    n,m : natural32 := 0;
    p : Link_to_Array_of_Strings;

  begin
    Read_Strings(n,m,p);
    Symbol_Table.Init(m);
    declare
      q : constant Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n))
        := Parse(m,p.all);
    begin
      put_line("The parsed polynomial system : "); put(q);
    end;
  end Parse_Standard_Complex_Strings_from_File;

  procedure Parse_Standard_Laurent_Strings_from_File is

    n,m : natural32 := 0;
    p : Link_to_Array_of_Strings;

  begin
    Read_Strings(n,m,p);
    Symbol_Table.Init(m);
    declare
      q : constant Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(n))
        := Parse(m,p.all);
    begin
      put_line("The parsed Laurent polynomial system : "); put(q);
    end;
  end Parse_Standard_Laurent_Strings_from_File;

  procedure Parse_Multprec_Complex_Strings_from_File is

    n,m,size : natural32 := 0;
    p : Link_to_Array_of_Strings;

  begin
    Read_Strings(n,m,p);
    Symbol_Table.Init(m);
    put("Give the size of the numbers : "); get(size);
    declare
      q : constant Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n))
        := Parse(m,size,p.all);
    begin
      put_line("The parsed polynomial system : "); put(q);
    end;
  end Parse_Multprec_Complex_Strings_from_File;

  procedure Parse_Multprec_Laurent_Strings_from_File is

    n,m,size : natural32 := 0;
    p : Link_to_Array_of_Strings;

  begin
    Read_Strings(n,m,p);
    Symbol_Table.Init(m);
    put("Give the size of the numbers : "); get(size);
    declare
      q : constant Multprec_Complex_Laur_Systems.Laur_Sys(1..integer32(n))
        := Parse(m,size,p.all);
    begin
      put_line("The parsed Laurent polynomial system : "); put(q);
    end;
  end Parse_Multprec_Laurent_Strings_from_File;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Conversion of strings into polynomials, and vice versa.");
    new_line;
    put_line("MENU to test parsing : ");
    put_line("  1. parse a string for a standard hardware number;");
    put_line("  2. parse a string for a multiprecision number;");
    put_line("  3. test parsing of one standard common polynomial;");
    put_line("  4. test parsing of one double double common polynomial;");
    put_line("  5. test parsing of one quad double common polynomial;");
    put_line("  6. test parsing of one multiprecision common polynomial;");
    put_line("  7. test parsing of one Laurent polynomial;");
    put_line("  8. test parsing of one double double Laurent polynomial;");
    put_line("  9. test parsing of one quad double Laurent polynomial;");
    put_line("  A. test parsing of one multiprecision Laurent polynomial;");
    put_line("  B. parse a polynomial system, given from file;");
    put_line("  C. parse a string from file into a polynomial system;");
    put_line("  D. read strings from file, parse standard complex system;");
    put_line("  E. read strings from file, parse standard Laurent system;");
    put_line("  F. read strings from file, parse multprec complex system;");
    put_line("  G. read strings from file, parse multprec Laurent system.");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, or G : ");
    Ask_Alternative(ans,"123456789ABCDEFG");
    new_line;
    case ans is
      when '1' => Parse_Standard_Number;
      when '2' => Parse_Multiprecision_Number;
      when '3' => Parse_Standard_Polynomial;
      when '4' => Parse_DoblDobl_Polynomial;
      when '5' => Parse_QuadDobl_Polynomial;
      when '6' => Parse_Multprec_Polynomial;
      when '7' => Parse_Standard_Laurent_Polynomial;
      when '8' => Parse_DoblDobl_Laurent_Polynomial;
      when '9' => Parse_QuadDobl_Laurent_Polynomial;
      when 'A' => Parse_Multprec_Laurent_Polynomial;
      when 'B' => Parse_String_from_System;
      when 'C' => Parse_String_from_File;
      when 'D' => Parse_Standard_Complex_Strings_from_File;
      when 'E' => Parse_Standard_Laurent_Strings_from_File;
      when 'F' => Parse_Multprec_Complex_Strings_from_File;
      when 'G' => Parse_Multprec_Laurent_Strings_from_File;
      when others => null;
    end case;
  end Main;

end Test_Parse_Polynomials;
