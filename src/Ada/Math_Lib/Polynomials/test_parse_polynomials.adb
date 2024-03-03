with Communications_with_User;         use Communications_with_User;
with File_Scanning;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;      use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;      use DoblDobl_Complex_Numbers_io;
with TripDobl_Complex_Numbers;
with TripDobl_Complex_Numbers_io;      use TripDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;      use QuadDobl_Complex_Numbers_io;
with PentDobl_Complex_Numbers;
with PentDobl_Complex_Numbers_io;      use PentDobl_Complex_Numbers_io;
with OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Numbers_io;      use OctoDobl_Complex_Numbers_io;
with DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Numbers_io;      use DecaDobl_Complex_Numbers_io;
with HexaDobl_Complex_Numbers;
with HexaDobl_Complex_Numbers_io;      use HexaDobl_Complex_Numbers_io;
with Standard_Parse_Numbers;
with DoblDobl_Complex_Numbers_cv;
with TripDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers_cv;
with PentDobl_Complex_Numbers_cv;
with OctoDobl_Complex_Numbers_cv;
with DecaDobl_Complex_Numbers_cv;
with HexaDobl_Complex_Numbers_cv;
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
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Polynomials_io;  use TripDobl_Complex_Polynomials_io;
with TripDobl_Complex_Laurentials;
with TripDobl_Complex_Laurentials_io;  use TripDobl_Complex_Laurentials_io;
with TripDobl_Complex_Poly_Strings;    use TripDobl_Complex_Poly_Strings;
with TripDobl_Complex_Laur_Strings;    use TripDobl_Complex_Laur_Strings;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;  use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Laurentials_io;  use QuadDobl_Complex_Laurentials_io;
with QuadDobl_Complex_Poly_Strings;    use QuadDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Laur_Strings;    use QuadDobl_Complex_Laur_Strings;
with PentDobl_Complex_Polynomials;
with PentDobl_Complex_Polynomials_io;  use PentDobl_Complex_Polynomials_io;
with PentDobl_Complex_Laurentials;
with PentDobl_Complex_Laurentials_io;  use PentDobl_Complex_Laurentials_io;
with PentDobl_Complex_Poly_Strings;    use PentDobl_Complex_Poly_Strings;
with PentDobl_Complex_Laur_Strings;    use PentDobl_Complex_Laur_Strings;
with OctoDobl_Complex_Polynomials;
with OctoDobl_Complex_Polynomials_io;  use OctoDobl_Complex_Polynomials_io;
with OctoDobl_Complex_Laurentials;
with OctoDobl_Complex_Laurentials_io;  use OctoDobl_Complex_Laurentials_io;
with OctoDobl_Complex_Poly_Strings;    use OctoDobl_Complex_Poly_Strings;
with OctoDobl_Complex_Laur_Strings;    use OctoDobl_Complex_Laur_Strings;
with DecaDobl_Complex_Polynomials;
with DecaDobl_Complex_Polynomials_io;  use DecaDobl_Complex_Polynomials_io;
with DecaDobl_Complex_Laurentials;
with DecaDobl_Complex_Laurentials_io;  use DecaDobl_Complex_Laurentials_io;
with DecaDobl_Complex_Poly_Strings;    use DecaDobl_Complex_Poly_Strings;
with DecaDobl_Complex_Laur_Strings;    use DecaDobl_Complex_Laur_Strings;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Polynomials_io;  use HexaDobl_Complex_Polynomials_io;
with HexaDobl_Complex_Laurentials;
with HexaDobl_Complex_Laurentials_io;  use HexaDobl_Complex_Laurentials_io;
with HexaDobl_Complex_Poly_Strings;    use HexaDobl_Complex_Poly_Strings;
with HexaDobl_Complex_Laur_Strings;    use HexaDobl_Complex_Laur_Strings;
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

  procedure Parse_DoblDobl_Number is

    s : constant string := Read_String;
    ddc : DoblDobl_Complex_Numbers.Complex_Number; 
    mpc : Multprec_Complex_Numbers.Complex_Number; 
    p : integer := s'first;
    size : constant natural32 := 5;

  begin
    Multprec_Parse_Numbers.Parse(s & " ",size,p,mpc);
    ddc := DoblDobl_Complex_Numbers_cv.Multprec_to_DoblDobl_Complex(mpc);
    put("The number : "); put(ddc); new_line;
  end Parse_DoblDobl_Number;

  procedure Parse_TripDobl_Number is

    s : constant string := Read_String;
    tdc : TripDobl_Complex_Numbers.Complex_Number; 
    mpc : Multprec_Complex_Numbers.Complex_Number; 
    p : integer := s'first;
    size : constant natural32 := 8;

  begin
    Multprec_Parse_Numbers.Parse(s & " ",size,p,mpc);
    tdc := TripDobl_Complex_Numbers_cv.Multprec_to_TripDobl_Complex(mpc);
    put("The number : "); put(tdc); new_line;
  end Parse_TripDobl_Number;

  procedure Parse_QuadDobl_Number is

    s : constant string := Read_String;
    qdc : QuadDobl_Complex_Numbers.Complex_Number; 
    mpc : Multprec_Complex_Numbers.Complex_Number; 
    p : integer := s'first;
    size : constant natural32 := 8;

  begin
    Multprec_Parse_Numbers.Parse(s & " ",size,p,mpc);
    qdc := QuadDobl_Complex_Numbers_cv.Multprec_to_QuadDobl_Complex(mpc);
    put("The number : "); put(qdc); new_line;
  end Parse_QuadDobl_Number;

  procedure Parse_PentDobl_Number is

    s : constant string := Read_String;
    pdc : PentDobl_Complex_Numbers.Complex_Number; 
    mpc : Multprec_Complex_Numbers.Complex_Number; 
    p : integer := s'first;
    size : constant natural32 := 10;

  begin
    Multprec_Parse_Numbers.Parse(s & " ",size,p,mpc);
    pdc := PentDobl_Complex_Numbers_cv.Multprec_to_PentDobl_Complex(mpc);
    put("The number : "); put(pdc); new_line;
  end Parse_PentDobl_Number;

  procedure Parse_OctoDobl_Number is

    s : constant string := Read_String;
    odc : OctoDobl_Complex_Numbers.Complex_Number; 
    mpc : Multprec_Complex_Numbers.Complex_Number; 
    p : integer := s'first;
    size : constant natural32 := 16;

  begin
    Multprec_Parse_Numbers.Parse(s & " ",size,p,mpc);
    odc := OctoDobl_Complex_Numbers_cv.Multprec_to_OctoDobl_Complex(mpc);
    put("The number : "); put(odc); new_line;
  end Parse_OctoDobl_Number;

  procedure Parse_DecaDobl_Number is

    s : constant string := Read_String;
    dac : DecaDobl_Complex_Numbers.Complex_Number; 
    mpc : Multprec_Complex_Numbers.Complex_Number; 
    p : integer := s'first;
    size : constant natural32 := 20;

  begin
    Multprec_Parse_Numbers.Parse(s & " ",size,p,mpc);
    dac := DecaDobl_Complex_Numbers_cv.Multprec_to_DecaDobl_Complex(mpc);
    put("The number : "); put(dac); new_line;
  end Parse_DecaDobl_Number;

  procedure Parse_HexaDobl_Number is

    s : constant string := Read_String;
    hdc : HexaDobl_Complex_Numbers.Complex_Number; 
    mpc : Multprec_Complex_Numbers.Complex_Number; 
    p : integer := s'first;
    size : constant natural32 := 32;

  begin
    Multprec_Parse_Numbers.Parse(s & " ",size,p,mpc);
    hdc := HexaDobl_Complex_Numbers_cv.Multprec_to_HexaDobl_Complex(mpc);
    put("The number : "); put(hdc); new_line;
  end Parse_HexaDobl_Number;

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

  procedure Parse_TripDobl_Polynomial is

    use TripDobl_Complex_Polynomials;

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
  end Parse_TripDobl_Polynomial;

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

  procedure Parse_PentDobl_Polynomial is

    use PentDobl_Complex_Polynomials;

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
  end Parse_PentDobl_Polynomial;

  procedure Parse_OctoDobl_Polynomial is

    use OctoDobl_Complex_Polynomials;

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
  end Parse_OctoDobl_Polynomial;

  procedure Parse_DecaDobl_Polynomial is

    use DecaDobl_Complex_Polynomials;

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
  end Parse_DecaDobl_Polynomial;

  procedure Parse_HexaDobl_Polynomial is

    use HexaDobl_Complex_Polynomials;

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
  end Parse_HexaDobl_Polynomial;

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

  procedure Parse_TripDobl_Laurent_Polynomial is

    use TripDobl_Complex_Laurentials;

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
  end Parse_TripDobl_Laurent_Polynomial;

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

  procedure Parse_PentDobl_Laurent_Polynomial is

    use PentDobl_Complex_Laurentials;

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
  end Parse_PentDobl_Laurent_Polynomial;

  procedure Parse_OctoDobl_Laurent_Polynomial is

    use OctoDobl_Complex_Laurentials;

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
  end Parse_OctoDobl_Laurent_Polynomial;

  procedure Parse_DecaDobl_Laurent_Polynomial is

    use DecaDobl_Complex_Laurentials;

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
  end Parse_DecaDobl_Laurent_Polynomial;

  procedure Parse_HexaDobl_Laurent_Polynomial is

    use HexaDobl_Complex_Laurentials;

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
  end Parse_HexaDobl_Laurent_Polynomial;

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

  function Prompt_for_Precision return character is

    res : character;

  begin
    new_line;
    put_line("MENU for the precision : ");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. triple double precision");
    put_line("  3. quad double precision");
    put_line("  4. penta double precision");
    put_line("  5. octo double precision");
    put_line("  6. deca double precision");
    put_line("  7. hexa double precision");
    put_line("  8. arbitrary multiprecision");
    put("Type 0, 1, 2, 3, 4, 5, 6, 7, or 8 to select the precision : ");
    Ask_Alternative(res,"012345678");
    return res;
  end Prompt_for_Precision;

  procedure Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in double precision : ");
    put_line("  1. parse a string for a double number");
    put_line("  2. parse for one polynomial");
    put_line("  3. parse for one Laurent polynomial");
    put_line("  4. parse a system, given from file");
    put_line("  5. parse a string from file into a system");
    put_line("  6. read strings from file, parse system");
    put_line("  7. read strings from file, parse Laurent system");
    put("Type 1, 2, 3, 4, 5, 6, or 7 : ");
    Ask_Alternative(ans,"1234567");
    new_line;
    case ans is 
      when '1' => Parse_Standard_Number;
      when '2' => Parse_Standard_Polynomial;
      when '3' => Parse_Standard_Laurent_Polynomial;
      when '4' => Parse_String_from_System;
      when '5' => Parse_String_from_File;
      when '6' => Parse_Standard_Complex_Strings_from_File;
      when '7' => Parse_Standard_Laurent_Strings_from_File;
      when others => null;
    end case;
  end Double_Test;

  procedure Double_Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in double double precision : ");
    put_line("  1. parse a string for a double double number");
    put_line("  2. parse for one double double polynomial");
    put_line("  3. parse for one double double Laurent polynomial");
    put("Type 1, 2, or 3 : "); Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Parse_DoblDobl_Number;
      when '2' => Parse_DoblDobl_Polynomial;
      when '3' => Parse_DoblDobl_Laurent_Polynomial;
      when others => null;
    end case;
  end Double_Double_Test;

  procedure Triple_Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in triple double precision : ");
    put_line("  1. parse a string for a triple double number");
    put_line("  2. parse for one triple double polynomial");
    put_line("  3. parse for one triple double Laurent polynomial");
    put("Type 1, 2, or 3 : "); Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Parse_TripDobl_Number;
      when '2' => Parse_TripDobl_Polynomial;
      when '3' => Parse_TripDobl_Laurent_Polynomial;
      when others => null;
    end case;
  end Triple_Double_Test;

  procedure Quad_Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in quad double precision : ");
    put_line("  1. parse a string for a quad double number");
    put_line("  1. parse for one quad double polynomial");
    put_line("  2. parse for one quad double Laurent polynomial");
    put("Type 1, 2, or 3 : "); Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Parse_QuadDobl_Number;
      when '2' => Parse_QuadDobl_Polynomial;
      when '3' => Parse_QuadDobl_Laurent_Polynomial;
      when others => null;
    end case;
  end Quad_Double_Test;

  procedure Penta_Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in penta double precision : ");
    put_line("  1. parse a string for a penta double number");
    put_line("  2. parse for one penta double polynomial");
    put_line("  3. parse for one penta double Laurent polynomial");
    put("Type 1, 2, or 3 : "); Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Parse_PentDobl_Number;
      when '2' => Parse_PentDobl_Polynomial;
      when '3' => Parse_PentDobl_Laurent_Polynomial;
      when others => null;
    end case;
  end Penta_Double_Test;

  procedure Octo_Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in octo double precision : ");
    put_line("  1. parse a string for a octo double number");
    put_line("  2. parse for one octo double polynomial");
    put_line("  3. parse for one octo double Laurent polynomial");
    put("Type 1, 2, or 3 : "); Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Parse_OctoDobl_Number;
      when '2' => Parse_OctoDobl_Polynomial;
      when '3' => Parse_OctoDobl_Laurent_Polynomial;
      when others => null;
    end case;
  end Octo_Double_Test;

  procedure Deca_Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in deca double precision : ");
    put_line("  1. parse a string for a deca double number");
    put_line("  2. parse for one deca double polynomial");
    put_line("  3. parse for one deca double Laurent polynomial");
    put("Type 1, 2, or 3 : "); Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Parse_DecaDobl_Number;
      when '2' => Parse_DecaDobl_Polynomial;
      when '3' => Parse_DecaDobl_Laurent_Polynomial;
      when others => null;
    end case;
  end Deca_Double_Test;

  procedure Hexa_Double_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in hexa double precision : ");
    put_line("  1. parse a string for a hexa double number");
    put_line("  2. parse for one hexa double polynomial");
    put_line("  3. parse for one hexa double Laurent polynomial");
    put("Type 1, 2, or 3 : "); Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => Parse_HexaDobl_Number;
      when '2' => Parse_HexaDobl_Polynomial;
      when '3' => Parse_HexaDobl_Laurent_Polynomial;
      when others => null;
    end case;
  end Hexa_Double_Test;

  procedure Multiprecision_Test is

    ans : character;

  begin
    put_line("MENU to test parsing in multiprecision : ");
    put_line("  1. parse a string for a multiprecision number");
    put_line("  2. parse for one multiprecision polynomial");
    put_line("  3. parse for one multiprecision Laurent polynomial");
    put_line("  4. read strings from file, parse system");
    put_line("  5. read strings from file, parse Laurent system");
    put("Type 1, 2, 3, 4, or 5 : ");
    Ask_Alternative(ans,"12345");
    new_line;
    case ans is
      when '1' => Parse_Multiprecision_Number;
      when '2' => Parse_Multprec_Polynomial;
      when '3' => Parse_Multprec_Laurent_Polynomial;
      when '4' => Parse_Multprec_Complex_Strings_from_File;
      when '5' => Parse_Multprec_Laurent_Strings_from_File;
      when others => null;
    end case;
  end Multiprecision_Test;

  procedure Main is

    prc : character;

  begin
    new_line;
    put_line("Parsing strings for polynomials ...");
    prc := Prompt_for_Precision;
    new_line;
    case prc is
      when '0' => Double_Test;
      when '1' => Double_Double_Test;
      when '2' => Triple_Double_Test;
      when '3' => Quad_Double_Test;
      when '4' => Penta_Double_Test;
      when '5' => Octo_Double_Test;
      when '6' => Deca_Double_Test;
      when '7' => Hexa_Double_Test;
      when '8' => Multiprecision_Test;
      when others => null;
    end case;
  end Main;

end Test_Parse_Polynomials;
