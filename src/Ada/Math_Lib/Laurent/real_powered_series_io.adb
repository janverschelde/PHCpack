with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Strings_and_Numbers;
with Standard_Write_Numbers;
with Standard_Parse_Numbers;

package body Real_Powered_Series_IO is

-- WRITE OUTPUT :

  function append_term ( k : integer32; s : string;
                         c : Standard_Complex_Vectors.Vector;
                         p : Standard_Floating_Vectors.Vector;
                         t : character := 't';
                         vrblvl : integer32 := 0 ) return string is

  -- DESCRIPTION :
  --   Appends the k-th term to the string s, if k <= c'last.

  begin
    if vrblvl > 0 then
      put("-> in Real_Powered_Series_IO.append_term, k = ");
      put(k,1); put_line(" ...");
    end if;
    if k > c'last then
      return s;
    else
      declare
        cff : constant string := Strings_and_Numbers.Convert(c(k));
        pwr : constant string := Strings_and_Numbers.Convert(p(k));
        trm : constant string := s & " + " & cff & '*' & t & "**" & pwr;
      begin
        return append_term(k+1,trm,c,p,t);
      end;
    end if;
  end append_term;

  function to_string ( c : Standard_Complex_Vectors.Vector;
                       p : Standard_Floating_Vectors.Vector;
                       t : character := 't';
                       vrblvl : integer32 := 0 ) return string is

    cst : constant string := Strings_and_Numbers.Convert(c(0));

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Series_IO.to_string ...");
    end if;
    return append_term(1,cst,c,p,t,vrblvl-1);
  end to_string;

  procedure put ( c : in Standard_Complex_Vectors.Vector;
                  p : in Standard_Floating_Vectors.Vector;
                  t : in character := 't' ) is
  begin
    put(standard_output,c,p,t);
  end put;

  procedure put ( file : in file_type;
                  c : in Standard_Complex_Vectors.Vector;
                  p : in Standard_Floating_Vectors.Vector;
                  t : in character := 't' ) is

    cnt : natural32 := 0;

  begin
    Standard_Write_Numbers.Write_Number(file,c(0),cnt);
    for i in p'range loop
      put(file," + ");
      Standard_Write_Numbers.Write_Coefficient(file,c(i),cnt);
      put(file,t);
      put(file,"**");
      put(file,p(i),1,14,3);
    end loop;
  end put;

  procedure put_line ( c : in Standard_Complex_Vectors.Vector;
                       p : in Standard_Floating_Vectors.Vector;
                       t : in character := 't' ) is
  begin
    put_line(standard_output,c,p,t);
  end put_line;

  procedure put_line ( file : in file_type;
                       c : in Standard_Complex_Vectors.Vector;
                       p : in Standard_Floating_Vectors.Vector;
                       t : in character := 't' ) is

    cnt : natural32 := 0;
  
  begin
    Standard_Write_Numbers.Write_Number(file,c(0),cnt); new_line(file);
    for i in p'range loop
      put(file," + ");
      Standard_Write_Numbers.Write_Coefficient(file,c(i),cnt);
      put(file,t);
      put(file,"**");
      put(file,p(i),1,14,3); new_line(file);
    end loop;
  end put_line;

-- PARSE INPUT :

  function Size_Count
             ( s : string; t : character := 't' ) return integer32 is

    res : integer32 := 0;

  begin
    for i in s'range loop
      if s(i) = t
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Size_Count;

  procedure parse_string
              ( s : in string;
                c : out Standard_Complex_Vectors.Link_to_Vector;
                p : out Standard_Floating_Vectors.Link_to_Vector;
                t : in character := 't'; vrblvl : in integer32 := 0 ) is

    cst : Complex_Number;
    re_cst,im_cst : double_float;
    pos : integer := s'first;
    size : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Series_IO.parse_string ...");
    end if;
    while pos <= s'last loop -- find the opening round bracket
      exit when s(pos) = '(';
      pos := pos + 1;
    end loop;
    pos := pos + 1; -- skip the opening bracket
    Standard_Parse_Numbers.Parse(s,pos,re_cst);
    Standard_Parse_Numbers.Parse(s,pos,im_cst);
    cst := Standard_Complex_Numbers.Create(re_cst,im_cst);
    size := Size_Count(s,t);
    if vrblvl > 0 then
      put("constant : "); put(cst); new_line;
      put("size of the series : "); put(size,1); new_line;
    end if;
    c := new Standard_Complex_Vectors.Vector(0..size);
    p := new Standard_Floating_Vectors.Vector(1..size);
    c(0) := cst;
    for i in 1..size loop
      while pos <= s'last loop -- find the opening round bracket
        exit when s(pos) = '(';
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip the opening bracket
      Standard_Parse_Numbers.Parse(s,pos,re_cst);
      Standard_Parse_Numbers.Parse(s,pos,im_cst);
      c(i) := Standard_Complex_Numbers.Create(re_cst,im_cst);
      while pos <= s'last loop -- find t
        exit when s(pos) = 't';
        pos := pos + 1;
      end loop;
      pos := pos + 1; -- skip t
      while pos <= s'last loop -- skip * and spaces ...
        exit when ((s(pos) /= ' ') and (s(pos) /= '*'));
        pos := pos + 1;
      end loop;
      Standard_Parse_Numbers.Parse(s,pos,p(i));
    end loop;
  end parse_string;

  procedure get ( c : out Standard_Complex_Vectors.Vector;
                  p : out Standard_Floating_Vectors.Vector;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Series_IO.get ...");
    end if;
    get(standard_input,c,p,t);
  end get;

  procedure get ( file : in file_type;
                  c : out Standard_Complex_Vectors.Vector;
                  p : out Standard_Floating_Vectors.Vector;
                  t : in character := 't'; vrblvl : in integer32 := 0 ) is

    ch : character := ' ';
    re_cst,im_cst : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in Real_Powered_Series_IO.get ...");
    end if;
    while not end_of_file(file) loop -- find opening bracket
      get(file,ch);
      exit when (ch = '(');
    end loop;
    get(file,ch); -- skip opening bracket
    Standard_Parse_Numbers.parse(file,ch,re_cst);
    Standard_Parse_Numbers.parse(file,ch,im_cst);
    c(0) := Standard_Complex_Numbers.create(re_cst,im_cst);
    if vrblvl > 0 then
      put("constant : "); put(c(0)); new_line;
    end if;
    for k in p'range loop
      while not end_of_file(file) loop -- find opening bracket
        get(file,ch);
        exit when (ch = '(');
      end loop;
      get(file,ch); -- skip opening bracket
      Standard_Parse_Numbers.parse(file,ch,re_cst);
      Standard_Parse_Numbers.parse(file,ch,im_cst);
      c(k) := Standard_Complex_Numbers.create(re_cst,im_cst);
      while not end_of_file(file) loop -- find t
        get(file,ch);
        exit when (ch = t);
      end loop;
      get(file,ch); -- skip t
      while not end_of_file(file) loop -- skip spaces and *
        get(file,ch);
        exit when ((ch /= ' ') and (ch /= '*'));
      end loop;
      Standard_Parse_Numbers.Parse(file,ch,p(k));
    end loop;
  end get;

end Real_Powered_Series_IO;
