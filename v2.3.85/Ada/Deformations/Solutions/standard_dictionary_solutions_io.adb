with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with generic_lists;
with Lists_of_Symbols;                   use Lists_of_Symbols;

package body Standard_Dictionary_Solutions_io is

  package Coefficient_Lists is new Generic_Lists(double_float);

-- AUXILIARIES :

  procedure Write_Complex_Number
              ( file : in file_type; x : in Complex_Number ) is

  -- DESCRIPTION :
  --   A complex number x = a + i*b in Python format is written to
  --   file in the format complex(a,b).

  begin
    put(file,REAL_PART(x));
    if IMAG_PART(x) >= 0.0 
     then put(file," + ");
    end if;
    put(file,IMAG_PART(x)); put(file,"*1j");
  end Write_Complex_Number;

  function Create ( n : natural32; cl : Coefficient_Lists.List;
                    sl : Symbol_List  ) return Solution is

  -- DESCRIPTION :
  --   Returns a solution of dimension n with coefficients in cl.
  --   The coefficients are pairs of floats, storing the real and
  --   imaginary parts of the complex coefficients one after the other.
  --   The list of symbols indicates the relative position of the
  --   coefficients in the symbol table.

  -- REQUIRED : Length_Of(cl) = 2*Length_Of(sl).

    s : Solution(integer32(n));
    cl_tmp : Coefficient_Lists.List := cl;
    sl_tmp : Symbol_List := sl;

  begin
    while not Is_Null(sl_tmp) loop
      declare
        rp,ip : double_float;  -- real and imaginary part of coefficient
        sb : constant Symbol := Head_Of(sl_tmp);
        ind : constant integer32 := integer32(Symbol_Table.Get(sb));
      begin
        rp := Coefficient_Lists.Head_Of(cl_tmp);
        cl_tmp := Coefficient_Lists.Tail_Of(cl_tmp);
        ip := Coefficient_Lists.Head_Of(cl_tmp);
        cl_tmp := Coefficient_Lists.Tail_Of(cl_tmp);
       -- put("position of the symbol "); Symbol_Table_io.put(sb);
       -- put(" is "); put(ind,1); new_line;
        s.v(ind) := Create(rp,ip);
      end;
      sl_tmp := Tail_Of(sl_tmp);
    end loop;
    s.m := 1;
    s.t := Create(0.0);
    s.err := 0.0; s.rco := 1.0; s.res := 0.0;
    return s;
  end Create;

-- TARGET ROUTINES :

  procedure put ( s : in Solution ) is
  begin
    put(standard_output,s);
  end put;
 
  procedure put ( file : in file_type; s : in Solution ) is
  begin
    put(file,'{');
    put(file,"'time':"); Write_Complex_Number(file,s.t); put(file,",");
    put(file,"'multiplicity':"); put(file,s.m,1); put(file,",");
    for i in s.v'range loop
      put(file,"'");
      Symbol_Table_io.put(file,Symbol_Table.Get(natural32(i)));
      put(file,"':");
      Write_Complex_Number(file,s.v(i)); put(file,",");
    end loop;
    put(file,"'err':"); put(file,s.err,3); put(file,",");
    put(file,"'rco':"); put(file,s.rco,3); put(file,",");
    put(file,"'res':"); put(file,s.res,3);
    put(file,'}');
  end put;

  procedure put ( sols : in Solution_List ) is
  begin
    put(standard_output,sols);
  end put;

  procedure put ( file : in file_type; sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    put_line(file,"[");
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      put(file,ls.all);
      tmp := Tail_Of(tmp);
      if not Is_Null(tmp)
       then put_line(file,",");
      end if;
    end loop;
    new_line(file);
    put_line(file,"]");
  end put;

--  function Write ( s : Solution ) return string is
--  begin
--    return "";
--  end Write;

--  function Write ( sols : Solution_List ) return string is
--  begin
--    return "";
--  end Write;

  procedure read_symbol ( file : in file_type; sb : out Symbol ) is

  -- DESCRIPTION :
  --   Skips all characters on file till the first opening quote
  --   and stores all characters in sb till the closing quote.

    p : integer := sb'first-1;
    c : character;

  begin
    sb := (sb'range => ' ');
    loop
      get(file,c);
      exit when (c = ''');
    end loop;
    loop
      get(file,c);
      exit when (c = ''');
      p := p + 1;
      sb(p) := c;
    end loop;
  end read_symbol;

--  procedure read_float ( file : in file_type; f : out double_float;
--                         fail : out boolean ) is
--
--  -- DESCRIPTION :
--  --   Reads a float from file, skipping the line if an exception happens.
--
--  begin
--    put_line("inside read_float...");
--    f := 0.0;
--    get(file,f);
--  exception
--    when others
--      => put_line("an exception happened when reading a float");
--         skip_line(file); f := 0.0; fail := true;
--  end read_float;

  procedure get ( file : in file_type; ls : out Link_to_Solution ) is

    cnt : natural32 := 0;
    c : character;
    sb : Symbol;
    sl,sl_last : Symbol_List;
    class : natural32;
    m : integer32 := 1;
    f,t_re,err,rco,res : double_float := 0.0;
    t : Complex_Number := Create(0.0);
    cl,cl_last : Coefficient_Lists.List;
   -- fail : boolean;

  begin
    loop get(file,c); exit when (c = '{'); end loop;
    loop
      read_symbol(file,sb);
      class := Classify_Symbol(sb);
      if class = 0 then
        cnt := cnt + 1;
        Append(sl,sl_last,sb);
      end if;
      loop get(file,c); exit when (c = ':'); end loop;
      get(file,f); -- read_float(file,f,fail); --get(file,f);
      case class is
         when 1 => t_re := f;
         when 2 => m := integer32(f);
         when 3 => err := f;
         when 4 => rco := f;
         when 5 => res := f;
         when others => Coefficient_Lists.Append(cl,cl_last,f);
      end case;
      loop get(file,c); exit when (c /= ' '); end loop;
      case c is
        when '+' => -- read_float(file,f,fail);
                    get(file,f);          -- put("  "); put(f); new_line;
                    if class = 0 then
                      Coefficient_Lists.Append(cl,cl_last,f);
                    elsif class = 1 then
                      t := Create(t_re,f);
                    end if;
        when '-' => -- read_float(file,f,fail);
                    get(file,f); f := -f; -- put("  "); put(f); new_line;
                    if class = 0 then
                      Coefficient_Lists.Append(cl,cl_last,f);
                    elsif class = 1 then
                      t := Create(t_re,f);
                    end if;
        when others => null;
      end case;
      if (c = ',') or (c = '}') then
        if class /= 2
         then Coefficient_Lists.Append(cl,cl_last,0.0);
        end if;
      else
        loop 
          get(file,c);
          exit when ((c = ',') or (c = '}'));
        end loop;
      end if;
      exit when (c = '}');
    end loop;
    if Symbol_Table.Number = 0
     then Create_Symbol_Table(sl);
    end if;
    ls := new Solution'(Create(cnt,cl,sl));
    ls.t := t; ls.m := m;
    ls.err := err; ls.rco := rco; ls.res := res;
    Coefficient_Lists.Clear(cl);
    Clear(sl);
  end get;

  procedure get ( sols : out Solution_List ) is
  begin
    get(standard_input,sols);
  end get;

  procedure get ( file : in file_type; sols : out Solution_List ) is

    c : character;
    cnt : natural32 := 0;
    res,res_last : Solution_List;

  begin
    loop
      get(file,c);
      exit when (c = '[');
    end loop;
    skip_line(file);
    loop
      cnt := cnt + 1;
      declare
        ls : Link_to_Solution;
      begin
       -- put("reading solution "); put(cnt,1); put_line(" ...");
        get(file,ls);
        Append(res,res_last,ls);
      end;
      loop
        get(file,c);
        exit when (c /= ' ');
      end loop;
      exit when (c /= ',');  -- no comma => end of solution list
    end loop;
    sols := res;
  end get;

end Standard_Dictionary_Solutions_io;
