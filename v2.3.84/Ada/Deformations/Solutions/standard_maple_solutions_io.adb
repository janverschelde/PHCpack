with Symbol_Table;                       use Symbol_Table;
with Symbol_Table_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with generic_lists;
with Lists_of_Symbols;                   use Lists_of_Symbols;

package body Standard_Maple_Solutions_io is

  package Coefficient_Lists is new Generic_Lists(double_float);

  procedure put_pairs ( sols : in Solution_List ) is
  begin
    put_pairs(Standard_Output,sols);
  end put_pairs;

  procedure put_pairs ( file : in file_type; sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution; -- := Head_Of(sols);
    first : boolean := true;

  begin
   -- Write_Symbols(file,ls.n);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      tmp := Tail_Of(tmp);
      for i in ls.v'range loop
        if i = 1 then
          if first
           then put(file,"[[["); first := false;
           else put(file," [[");
          end if;
        else
          put(file,"  [");
        end if;
        put(file,REAL_PART(ls.v(i))); put(file," , ");
        put(file,IMAG_PART(ls.v(i)));
        if i < ls.n then
          put_line(file," ] ,");
        elsif not Is_Null(tmp) then
          put_line(file," ] ] ,");
          else put_line(file," ] ] ];");
        end if;
      end loop;
    end loop;
  end put_pairs;

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

  procedure put ( s : in Solution ) is
  begin
    put(Standard_Output,s);
  end put;

  procedure put ( sols : in Solution_List ) is
  begin
    put(Standard_Output,sols);
  end put;

 -- function Write ( s : Solution ) return string is
 -- begin
 --   return "";
 -- end Write;

 -- function Write ( sols : Solution_List ) return string is
 -- begin
 --   return "";
 -- end Write;

  procedure Write_Complex_Number
              ( file : in file_type; x : in Complex_Number )  is
  begin
    put(file,REAL_PART(x));
    if IMAG_PART(x) >= 0.0
     then put(file,"  + ");
     else put(file,"  - ");
    end if;
    put(file,Abs(IMAG_PART(x)));
    put(file,"*I");
  end Write_Complex_Number;

  procedure put ( file : in file_type; s : in Solution ) is
  begin
    put(file,"time = ");
    Write_Complex_Number(file,s.t);
    put_line(file," ,");
    put(file,"  multiplicity = "); put(file,s.m,1);
    put_line(file," ,");
    for i in s.v'range loop
      put(file,"  ");
      Symbol_Table_io.put(file,Symbol_Table.Get(natural32(i)));
      put(file," = ");
      Write_Complex_Number(file,s.v(i));
      put_line(file," ,");
    end loop;
    put(file,"  err ="); put(file,s.err,3); put(file," ,");
    put(file,"  rco ="); put(file,s.rco,3); put(file," ,");
    put(file,"  res ="); put(file,s.res,3);
  end put;

--  function Contract ( s : Solution ) return Solution is
--
--  -- DESCRIPTION :
--  --   The fields t and m are in s.v(1) and s.v(2),
--  --   whereas the last three numbers in s.v are the
--  --   error, rcond and residual of the solution on return.
--
--    res : Solution(s.n-5);
--
--  begin
--    res.t := s.v(1);
--    res.m := natural(REAL_PART(s.v(2)));
--    for i in res.v'range loop
--      res.v(i) := s.v(i+2);
--    end loop;
--    res.err := REAL_PART(s.v(s.v'last-2));
--    res.rco := REAL_PART(s.v(s.v'last-1));
--    res.res := REAL_PART(s.v(s.v'last));
--    return res;
--  end Contract;

  procedure put ( file : in file_type; sols : in Solution_List ) is

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    first : boolean := true;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if first 
       then put(file,"[["); first := false;
       else put(file," [");
      end if;
      put(file,ls.all);
      tmp := Tail_Of(tmp);
      if not Is_Null(tmp)
       then put_line(file," ] ,");
       else put_line(file," ] ];");
      end if;
    end loop;
  end put;

  procedure get ( sols : out Solution_List ) is
  begin
    get(Standard_Input,sols);
  end get;

  procedure Skip_Spaces ( file : in file_type; c : in out character ) is

  -- DESCRIPTION :
  --   Reads next character and continues reading if c is space.

  begin
    loop
      get(file,c);
      exit when (c /= ' ');
    end loop;
  end Skip_Spaces;

  procedure get ( file : in file_type; ls : out Link_to_Solution ) is

    cnt : natural32 := 0;
    c : character := '+';
    f : double_float := 0.0;
    cl,cl_last : Coefficient_Lists.List;
    sl,sl_last : Symbol_List;
    class : natural32;
    m : integer32 := 1;
    t_re,err,rco,res : double_float := 0.0;
    t : Complex_Number := Create(0.0);

  begin
    loop                              -- read opening bracket
      get(file,c);
      exit when (c = '[');
    end loop;
    loop
      declare
        sb : Symbol;
      begin
        Symbol_Table_io.get(file,sb);
        class := Classify_Symbol(sb);
        if class = 0 then
          cnt := cnt + 1;
          Append(sl,sl_last,sb);
        end if;
      end;
      Skip_Spaces(file,c);            -- skip the "="
      get(file,f); -- put(f);
      case class is
        when 1 => t_re := f;
        when 2 => m := integer32(f);
        when 3 => err := f;
        when 4 => rco := f;
        when 5 => res := f;
        when others => Coefficient_Lists.Append(cl,cl_last,f);
      end case;
      Skip_Spaces(file,c);
      case c is
        when '+' => get(file,f);          -- put("  "); put(f); new_line;
                    if class = 0 then
                      Coefficient_Lists.Append(cl,cl_last,f);
                    elsif class = 1 then
                      t := Create(t_re,f);
                    end if;
        when '-' => get(file,f); f := -f; -- put("  "); put(f); new_line;
                    if class = 0 then
                      Coefficient_Lists.Append(cl,cl_last,f);
                    elsif class = 1 then
                      t := Create(t_re,f);
                    end if;
        when others => null;
      end case;
      if (c = ',') or (c = ']') then
        if class /= 2
         then Coefficient_Lists.Append(cl,cl_last,0.0);
        end if;
      else
        loop
          Skip_Spaces(file,c);
          exit when (c = ',' or c = ']');
        end loop;
      end if;
      exit when (c = ']');
    end loop;
    if c /= ']' then
      loop                        -- read closing bracket
        get(file,c);
        exit when (c = ']');
      end loop;
    end if;
    if Symbol_Table.Number = 0
     then Create_Symbol_Table(sl);
    end if;
    ls := new Solution'(Create(cnt,cl,sl));
    ls.t := t; ls.m := m;
    ls.err := err; ls.rco := rco; ls.res := res;
    Coefficient_Lists.Clear(cl);
    Clear(sl);
  end get;

  procedure get ( file : in file_type; sols : out Solution_List ) is

    res,res_last : Solution_List;
    c : character;
    cnt : natural := 0;

  begin
    loop
      get(file,c);
      exit when (c = '[');  
    end loop;
    loop
      cnt := cnt + 1;
     -- put("reading solution "); put(cnt,1); put_line("...");
      declare
        ls : Link_to_Solution;
      begin
        get(file,ls);
        Append(res,res_last,ls);
      end;
      loop
        get(file,c);
        exit when (c /= ' ');
      end loop;
      exit when (c /= ',');
    end loop;
    sols := res;
  end get;

end Standard_Maple_Solutions_io;
