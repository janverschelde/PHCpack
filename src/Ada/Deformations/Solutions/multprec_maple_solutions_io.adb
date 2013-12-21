with Symbol_Table;                       use Symbol_Table;
with Symbol_Table_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with generic_lists;
with Lists_of_Symbols;                   use Lists_of_Symbols;

package body Multprec_Maple_Solutions_io is

  package Coefficient_Lists is new Generic_Lists(Floating_Number);

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
        declare
          re : Floating_Number := REAL_PART(ls.v(i));
          im : Floating_Number := IMAG_PART(ls.v(i));
        begin
          put(file,re); put(file,", "); put(file,im);
          Clear(re); Clear(im);
        end;
        if i < ls.n then
          put_line(file,"],");
        elsif not Is_Null(tmp) then
          put_line(file,"]],");
        else
          put_line(file,"]]];");
        end if;
      end loop;
    end loop;
  end put_pairs;

  function Create ( n : natural32; cl : Coefficient_Lists.List;
                    sl : Symbol_List ) return Solution is

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
        rp,ip : Floating_Number;  -- real and imaginary part of coefficient
        sb : constant Symbol := Head_Of(sl_tmp);
        ind : constant natural32 := Symbol_Table.Get(sb);
      begin
        rp := Coefficient_Lists.Head_Of(cl_tmp);
        cl_tmp := Coefficient_Lists.Tail_Of(cl_tmp);
        ip := Coefficient_Lists.Head_Of(cl_tmp);
        cl_tmp := Coefficient_Lists.Tail_Of(cl_tmp);
       -- put("position of the symbol "); Symbol_Table_io.put(sb);
       -- put(" is "); put(ind,1); new_line;
        s.v(integer32(ind)) := Create(rp,ip);
      end;
      sl_tmp := Tail_Of(sl_tmp);
    end loop;
    s.m := 1;
    s.t := Create(integer(0));
    s.err := Create(integer(0));
    s.rco := Create(integer(1));
    s.res := Create(integer(0));
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

  procedure Write_Complex_Number
              ( file : in file_type; x : in Complex_Number )  is

    re : Floating_Number := REAL_PART(x);
    im : Floating_Number := IMAG_PART(x);
    abs_im : constant Floating_Number := AbsVal(im);

  begin
    put(file,re);
    if im > 0.0 or Equal(im,0.0)
     then put(file," + ");
     else put(file," - ");
    end if;
    put(file,abs_im);
    put(file,"*I");
    Clear(re); Clear(im);
  end Write_Complex_Number;

  procedure put ( file : in file_type; s : in Solution ) is
  begin
    put(file,"time = ");
    Write_Complex_Number(file,s.t);
    put_line(file,",");
    put(file,"  multiplicity = "); put(file,s.m,1);
    put_line(file,",");
    for i in s.v'range loop
      put(file,"  ");
      Symbol_Table_io.put(file,Symbol_Table.Get(natural32(i)));
      put(file," = ");
      Write_Complex_Number(file,s.v(i));
      put_line(file,",");
    end loop;
    put(file,"  err ="); put(file,s.err,3); put(file,",");
    put(file,"  rco ="); put(file,s.rco,3); put(file,",");
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
--    re : Floating_Number;
--
--  begin
--    res.t := s.v(1);
--    re := REAL_PART(s.v(2));
--    res.m := natural(Trunc(re));
--    Clear(re);
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
       then put_line(file,"],");
       else put_line(file,"]];");
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
  --   Important: the character on input is ignored!

  begin
    loop
      get(file,c);
      exit when (c /= ' ');
    end loop;
  end Skip_Spaces;

  procedure Append_Zero ( first,last : in out Coefficient_Lists.List ) is

  -- DESCRIPTION :
  --   Appends zero to the list first,
  --   where last is a pointer to the last element in the list.

    zero : constant Floating_Number := Create(integer(0));

  begin
    Coefficient_Lists.Append(first,last,zero);
  end Append_Zero;

  procedure Append_Float ( first,last : in out Coefficient_Lists.List;
                           f : Floating_Number ) is

  -- DESCRIPTION :
  --   Appends a copy of the floating-number f to the list first,
  --   where last is a pointer to the last element in the list.

    nf : Floating_Number;

  begin
    Copy(f,nf);
    Coefficient_Lists.Append(first,last,nf);
  end Append_Float;

  procedure get ( file : in file_type; ls : out Link_to_Solution ) is

    cnt : natural32 := 0;
    c : character;
    f : Floating_Number := Create(integer(0));
    cl,cl_last : Coefficient_Lists.List;
    sl,sl_last : Symbol_List;
    class : natural32;
    m : natural32 := 1;
    t_re,err,rco,res : Floating_Number := Create(integer(0));
    t : Complex_Number := Create(integer(0));
    minus : boolean;

  begin
    loop                              -- read opening bracket
      get(file,c);
      exit when (c = '[');
    end loop;
    loop
      c := ' ';
      declare
        sb : Symbol;
      begin
        Symbol_Table_io.get(file,sb,'=');
       -- put("Symbol read : "); Symbol_Table_io.put(sb); new_line;
        class := Classify_Symbol(sb);
        if class = 0 then
          cnt := cnt + 1;
          Append(sl,sl_last,sb);
        end if;
      end;
      Skip_Spaces(file,c);      -- skip the "="
     -- put_line("Character after reading the symbol : " & c);
      if c = '=' then c := ' '; end if;
      if c = '-' then minus := true; else minus := false; end if;
      Clear(f); get(file,f,c); -- put(f);
     -- put_line("  Terminating symbol : " & c);
     -- put("The class of the symbol is : "); put(class,1); new_line;
      if c = 'I' then
        f := Create(1.0);
        if minus then Min(f); end if;
        c := '*';
      end if;
      case class is
        when 1 => if c = '*' then
                    Clear(t);
                    t := Create(t_re,f);
                    -- put("created pure imaginary t = ");
                    -- put(f); put_line("*I");
                  else 
                     Copy(f,t_re);
                  end if;
        when 2 => m := natural32(Round(f));
                 -- put("read multiplicity : "); put(f); new_line;
                 -- put("read multiplicity : "); put(m,1); new_line;
        when 3 => Copy(f,err);
                 -- put("read the error : "); put(f); new_line;
        when 4 => Copy(f,rco);
                 -- put("read rco : "); put(f); new_line;
        when 5 => Copy(f,res);
                 -- put("read res : "); put(f); new_line;
        when others => if c = '*'  -- case of pure imaginary number
                        then Append_Zero(cl,cl_last);
                       end if;
                       Append_Float(cl,cl_last,f);
      end case;
      if c = ' ' then Skip_Spaces(file,c); end if;
     -- put_line("Character after first float read : " & c);
      Clear(f);
      case c is
        when '+' => c := ' ';
                    get(file,f,c); -- put(" + "); put(f); new_line;
                    if class = 0
                     then Append_Float(cl,cl_last,f);
                     elsif class = 1
                         then t := Create(t_re,f);
                              Clear(f);
                    end if;
        when '-' => c := ' ';
                    get(file,f,c); Min(f); -- put(" - "); put(f); new_line;
                    if class = 0 then
                      Append_Float(cl,cl_last,f);
                    elsif class = 1 then
                      t := Create(t_re,f); Clear(f);
                    end if;
        when others => null; -- put_line("no + and no -, others case...");
      end case;
     -- put_line("After the case with c = " & c);
      if (c = ',') or (c = ']')
       then if class = 0 then
              Append_Zero(cl,cl_last);
             -- put_line("Appended a zero imaginary part...");
            elsif class = 1 then
              t := Create(t_re);
            end if;
       else loop
              Skip_Spaces(file,c);
              exit when (c = ',' or c = ']');
            end loop;
      end if;
     -- put_line("Before the exit test with c = " & c);
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
    ls.t := t; ls.m := integer32(m);
    ls.err := err; ls.rco := rco; ls.res := res;
    Coefficient_Lists.Clear(cl);
    Clear(sl);
  end get;

  procedure get ( file : in file_type; sols : out Solution_List ) is

    res,res_last : Solution_List;
    c : character;
    cnt : natural32 := 0;

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

end Multprec_Maple_Solutions_io;
