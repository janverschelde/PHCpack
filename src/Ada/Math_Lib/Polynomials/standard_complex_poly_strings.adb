with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Strings_and_Numbers;                use Strings_and_Numbers;
with Standard_Parse_Numbers;             use Standard_Parse_Numbers;
with Symbol_Table_io;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions;
with Standard_Complex_Term_Lists_io;     use Standard_Complex_Term_Lists_io;

package body Standard_Complex_Poly_Strings is

-- STRING MANIPULATORS :

  procedure Read_Exponent
              ( s : in string; k : in out integer; e : out natural32 ) is

    exp : string(1..16);
    cnt : natural := 0;

  begin
    while Convert(s(k)) < 10 loop
      cnt := cnt+1;
      exp(cnt) := s(k);
      k := k+1;
    end loop;
    e := Convert(exp(1..cnt));
  end Read_Exponent;

  function Delimiters ( n : natural32; s : string )
                      return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Searches the string s for n semicolons and returns the positions
  --   of those semicolons.

    res : Standard_Natural_Vectors.Vector(1..integer32(n));
    ind : integer32 := 1;

  begin
    for i in s'range loop
      if s(i) = ';' then
        res(ind) := natural32(i);
        ind := ind+1;
      end if;
    end loop;
    return res;
  end Delimiters;

  function Concat_Symbol0 ( s : string; sb : Symbol ) return string is

    nb : natural := 0;

  begin
    for i in sb'range loop
      if sb(i) = ' '
       then nb := i-1; exit;
      end if;
    end loop;
    declare
      res : string(1..s'last+nb);
    begin
      res(s'range) := s;
      for i in 1..nb loop
        res(s'last+i) := sb(i);
      end loop;
      return res;
    end;
  end Concat_Symbol0;

  function Concat_Symbol1 ( s : string; sb : Symbol ) return string is

    nb : natural := 0;

  begin
    for i in sb'range loop
      if sb(i) = ' '
       then nb := i-1; exit;
      end if;
    end loop;
    declare
      res : string(1..s'last+1+nb);
    begin
      res(s'range) := s;
      res(s'last+1) := '*';
      for i in 1..nb loop
        res(s'last+1+i) := sb(i);
      end loop;
      return res;
    end;
  end Concat_Symbol1;

-- AUXILIARIES FOR INPUT PARSING OPERATORS :

  procedure Parse_Term
              ( s : in string; bc : in out integer32; p : in out integer;
                n : in natural32; termp : in out Poly );
  procedure Parse_Term
              ( s : in string; bc : in out integer32; p : in out integer;
                n : in natural32; termp,termp_last : in out Term_List );

  -- DESCRIPTION :
  --   Parses the string for a term, starting at s(p).

  -- ON ENTRY :
  --   s        string of characters to contain a polynomial;
  --   bc       counts #open - #closed brackets;
  --   p        current position in the string s;
  --   n        number of variables;
  --   termp    accumulating polynomial for the term.

  -- ON RETURN :
  --   bc       updated bracket counter;
  --   p        updated position in the string;
  --   termp    resulting term read from string.

  procedure Parse_Factor
              ( s : in string; bc : in out integer32; p : in out integer;
                n : in natural32; d : in out Degrees; pb : in out Poly );
  procedure Parse_Factor
              ( s : in string; bc : in out integer32; p : in out integer;
                n : in natural32; d : in out Degrees;
                pb,pb_last : in out Term_List );

  -- DESCRIPTION :
  --   Parses the string s for a factor, starting at s(p).

  -- ON ENTRY :
  --   s        string of characters to contain a polynomial;
  --   bc       counts #open - #closed brackets;
  --   p        current position in the string s;
  --   n        number of variables;
  --   d        current degrees of the factor;
  --   pb       accumulating polynomial for the factor.

  -- ON RETURN :
  --   bc       updated bracket counter;
  --   p        updated position in the string;
  --   pb       resulting factor read from string.

  procedure Parse_Power_Factor
              ( s : in string; k : in out integer; p : in out Poly ) is

  -- DESCRIPTION :
  --   Reads the exponent of a factor stored in s and returns
  --   in p the polynomial p raised to the read exponent.

  -- ON ENTRY :
  --   s        string to contain a polynomial
  --   k        s(k) equals the exponentiation symbol '^';
  --   p        current factor.

  -- ON RETURN :
  --   k        updated position in the string;
  --   p        p raised to the read exponent.
 
    d : natural32 := 0;
    bracket : boolean := false;

  begin
    k := k + 1;                -- skip the '^'
    if k > s'last
     then return;
    end if;
    if s(k) = '('
     then k := k + 1; bracket := true;
    end if;
    Read_Exponent(s,k,d);
    if bracket then
      Skip_Spaces_and_CR(s,k);
      if k > s'last
       then return;
      end if;
      if s(k) = ')'
       then k := k + 1;        -- skip closing bracket
       else raise BAD_BRACKET; -- no closing bracket found
      end if;
    end if;
    Pow(p,d);
    Skip_Spaces_and_CR(s,k);   -- skip last digit of exponent 
  exception
    when others =>
      put("Exception raised at character "); put(integer32(k),1);
      put_line(" of " & s & " in Parse_Power_Factor."); raise;
  end Parse_Power_Factor;

  procedure Mul ( p,p_last : in out Term_List; q : in Term_List ) is

  -- DESCRIPTION :
  --   Multiplies every term in p by every term in q.

    ptmp : Term_List := p;
    qtmp : Term_List;
    pt,qt : Term;

  begin
    while not Is_Null(ptmp) loop
      pt := Head_Of(ptmp);
      qtmp := q;
      while not Is_Null(qtmp) loop
        qt := Head_Of(qtmp);
        pt.cf := pt.cf*qt.cf;
        for i in pt.dg'range loop
          pt.dg(i) := pt.dg(i) + qt.dg(i);
        end loop;
        qtmp := Tail_Of(qtmp);
      end loop;
      ptmp := Tail_Of(ptmp);
    end loop;
  end Mul;

  procedure Pow ( p,p_last : in out Term_List; d : in natural32 ) is

  -- DESCRIPTION :
  --   Raises the terms in the list p to the power d.
  --   For lists of terms, this procedure is just a stub because
  --   the input format does not support powers of polynomials.

  begin
    null;
  end Pow;

  procedure Min ( p : in out Term_List ) is

  -- DESCRIPTION :
  --   Flips the sign of all terms in the list p.

    tmp : Term_List;
    t : Term;

  begin
    while not Is_Null(tmp) loop
      t := Head_Of(tmp);
      t.cf := -t.cf;
      tmp := Tail_Of(tmp);
    end loop;
  end Min;

  procedure Sub ( p,p_last : in out Term_List; q : in Term_List ) is

  -- DESCRIPTION :
  --   Flips the sign of the coefficients of the terms in q
  --   before appending them to the list p.

    tmp : Term_List := q;
    t : Term;

  begin 
    while not Is_Null(tmp) loop
      t := Head_Of(tmp);
      declare
        mint : Term;
      begin
        mint.cf := -t.cf;
        mint.dg := t.dg;
        Merge_Append(p,p_last,mint);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Sub;

  procedure Parse_Power_Factor
              ( s : in string; k : in out integer;
                p,p_last : in out Term_List ) is

  -- DESCRIPTION :
  --   Reads the exponent of a factor stored in s and returns
  --   in p the polynomial p raised to the read exponent.

  -- ON ENTRY :
  --   s        string to contain a polynomial
  --   k        s(k) equals the exponentiation symbol '^';
  --   p        current factor.

  -- ON RETURN :
  --   k        updated position in the string;
  --   p        p raised to the read exponent.
 
    d : natural32 := 0;
    bracket : boolean := false;

  begin
   -- put_line("In Parse_Power_Factor ...");
    k := k + 1;                -- skip the '^'
    if k > s'last
     then return;
    end if;
    if s(k) = '('
     then k := k + 1; bracket := true;
    end if;
    Read_Exponent(s,k,d);
    if bracket then
      Skip_Spaces_and_CR(s,k);
      if k > s'last
       then return;
      end if;
      if s(k) = ')'
       then k := k + 1;        -- skip closing bracket
       else raise BAD_BRACKET; -- no closing bracket found
      end if;
    end if;
    Pow(p,p_last,d);
    Skip_Spaces_and_CR(s,k);   -- skip last digit of exponent 
  exception
    when others =>
      put("Exception raised at character "); put(integer32(k),1);
      put_line(" of " & s & " in Parse_Power_Factor."); raise;
  end Parse_Power_Factor;

  procedure Parse_Polynomial
              ( s : in string; bc : in out integer32; k : in out integer;
                n : in natural32; p : in out Poly ) is

  -- DESCRIPTION :
  --   Main parsing routine for a polynomial from a string,
  --   keeps track of the bracket count in bc.

    oper : character;
    term,res,acc : Poly := Null_Poly;

  begin
    oper := '+';
    Skip_Spaces_and_CR(s,k);
    if k > s'last
     then return;
    end if;
    if s(k) = '-'
     then oper := '-';
    end if;                         -- the first term can have no sign
    Parse_Term(s,bc,k,n,res);       -- therefore read it first
    loop
      case s(k) is
        when ' ' | ASCII.CR | ASCII.LF => k := k + 1;   -- skip blanks
        when '+' | '-' =>
          oper := s(k);
          Parse_Term(s,bc,k,n,term);
          Add(res,term); Clear(term);
        when delimiter => 
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' =>
          bc := bc + 1;
          k := k + 1;
          Parse_Polynomial(s(k..s'last),bc,k,n,term);
          Skip_Spaces_and_CR(s,k);
          if s(k) = '(' -- or s(k) = ')'
           then raise BAD_BRACKET;
          end if;
          if s(k) = '^' then
            Parse_Power_Factor(s,k,term);
          elsif s(k) = '*' then
            if s(k+1) = '*' then
              k := k+1;
              Parse_Power_Factor(s,k,term);
            end if;
          end if;
          case oper is
            when '+' => Add(acc,res); Clear(res); Copy(term,res);
            when '-' => Add(acc,res); Clear(res); Copy(term,res); Min(res);
            when '*' => Mul(res,term);
            when others => raise ILLEGAL_OPERATION;
          end case;
          Clear(term);
        when ')' =>
          bc := bc - 1;
          k := k + 1;
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '*' =>
         -- if res = Null_Poly then    -- the case (0.0)*(2 + x)
         --   raise ILLEGAL_CHARACTER; -- should not raise exception
         -- else -- the case " ) * " :
            oper := s(k); k := k + 1;  -- skip '*'
            Parse_Term(s,bc,k,n,term);
            Skip_Spaces_and_CR(s,k);
            if s(k) = '^' then
              Parse_Power_Factor(s,k,term);
            elsif s(k) = '*' then
              if s(k+1) = '*' then
                k := k+1;
                Parse_Power_Factor(s,k,term);
              end if;
            end if;
            if s(k) /= '(' then
              case oper is
                when '+' => Add(res,term);
                when '-' => Sub(res,term);
                when '*' => Mul(res,term);
                when others => raise ILLEGAL_OPERATION;
              end case;
            end if;
            Clear(term);
         -- end if;
        when '^' =>
          if res = Null_Poly
           then raise ILLEGAL_CHARACTER;
           else Parse_Power_Factor(s,k,res);
          end if;
        when others => raise ILLEGAL_CHARACTER;
      end case;
    end loop;
    p := acc + res;
    Clear(acc); Clear(res);
  exception
    when others =>
      put("Exception raised at character "); put(integer32(k),1);
      put_line(" of " & s & " in Parse_Polynomial."); raise;
  end Parse_Polynomial;

  procedure Parse_Polynomial
              ( s : in string; bc : in out integer32; k : in out integer;
                n : in natural32; p,p_last : in out Term_List ) is

  -- DESCRIPTION :
  --   Main parsing routine for a polynomial from a string,
  --   keeps track of the bracket count in bc.

    oper : character;
    term,term_last,res,res_last,acc,acc_last : Term_List;

  begin
    oper := '+';
    Skip_Spaces_and_CR(s,k);
    if k > s'last
     then return;
    end if;
    if s(k) = '-'
     then oper := '-';
    end if;                            -- the first term can have no sign
    Parse_Term(s,bc,k,n,res,res_last); -- therefore read it first
    loop
      case s(k) is
        when ' ' | ASCII.CR | ASCII.LF => k := k + 1;   -- skip blanks
        when '+' | '-' =>
          oper := s(k);
          Parse_Term(s,bc,k,n,term,term_last);
          if not Is_Null(term)
           then Merge_Concat(res,res_last,term); Clear(term);
          end if;
        when delimiter => 
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' =>
          bc := bc + 1;
          k := k + 1;
          Parse_Polynomial(s(k..s'last),bc,k,n,term,term_last);
          Skip_Spaces_and_CR(s,k);
          if s(k) = '(' -- or s(k) = ')'
           then raise BAD_BRACKET;
          end if;
          if s(k) = '^' then
            Parse_Power_Factor(s,k,term,term_last);
          elsif s(k) = '*' then
            if s(k+1) = '*' then
              k := k + 1;
              Parse_Power_Factor(s,k,term,term_last);
            end if;
          end if;
          case oper is
            when '+' => Merge_Concat(acc,acc_last,res); 
                        Clear(res); Copy(term,res);
            when '-' => Merge_Concat(acc,acc_last,res);
                        Clear(res); Copy(term,res); Min(res);
            when '*' => Mul(res,res_last,term);
            when others => raise ILLEGAL_OPERATION;
          end case;
          Clear(term);
        when ')' =>
          bc := bc - 1;
          k := k + 1;
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '*' =>
          -- if Is_Null(res) then        -- case like (0.0)*(2 + x)
          --   raise ILLEGAL_CHARACTER;  -- should not trigger exception
          -- else -- the case " ) * " :
            oper := s(k); k := k + 1;  -- skip '*'
            Parse_Term(s,bc,k,n,term,term_last);
            Skip_Spaces_and_CR(s,k);
            if s(k) = '^' then
              Parse_Power_Factor(s,k,term,term_last);
            elsif s(k) = '*' then
              if s(k+1) = '*' then
                k := k + 1;
                Parse_Power_Factor(s,k,term,term_last);
              end if;
            end if;
            if s(k) /= '(' then
              case oper is
                when '+' => Merge_Concat(res,res_last,term);
                when '-' => Sub(res,res_last,term);
                when '*' => Mul(res,res_last,term);
                when others => raise ILLEGAL_OPERATION;
              end case;
            end if;
            Clear(term);
          -- end if;
        when '^' =>
          if Is_Null(res)
           then raise ILLEGAL_CHARACTER;
           else Parse_Power_Factor(s,k,res,res_last);
          end if;
        when others => raise ILLEGAL_CHARACTER;
      end case;
    end loop;
    Merge_Concat(p,p_last,acc);
    Merge_Concat(p,p_last,res);
    Clear(acc); Clear(res);
  exception
    when others =>
      put("Exception raised at character "); put(integer32(k),1);
      put_line(" of " & s & " in Parse_Polynomial."); raise;
  end Parse_Polynomial;

  procedure Parse_Factor
              ( s : in string; bc : in out integer32; p : in out integer;
                n : in natural32; d : in out Degrees; pb : in out Poly ) is

    sb : symbol;
    expo : natural32;
    k : integer32;
 
  begin
    if p > s'last
     then return;
    end if;
    Skip_Spaces_and_CR(s,p);
    if s(p) = '(' then 
      bc := bc + 1;
      p := p + 1;        -- get a new symbol, skip '('
      if p > s'last
       then return;
      end if;
      Parse_Polynomial(s(p..s'last),bc,p,n,pb);
      Skip_Spaces_and_CR(s,p);
      if s(p) = '^' then
        Parse_Power_Factor(s,p,pb);
      elsif s(p) = '*' then
        if s(p+1) = '*' then
          p := p + 1;
          Parse_Power_Factor(s,p,pb);
        end if;
      end if;
      return;
    end if;
    Symbol_Table_io.Parse_Symbol(s,p,sb);
    k := integer32(Symbol_Table.Check_Symbol(n,sb));
    Skip_Spaces_and_CR(s,p);
    if s(p) = '^' then
      p := p + 1;                               -- skip the '^'
      if p > s'last
       then return;
      end if;
      Read_Exponent(s,p,expo);
      d(k) := d(k) + expo;
      Skip_Spaces_and_CR(s,p);
      if s(p) /= '*'                            -- the case x^2*...
       then return;                             -- end of factor
       else p := p + 1;                         -- skip the '*'
      end if; 
    elsif s(p) = '*' then
      p := p + 1;
      if p > s'last
       then return;
      end if;
      if s(p) = '*' then
        p := p + 1;                             -- the case " x ** expo "
        Read_Exponent(s,p,expo);
        d(k) := d(k) + expo;
        Skip_Spaces_and_CR(s,p);
        if s(p) /= '*'
         then return;                   -- end of factor
         else p := p + 1;               -- skip the '*'
        end if;
      else
        d(k) := d(k) + 1;               -- the case " x * ? "
      end if;
    else -- the case " x ?", with ? /= '*' or ' ' or '^' :
      d(k) := d(k) + 1;
      return;
    end if;
    Skip_Spaces_and_CR(s,p);
    if (s(p) = '-') or (s(p) = '+') 
     then return;
    end if;
    if Convert(s(p)) < 10
     then Parse_Term(s,bc,p,n,pb); -- the case " x * c " or " x ** c * c "
     else Parse_Factor(s,bc,p,n,d,pb); -- the case " x * y " 
    end if;
  exception
    when others =>
      put("Exception raised at character "); put(integer32(p),1);
      put_line(" of " & s & " in Parse_Factor."); raise;
  end Parse_Factor;

  procedure Parse_Factor
              ( s : in string; bc : in out integer32; p : in out integer;
                n : in natural32; d : in out Degrees;
                pb,pb_last : in out Term_List ) is

    sb : symbol;
    expo : natural32;
    k : integer32;
 
  begin
    if p > s'last
     then return;
    end if;
    Skip_Spaces_and_CR(s,p);
    if s(p) = '(' then 
      bc := bc + 1;
      p := p + 1;        -- get a new symbol, skip '('
      if p > s'last
       then return;
      end if;
      Parse_Polynomial(s(p..s'last),bc,p,n,pb,pb_last);
      Skip_Spaces_and_CR(s,p);
      if s(p) = '^' then
        Parse_Power_Factor(s,p,pb,pb_last);
      elsif s(p) = '*' then
        if s(p+1) = '*' then
          p := p + 1;
          Parse_Power_Factor(s,p,pb,pb_last);
        end if;
      end if;
      return;
    end if;
    Symbol_Table_io.Parse_Symbol(s,p,sb);
    k := integer32(Symbol_Table.Check_Symbol(n,sb));
    Skip_Spaces_and_CR(s,p);
    if s(p) = '^' then
      p := p + 1;                               -- skip the '^'
      if p > s'last
       then return;
      end if;
      Read_Exponent(s,p,expo);
      d(k) := d(k) + expo;
      Skip_Spaces_and_CR(s,p);
      if s(p) /= '*'                            -- the case x^2*...
       then return;                             -- end of factor
       else p := p + 1;                         -- skip the '*'
      end if; 
    elsif s(p) = '*' then
      p := p + 1;
      if p > s'last
       then return;
      end if;
      if s(p) = '*' then
        p := p + 1;                             -- the case " x ** expo "
        Read_Exponent(s,p,expo);
        d(k) := d(k) + expo;
        Skip_Spaces_and_CR(s,p);
        if s(p) /= '*'
         then return;                   -- end of factor
         else p := p + 1;               -- skip the '*'
        end if;
      else
        d(k) := d(k) + 1;               -- the case " x * ? "
      end if;
    else -- the case " x ?", with ? /= '*' or ' ' or '^' :
      d(k) := d(k) + 1;
      return;
    end if;
    Skip_Spaces_and_CR(s,p);
    if (s(p) = '-') or (s(p) = '+') 
     then return;
    end if;
    if Convert(s(p)) < 10
     then Parse_Term(s,bc,p,n,pb,pb_last);
           -- the case " x * c " or " x ** c * c "
     else Parse_Factor(s,bc,p,n,d,pb,pb_last); -- the case " x * y " 
    end if;
  exception
    when others =>
      put("Exception raised at character "); put(integer32(p),1);
      put_line(" of " & s & " in Parse_Factor."); raise;
  end Parse_Factor;
 
  procedure Parse_Term ( s : in string; bc : in out integer32;
                         p : in out integer;
                         n : in natural32; termp : in out Poly ) is

    c : Complex_Number := Create(0.0);
    d : Degrees := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    pb,res : Poly := Null_Poly;
    tmp : Term;
    stop_kkloop : boolean;

    procedure Collect_Factor_Polynomial is
    begin
      if pb /= Null_Poly then
        if res = Null_Poly
         then Copy(pb,res); Clear(pb);
         else Mul(res,pb); Clear(pb);
        end if;
      end if;
    end Collect_Factor_Polynomial;

  begin
    stop_kkloop := true;  -- introduced to handle cases like (0 + 2*i)
    for kk in 1..2 loop
      Parse(s,p,c); -- look for 'i' :
      Skip_Spaces_and_CR(s,p);
      if ( c = Create(0.0) ) then
        if ((s(p) = 'i') or (s(p) = 'I')) then -- the case "+ i" :
          c := Create(0.0,1.0); 
          p := p + 1;        -- skip 'i'
        elsif s(p) = '+' or s(p) = '-' then
          stop_kkloop := false;
        elsif s(p) = '*' then  -- the case 0*i
          p := p + 1;          -- skip the '*'
          Skip_Spaces_and_CR(s,p);
          if (s(p) = 'i') or (s(p) = 'I') then -- the case ".. c * i.." :
            p := p + 1;    -- skip 'i' or 'I"
          else -- the case 0*factor
            Parse_Factor(s,bc,p,n,d,pb);
            Clear(pb); return; -- ignore the factor
          end if;
        end if;
      elsif ( c = Create(-1.0) ) and then ((s(p) = 'i') or (s(p) = 'I')) then
      -- the case "- i" :
        c := Create(0.0,-1.0);
        p := p + 1;      -- skip 'i'
      elsif s(p) = '*' then -- the case ".. c *.." :
        p := p + 1;  -- skip '*'
        Skip_Spaces_and_CR(s,p);
        if (s(p) = 'i') or (s(p) = 'I') then -- the case ".. c * i.." :
          c := Create(0.0,REAL_PART(c));
          p := p + 1;    -- skip 'i'
        else                                 -- the case ".. c * x.." :
          Parse_Factor(s,bc,p,n,d,pb);
          if pb /= Null_Poly
           then Clear(res); Copy(pb,res); Clear(pb);
          end if;
        end if;
      end if;
      exit when stop_kkloop;
    end loop; -- the case ".. c ?"  will be treated in the loop
    loop
      case s(p) is
        when ' ' | ASCII.CR | ASCII.LF => p := p + 1;
        when '*' =>
          p := p + 1;
          Parse_Factor(s,bc,p,n,d,pb);
          Collect_Factor_Polynomial;
        when '+' | '-' => 
          -- if c = Create(0.0)
          --  then raise ILLEGAL_CHARACTER;
          --  else exit;
          -- end if;
          exit; -- term with zero coefficient is not invalid input
        when delimiter =>
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' => 
          if c = Create(0.0) or else c = Create(-1.0)
           then c := Create(0.0);
                exit; -- the case "+ (" or "- (" :
           else raise BAD_BRACKET;      -- the case "c  (" :
          end if;
        when ')' =>
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when others  =>
          if c = Create(0.0) then
            c := Create(1.0);
            Parse_Factor(s,bc,p,n,d,pb);
          elsif c = Create(-1.0) then
            Parse_Factor(s,bc,p,n,d,pb);
          elsif s(p) = '^' then
            Parse_Power_Factor(s,p,res);
          elsif s(p) = '*' then
            if s(p+1) = '*' then
              p := p + 1;
              Parse_Power_Factor(s,p,res);
            end if;
          else
            raise ILLEGAL_CHARACTER;
          end if;
          Collect_Factor_Polynomial;
      end case;
    end loop;
    if c /= create(0.0) then
      tmp.cf := c;
      tmp.dg := d;
      termp := create(tmp);
      Clear(tmp);
      if Number_Of_Unknowns(res) > 0
       then Mul(termp,res); Clear(res);
      end if;
    end if;
  exception
    when BAD_BRACKET =>
      put("Exception BAD_BRACKET raised at character "); put(integer32(p),1);
      put_line(" of " & s & " in Parse_Term."); 
      put_line("The character is " & s(p)); raise;
    when ILLEGAL_CHARACTER =>
      put("Exception ILLEGAL_CHARACTER raised at character ");
      put(integer32(p),1);
      put_line(" of " & s & " in Parse_Term.");
      put_line("The character is " & s(p));
      put_line("Surrounding characters are " & s(s'first..p+1));
      put_line("Surrounding characters are " & s(p-3..p+1));
      put_line("Surrounding characters are " & s(p-5..p+1));
      put_line("Surrounding characters are " & s(p-7..p+1));
      raise;
    when others =>
      put("Exception raised at character "); put(integer32(p),1);
      put_line(" of " & s & " in Parse_Term.");
      put_line("The character is " & s(p)); raise;
  end Parse_Term;
 
  procedure Parse_Term ( s : in string; bc : in out integer32;
                         p : in out integer; n : in natural32;
                         termp,termp_last : in out Term_List ) is

    zero : constant Complex_Number := Create(0.0);
    c : Complex_Number := zero;
    d : Degrees := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    pb,pb_last,res,res_last : Term_List;
    tmp : Term;

    procedure Collect_Factor_Polynomial is
    begin
      if not Is_Null(pb) then
        if Is_Null(res)
         then Copy(pb,res); Clear(pb);
         else Mul(res,res_last,pb); Clear(pb);
        end if;
      end if;
    end Collect_Factor_Polynomial;

  begin
    Parse(s,p,c); -- look for 'i' :
    Skip_Spaces_and_CR(s,p);
    if ( c = Create(0.0) ) and then ((s(p) = 'i') or (s(p) = 'I')) then
      -- the case "+ i" :
      c := Create(0.0,1.0); 
      p := p + 1;        -- skip 'i'
    elsif ( c = Create(-1.0) ) and then ((s(p) = 'i') or (s(p) = 'I')) then
      -- the case "- i" :
      c := Create(0.0,-1.0);
      p := p + 1;      -- skip 'i'
    elsif s(p) = '*' then -- the case ".. c *.." :
      Skip_Spaces_and_CR(s,p);
      p := p + 1;  -- skip '*'
      Skip_Spaces_and_CR(s,p);
      if (s(p) = 'i') or (s(p) = 'I') then -- the case ".. c * i.." :
        c := Create(0.0,REAL_PART(c));
        p := p + 1;    -- skip 'i'
      else                                 -- the case ".. c * x.." :
        Parse_Factor(s,bc,p,n,d,pb,pb_last);
        if not Is_Null(pb)
         then Clear(res); Copy(pb,res); Clear(pb);
        end if;
      end if;
    end if; -- the case ".. c ?"  will be treated in the loop
    loop
      case s(p) is
        when ' ' | ASCII.CR | ASCII.LF => p := p + 1;
        when '*' =>
          p := p + 1;
          Parse_Factor(s,bc,p,n,d,pb,pb_last);
          Collect_Factor_Polynomial;
        when '+' | '-' => 
          -- if c = Create(0.0)
          --  then raise ILLEGAL_CHARACTER;
          --  else exit;
          -- end if;
          exit; -- term with zero coefficient is not invalid input
        when delimiter =>
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' => 
          if c = Create(0.0) or else c = Create(-1.0)
           then c := Create(0.0); exit; -- the case "+ (" or "- (" :
           else raise BAD_BRACKET;      -- the case "c  (" :
          end if;
        when ')' =>
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when others  =>
          if c = Create(0.0) then
            c := Create(1.0);
            Parse_Factor(s,bc,p,n,d,pb,pb_last);
          elsif c = Create(-1.0) then
            Parse_Factor(s,bc,p,n,d,pb,pb_last);
          elsif s(p) = '^' then
            Parse_Power_Factor(s,p,res,res_last);
          elsif s(p) = '*' then
            if s(p+1) = '*' then
              p := p + 1;
              Parse_Power_Factor(s,p,res,res_last);
            end if;
          else
            raise ILLEGAL_CHARACTER;
          end if;
          Collect_Factor_Polynomial;
      end case;
    end loop;
    if not Equal(c,zero) then
      tmp.cf := c;
      tmp.dg := d;
      Merge_Append(termp,termp_last,tmp);
      Clear(tmp);
      if Length_Of(res) > 0
       then Mul(termp,termp_last,res); Clear(res);
      end if;
    end if;
  exception
    when others =>
      put("Exception raised at character "); put(integer32(p),1);
      put_line(" of " & s & " in Parse_Term."); raise;
  end Parse_Term;

-- AUXILIARIES FOR OUTPUT :

  function Write ( t : Term; lead : boolean;
                   s : Array_of_Symbols ) return string is

  -- DESCRIPTION :
  --   Returns the string representation of the term.
  --   When the term is the leading term of the polynomial,
  --   then lead is true.
  --   The array of symbols in s represent the variables.

    function Rec_Write ( k : integer32; accu : string; first : boolean )
                       return string is
    begin
      if k > t.dg'last then
        return accu;
      elsif t.dg(k) = 0 then
        return Rec_Write(k+1,accu,first);
      else
        declare
          sb : constant Symbol := s(k); -- Symbol_Table.Get(natural32(k));
        begin
          if t.dg(k) > 1 then
            if first then
              declare
                new_accu : constant string
                         := Concat_Symbol0(accu,sb)
                            & '^' & Convert(integer32(t.dg(k)));
              begin
                return Rec_Write(k+1,new_accu,false);
              end;
            else
              declare
                new_accu : constant string
                         := Concat_Symbol1(accu,sb)
                             & '^' & Convert(integer32(t.dg(k)));
              begin
                return Rec_Write(k+1,new_accu,first);
              end;
            end if;
          elsif first then
            declare
              new_accu : constant string := Concat_Symbol0(accu,sb);
            begin
              return Rec_Write(k+1,new_accu,false);
            end;
          else
            declare
              new_accu : constant string := Concat_Symbol1(accu,sb);
            begin
              return Rec_Write(k+1,new_accu,first);
            end;
          end if;
        end;
      end if;
    end Rec_Write;
 
  begin
    declare
      factor : constant string := Rec_Write(1,"",true);
    begin
      if factor = "" then   -- then we have constant term
        if lead
         then return Unsigned_Constant(t.cf);
         else return Signed_Constant(t.cf);
        end if;
      elsif lead then
        if Is_Unit(t.cf) then
          if Sign(t.cf) = +1
           then return factor;
           else return " - " & factor;
          end if;
        else
          return Unsigned_Coefficient(t.cf) & "*" & factor;
        end if;
      elsif Is_Unit(t.cf) then
        if Sign(t.cf) = +1
         then return " + " & factor;
         else return " - " & factor;
        end if;
      else
        return Signed_Coefficient(t.cf) & "*" & factor;
      end if;
    end;
  end Write;

  function Write ( t : Term; lead : boolean ) return string is

  -- DESCRIPTION :
  --   Returns the string representation of the term.
  --   When the term is the leading term of the polynomial,
  --   then lead is true.

    s : Array_of_Symbols(t.dg'range);

  begin
    for k in s'range loop
      s(k) := Symbol_Table.Get(natural32(k));
    end loop;
    return Write(t,lead,s);
  end Write;

-- TARGET FUNCTIONS for PARSE OPERATORS :

  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p : in out Poly ) is

    bc : integer32 := 0;

  begin
    Parse_Polynomial(s,bc,k,n,p);
  end Parse;

  procedure Parse ( s : in string; k : in out integer;
                    n : in natural32; p,p_last : in out Term_List ) is

    bc : integer32 := 0;

  begin
   -- put_line("In procedure Parse_Polynomial ...");
    Parse_Polynomial(s,bc,k,n,p,p_last);
  end Parse;

  function Parse ( n : natural32; s : string ) return Poly is

    res : Poly;
    p : integer := s'first;

  begin
    Parse(s,p,n,res);
    return res;
  end Parse;

  function Parse ( n : natural32; s : string ) return Term_List is

    res,res_last : Term_List;
    p : integer := s'first;

  begin
   -- put_line("Calling procedure Parse ...");
    Parse(s,p,n,res,res_last);
    return res;
  exception
    when others =>
      put_line("Exception in Parse of string for term list"); raise;
  end Parse;

  function Parse ( n,m : natural32; s : string ) return Poly_Sys is

    res : Poly_Sys(1..integer32(n));
    ind : constant Standard_Natural_Vectors.Vector(1..integer32(n))
        := Delimiters(n,s);

  begin
   -- This function may be called when the symbol table has not yet
   -- been initialized and then it should not crash.
    if Symbol_Table.Number < m then
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    res(1) := Parse(m,s(s'first..integer(ind(1))));
    for i in 2..integer32(n) loop
      res(i) := Parse(m,s(integer(ind(i-1)+1)..integer(ind(i))));
    end loop;
    return res;
  end Parse;

  function Parse ( n,m : natural32; s : string ) return Array_of_Term_Lists is

    res : Array_of_Term_Lists(1..integer32(n));
    ind : constant Standard_Natural_Vectors.Vector(1..integer32(n))
        := Delimiters(n,s);

  begin
    if Symbol_Table.Number < m then -- same comment as other Parse
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    res(1) := Parse(m,s(s'first..integer(ind(1))));
    for i in 2..integer32(n) loop
      res(i) := Parse(m,s(integer(ind(i-1)+1)..integer(ind(i))));
    end loop;
    return res;
  end Parse;

  function Parse ( m : natural32; s : Array_of_Strings ) return Poly_Sys is

    res : Poly_Sys(integer32(s'first)..integer32(s'last));
 
  begin
    if Symbol_Table.Number < m then -- same comment as other Parse
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    for i in s'range loop
      declare
      begin
        res(integer32(i)) := Parse(m,s(i).all);
      exception
        when others =>
          put("Exception raised at string "); put(integer32(i),1);
          put_line(" in Parse to polynomial system."); raise;
      end;
    end loop;
    return res;
  end Parse;

  function Parse ( m : natural32; s : Array_of_Strings )
                 return Array_of_Term_Lists is

    res : Array_of_Term_Lists(integer32(s'first)..integer32(s'last));
 
  begin
    if Symbol_Table.Number < m then -- same comment as other Parse
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    for i in s'range loop
      declare
      begin
        res(integer32(i)) := Parse(m,s(i).all);
      exception
        when others =>
          put("Exception raised at string "); put(integer32(i),1);
          put_line(" in Parse to polynomial system."); raise;
      end;
    end loop;
    return res;
  end Parse;

  function Size_Limit ( p : Poly ) return natural32 is

    nbtrm : constant natural64 := natural64(Number_of_Terms(p));
    nbvar : constant natural64 := natural64(Number_of_Unknowns(p));
    symsz : constant natural64 := 5;
    cffsz : constant natural64 := 40;
    bound : constant natural64 := 2**31 - 1;
    res : constant natural64 := nbtrm*nbvar*symsz*cffsz;

  begin
    if res > bound
     then return natural32(bound);
     else return natural32(res);
    end if;
  end Size_Limit;

  function Write ( p : Poly; s : Array_of_Symbols ) return string is

    nb : constant natural32 := Number_of_Terms(p);
    terms : array(1..nb) of Term;
    ind : natural32 := 0;

    procedure Store_Term ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      terms(ind) := t;
      continue := true;
    end Store_Term;
    procedure Store_Terms is new Visiting_Iterator(Store_Term);
 
    function Write_Terms ( tfirst,tlast : natural32 ) return string is

    -- DESCRIPTION :
    --   Returns the string representation of all terms in the
    --   range tfirst..tlast of terms.

    begin
      if tfirst > tlast then
        return "";
      elsif tfirst = tlast then
        if tfirst = 1
         then return Write(terms(tfirst),true,s);  -- lead term
         else return Write(terms(tfirst),false,s); -- trailing term
        end if;
      else
        declare
          middle : constant natural32 := (tfirst+tlast)/2;
          midterm : constant string := Write(terms(middle),false,s);
          firsthalf : constant string := Write_Terms(tfirst,middle-1);
          secondhalf : constant string := Write_Terms(middle+1,tlast);
        begin
          return firsthalf & midterm & secondhalf;
        end;
      end if;
    end Write_Terms;

  begin
    Store_Terms(p);
    return Write_Terms(1,nb) & ";";
  end Write;

  function Write ( p : Poly ) return string is

    n : constant natural32 := Number_of_Unknowns(p);
    s : Array_of_Symbols(1..integer32(n));

  begin
    for k in s'range loop
      s(k) := Symbol_Table.Get(natural32(k));
    end loop;
    return Write(p,s);
  end Write;

  function Write ( p : Poly_Sys ) return string is
  begin
    if p'first = p'last
     then return Write(p(p'first));
     else return Write(p(p'first)) & Write(p(p'first+1..p'last));
    end if;
  end Write;

  function Write ( p : Poly_Sys; s : Array_of_Symbols ) return string is
  begin
    if p'first = p'last
     then return Write(p(p'first),s);
     else return Write(p(p'first),s) & Write(p(p'first+1..p'last),s);
    end if;
  end Write;

  function Write ( p : Poly_Sys ) return Array_of_Strings is

    res : Array_of_Strings(integer(p'first)..integer(p'last));

  begin
    for i in res'range loop
      declare
        s : constant string := Write(p(integer32(i)));
      begin
        res(i) := new string'(s);
      end;
    end loop;
    return res;
  end Write;

  function Write ( p : Poly_Sys; s : Array_of_Symbols )
                 return Array_of_Strings is

    res : Array_of_Strings(integer(p'first)..integer(p'last));

  begin
    for i in res'range loop
      declare
        pstr : constant string := Write(p(integer32(i)),s);
      begin
        res(i) := new string'(pstr);
      end;
    end loop;
    return res;
  end Write;

end Standard_Complex_Poly_Strings;
