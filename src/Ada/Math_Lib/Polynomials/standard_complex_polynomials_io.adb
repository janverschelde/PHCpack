with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Parse_Numbers;             use Standard_Parse_Numbers;
with Standard_Write_Numbers;             use Standard_Write_Numbers;
with Standard_Natural_Vectors;
with Symbol_Table_io;
with Line_Breaks;                        use Line_Breaks;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions;
with Write_Factors;                      use Write_Factors;

package body Standard_Complex_Polynomials_io is

-- AUXILIARIES FOR THE INPUT ROUTINES :

  procedure Read_Term ( file : in file_type; bc : in out integer32;
                        char : in out character;
                        n : in natural32; termp : in out Poly );
  -- DESCRIPTION :
  --   Reads a term from file, char is the first character of the term.

  -- ON ENTRY :
  --   file     must be opened for input;
  --   bc       counts #open - #closed brackets;
  --   char     first character of the term;
  --   n        number of variables;
  --   termp    accumulating polynomial for the term.

  -- ON RETURN :
  --   bc       updated bracket counter;
  --   char     first character following the term;
  --   termp    resulting term read from file.

  procedure Read_Factor ( file : in file_type; bc : in out integer32;
                          char : in out character; n : in natural32;
                          d : in out Degrees; pb : in out Poly );
  -- DESCRIPTION :
  --   Reads a factor from file, char is the first character of the factor.

  -- ON ENTRY :
  --   file     must be opened for input;
  --   bc       counts #open - #closed brackets;
  --   char     first character of the factor;
  --   n        number of variables;
  --   d        accumulates the degrees of the factor;
  --   pb       accumulating polynomial for the factor.

  -- ON RETURN :
  --   bc       updated bracket counter;
  --   char     first character following the factor;
  --   pb       resulting factor read from file.

  procedure Read_Power_Factor
              ( file : in file_type; char : in out character;
                p : in out Poly ) is

  -- DESCRIPTION :
  --   Reads the exponent of a factor stored in p and returns
  --   in p the polynomial p raised to the read exponent.

  -- ON ENTRY :
  --   file     must be opened for input;
  --   char     equals the exponentiation symbol '^';
  --   p        current factor.

  -- ON RETURN :
  --   char     character following the exponent;
  --   p        p raised to the read exponent.
 
    d,ne,ne2 : natural32 := 0;
    expo,expo2 : integer32 := 1;
    sign : character;
    bracket : boolean := false;

  begin
    get(file,char);                           -- skip the '^'
    if char = '('
     then get(file,char); bracket := true;
    end if;
    Parse(file,char,expo,expo2,ne,ne2,sign);
    if bracket then
      Skip_Spaces_and_CR(file,char);
      if char = ')'
       then get(file,char);     -- skip closing bracket
       else raise BAD_BRACKET;  -- no closing bracket found
      end if;
    end if;
    d := natural32(expo);
    Pow(p,d);
  end Read_Power_Factor;

  procedure Read_Polynomial
              ( file : in file_type; bc : in out integer32;
                p : out Poly ) is

  -- DESCRIPTION :
  --   Reads a polynomial from file, returned in p, raising BAD_BRACKET
  --   if the bracket counter bc is nonzero when the delimiter is encountered.

    n : constant natural32 := Symbol_Table.Maximal_Size;
    char,oper,nextchar : character;
    term,res,acc : Poly := Null_Poly;
    eol : boolean;

  begin
    oper := '+';
    get(file,char);
    Skip_Spaces_and_CR(file,char);
    if char = '-'
     then oper := '-';
    end if;                         -- the first term can have no sign
    Read_Term(file,bc,char,n,res);  -- therefore read it first
    loop
      case char is
        when ' ' | ASCII.CR | ASCII.LF => get(file,char);  -- skip blanks
        when '+' | '-' =>
          oper := char;
          Read_Term(file,bc,char,n,term);
          Add(res,term); Clear(term);
        when delimiter =>
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' =>
          bc := bc + 1;
          Read_Polynomial(file,bc,term);
          get(file,char);
          Skip_Spaces_and_CR(file,char);
          if char = '(' then
            raise BAD_BRACKET;
         -- elsif char = ')' then
         --   bc := bc - 1;
          end if;
          if char = '^' then
            Read_Power_Factor(file,char,term);
          elsif char = '*' then  -- look for the **2
            look_ahead(file,nextchar,eol);
            if nextchar = '*' then
              get(file,char);
              Read_Power_Factor(file,char,term);
            end if;
          end if;
          case oper is
            when '+' => Add(acc,res); Clear(res); Copy(term,res);
            when '-' => Add(acc,res);Clear(res);
                        Copy(term,res); Min(res);
            when '*' => Mul(res,term);
            when others => raise ILLEGAL_OPERATION;
          end case;
          Clear(term);
        when ')' =>
          bc := bc - 1;
          if bc < 0 then
            raise BAD_BRACKET;
          end if;
          exit;
        when '*' =>
         -- if res = Null_Poly then    -- input case (0.0)*(2+x)
         --   raise ILLEGAL_CHARACTER; -- should not be an error
         -- else -- the case " ) * " :
            oper := char; get(file,char);  -- skip '*'
            Read_Term(file,bc,char,n,term);
            Skip_Spaces_and_CR(file,char);
            if char = '^' then
              Read_Power_Factor(file,char,term);
            elsif char = '*' then
              look_ahead(file,nextchar,eol);
              if nextchar = '*' then
                get(file,char);
                Read_Power_Factor(file,char,term);
              end if;
            end if;
            if char /= '(' then
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
           else Read_Power_Factor(file,char,res);
          end if;
        when others => raise ILLEGAL_CHARACTER;
      end case;
    end loop;
    p := acc + res;
    Clear(acc); Clear(res);
  end Read_Polynomial;

  procedure Read_Factor ( file : in file_type; bc : in out integer32;
                          char : in out character; n : in natural32;
                          d : in out Degrees; pb : in out Poly ) is

    sb : symbol;
    ne,ne2 : natural32 := 0;
    k : integer32;
    expo,expo2 : integer32 := 1;
    sign : character;
 
  begin
    Skip_Spaces_and_CR(file,char);
    if char = '(' then 
      bc := bc + 1;
      Read_Polynomial(file,bc,pb);
      get(file,char);
      Skip_Spaces_and_CR(file,char);
      if char = '^'
       then Read_Power_Factor(file,char,pb);
      end if;
      return;
    end if;
    Symbol_Table_io.Parse_Symbol(file,char,sb);
    k := integer32(Symbol_Table.Check_Symbol(n,sb));
    Skip_Spaces_and_CR(file,char);
    if char = '^' then
      get(file,char);                                    -- skip the '^'
      Parse(file,char,expo,expo2,ne,ne2,sign);
      d(k) := d(k) + natural32(expo);
      Skip_Spaces_and_CR(file,char);
      if char /= '*'                            -- the case x^2*...
       then return;                             -- end of factor
       else get(file,char);                     -- skip the '*'
      end if; 
    elsif char = '*' then
      get(file,char);
      if char = '*' then
        get(file,char);                 -- the case " x ** expo "
        Parse(file,char,expo,expo2,ne,ne2,sign);
        d(k) := d(k) + natural32(expo);
        Skip_Spaces_and_CR(file,char);
        if char /= '*'
         then return;                   -- end of factor
         else get(file,char);           -- skip the '*'
        end if;
      else
        d(k) := d(k) + 1;               -- the case " x * ? "
      end if;
    else -- the case " x ?", with ? /= '*' or ' ' or '^' :
      d(k) := d(k) + 1;
      return;
    end if;
    Skip_Spaces_and_CR(file,char);
    if (char = '-') or (char = '+') 
     then return;
    end if;
    if Convert(char) < 10
     then Read_Term(file,bc,char,n,pb);     -- case " x * c " or " x ** c * c " 
     else Read_Factor(file,bc,char,n,d,pb); -- case " x * y "
    end if;
  end Read_Factor;
 
  procedure Read_Term ( file : in file_type; bc : in out integer32;
                        char : in out character;
                        n : in natural32; termp : in out Poly ) is

    c : Complex_Number;
    d : Degrees := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    pb,res : Poly := Null_Poly;
    tmp : Term;

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
    Parse(file,char,c); -- look for 'i' :
    Skip_Spaces_and_CR(file,char);
    if ( c = Create(0.0) ) and then ((char = 'i') or (char = 'I')) then
      -- the case "+ i" :
      c := Create(0.0,1.0); 
      get(file,char);        -- skip 'i'
    elsif ( c = Create(-1.0) ) and then ((char = 'i') or (char = 'I')) then
      -- the case "- i" :
      c := Create(0.0,-1.0);
      get(file,char);    -- skip 'i'
    elsif char = '*' then -- the case ".. c *.." :
      Skip_Spaces_and_CR(file,char);
      get(file,char);  -- skip '*'
      Skip_Spaces_and_CR(file,char);
      if (char = 'i') or (char = 'I') then -- the case ".. c * i.." :
        c := Create(0.0,REAL_PART(c));
        get(file,char);    -- skip 'i'
      else                                 -- the case ".. c * x.." :
        Read_Factor(file,bc,char,n,d,pb);
        if pb /= Null_Poly
         then Clear(res); Copy(pb,res); Clear(pb);
        end if;
      end if;
    end if; -- the case ".. c ?"  will be treated in the loop
    loop
      case char is
        when ' ' | ASCII.CR | ASCII.LF => get(file,char);
        when '*' =>
          get(file,char); Read_Factor(file,bc,char,n,d,pb);
          Collect_Factor_Polynomial;
        when '+' | '-' => 
          -- if c = Create(0.0)
          --  then raise ILLEGAL_CHARACTER;
          --  else exit;
          -- end if;
          exit; -- zero coefficient no longer wrong input
        when delimiter => 
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' =>  -- no modification of bc by Read_Term
          if c = Create(0.0) or else c = Create(-1.0)
           then c := Create(0.0); exit; -- the case "+ (" or "- (" :
           else raise BAD_BRACKET;      -- the case "c  (" :
          end if;
        when ')' =>  -- no modification of bc by Read_Term
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when others =>
          if c = Create(0.0) then
            c := Create(1.0);
            Read_Factor(file,bc,char,n,d,pb);
          elsif c = Create(-1.0) then
            Read_Factor(file,bc,char,n,d,pb);
          elsif char = '^' then
            Read_Power_Factor(file,char,res);
          else
            raise ILLEGAL_CHARACTER;
          end if;
          Collect_Factor_Polynomial;
      end case;
    end loop;
    if c /= Create(0.0) then
      tmp.cf := c;
      tmp.dg := d;
      termp := create(tmp);
      Clear(tmp);
      if Number_Of_Unknowns(res) > 0
       then Mul(termp,res); Clear(res);
      end if;
    end if;
  end Read_Term;

-- THE INPUT OPERATIONS :

  procedure get ( p : in out Poly ) is
  begin
    get(standard_input,p);
  end get;

  procedure get ( file : in file_type; p : in out Poly ) is

    bc : integer32 := 0;

  begin
    Read_Polynomial(file,bc,p);
  end get;
 
  procedure get ( n : in out natural32; p : in out Poly ) is
  begin
    get(standard_input,n,p);
  end get;

  procedure get ( file : in file_type;
                  n : in out natural32; p : in out Poly ) is
  begin
    get(file,n);
    if Symbol_Table.Empty
     then Symbol_Table.Init(n);
    end if;
    get(file,p);
  end get;

-- THE OUTPUT OPERATIONS :

  procedure put_terms ( p : in Poly; pow : in power ) is
  begin
    put_terms(Standard_Output,p,pow);
  end put_terms;

  procedure put ( p : in Poly; pow : in power ) is
  begin
    put(Standard_Output,p,pow);
  end put;

  procedure put_terms ( file : in file_type; p : in Poly; pow : in power ) is

    nn : constant natural32 := Number_of_Unknowns(p);
    standard : constant boolean := ( Symbol_Table.Number < nn );
    first_time : boolean := true;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
 
    -- DESCRIPTION : 
    --   Writes a term is written on file.
 
      passed : boolean;
      sc : natural32 := 0;
      nc : natural32;
      degtg : natural32 := 0;

    begin
      if first_time 
       then first_time := false;
       else Write_Plus(file,t.cf,nc); sc := sc + nc;
      end if;
      for i in t.dg'range loop
        degtg := degtg + t.dg(i);
      end loop;
      if degtg = 0 then -- if Sum(t.dg) = 0 then -- crash on x + (-1-2*i)*y;
        Write_Number(file,t.cf,nc); sc := sc + nc;
      else
        Write_Coefficient(file,t.cf,nc); sc := sc + nc;
        passed := false;
        for i in t.dg'range loop
          if t.dg(i) > 0 then
            if passed
             then put(file,'*'); sc := sc + 1;
             else passed := true;
            end if;
            sc := sc + Length(integer32(t.dg(i)),i,standard,pow);
            Write_Factor(file,t.dg(i),natural32(i),standard,pow);
          end if;
        end loop;
        Line(file,sc);
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Init_Line;
    Write_Terms(p);
  end put_terms;

  procedure put ( file : in file_type; p : in Poly; pow : in power ) is
  begin
    put_terms(file,p,pow);
    put(file,delimiter); Line(file,1);
  end put;

  procedure put_terms ( n : in natural32; p : in Poly; pow : in power ) is
  begin
    put_terms(Standard_Output,n,p,pow);
  end put_terms;

  procedure put ( n : in natural32; p : in Poly; pow : in power ) is
  begin
    put(Standard_Output,n,p,pow);
  end put;

  procedure put_terms ( file : in file_type; n : in natural32;
                        p : in Poly; pow : in power ) is
  begin
    put(file,n,1);
    put_line(file," ");
    put_terms(file,p,pow);
  end put_terms;

  procedure put ( file : in file_type; n : in natural32;
                  p : in Poly; pow : in power ) is
  begin
    put(file,n,1);
    put_line(file," ");
    put(file,p,pow);
  end put;

  procedure put_terms ( p : in Poly ) is
  begin
    put_terms(Standard_Output,p,'^');
  end put_terms;

  procedure put ( p : in Poly ) is
  begin
    put(Standard_Output,p,'^');
  end put;
 
  procedure put_terms ( file : in file_type; p : in Poly ) is
  begin
    put_terms(file,p,'^');
  end put_terms;
 
  procedure put ( file : in file_type; p : in Poly ) is
  begin
    put(file,p,'^');
  end put;

  procedure put_terms ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_terms(Standard_Output,p,s);
  end put_terms;

  procedure put ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put(Standard_Output,p,s);
  end put;

  procedure put_terms ( file : in file_type;
                        p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_terms(file,p,s,'^');
  end put_terms;

  procedure put ( file : in file_type;
                  p : in Poly; s : in Array_of_Symbols ) is
  begin
    put(file,p,s,'^');
  end put;

  procedure put_terms ( p : in Poly; s : in Array_of_Symbols;
                        pow : in Power ) is
  begin
    put_terms(Standard_Output,p,s,pow);
  end put_terms;

  procedure put ( p : in Poly; s : in Array_of_Symbols; pow : in Power ) is
  begin
    put(Standard_Output,p,s,pow);
  end put;

  procedure put_terms ( file : in file_type; p : in Poly;
                        s : in Array_of_Symbols; pow : in Power ) is

    first_time : boolean := true;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
 
    -- DESCRIPTION : 
    --   Writes a term is written on file.
 
      passed : boolean;
      sc : natural32 := 0;
      nc : natural32;

    begin
      if first_time 
       then first_time := false;
       else Write_Plus(file,t.cf,nc); sc := sc + nc;
      end if;
      if Sum(t.dg) = 0 then
        Write_Number(file,t.cf,nc); sc := sc + nc;
      else 
        Write_Coefficient(file,t.cf,nc); sc := sc + nc;
        passed := false;
        for i in t.dg'range loop
          if t.dg(i) > 0 then
            if passed
             then put(file,'*'); sc := sc + 1;
             else passed := true;
            end if;
            sc := sc + Length(integer32(t.dg(i)),i,true,pow);
            Write_Factor(file,t.dg(i),s(i),pow);
          end if;
        end loop;
        Line(file,sc);
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Init_Line;
    Write_Terms(p);
  end put_terms;

  procedure put ( file : in file_type;
                  p : in Poly; s : in Array_of_Symbols; pow : in Power ) is
  begin
    put_terms(file,p,s,pow);
    put(file,delimiter); Line(file,1);
  end put;

  procedure put_terms ( n : in natural32; p : in Poly;
                        s : in Array_of_Symbols; pow : in Power ) is
  begin
    put_terms(Standard_Output,n,p,s,pow);
  end put_terms;

  procedure put ( n : in natural32; p : in Poly;
                  s : in Array_of_Symbols; pow : in Power ) is
  begin
    put(Standard_Output,n,p,s,pow);
  end put;

  procedure put_terms ( file : in file_type;
                        n : in natural32; p : in Poly;
                        s : in Array_of_Symbols; pow : in Power ) is
  begin
    put(file,n,1);
    put_line(file," ");
    put_terms(file,p,s,pow);
  end put_terms;

  procedure put ( file : in file_type;
                  n : in natural32; p : in Poly;
                  s : in Array_of_Symbols; pow : in Power ) is
  begin
    put(file,n,1);
    put_line(file," ");
    put(file,p,s,pow);
  end put;

  procedure put_terms ( p : in Poly; dp : in natural32 ) is
  begin
    put_terms(Standard_Output,p,dp);
  end put_terms;

  procedure put ( p : in Poly; dp : in natural32 ) is
  begin
    put(Standard_Output,p,dp);
  end put;

  procedure put_terms ( file : in file_type;
                         p : in Poly; dp : in natural32 ) is
  begin
    put_terms(file,p);
  end put_terms;

  procedure put ( file : in file_type; p : in Poly; dp : in natural32 ) is
  begin
    put(file,p);
  end put;

  procedure put_terms_line ( p : in Poly; pow : in power ) is
  begin
    put_terms_line(Standard_Output,p,pow);
  end put_terms_line;

  procedure put_terms_pair ( p : in Poly; pow : in power ) is
  begin
    put_terms_pair(Standard_Output,p,pow);
  end put_terms_pair;

  procedure put_line ( p : in Poly; pow : in power ) is
  begin
    put_line(Standard_Output,p,pow);
  end put_line;

  procedure put_pair ( p : in Poly; pow : in power ) is
  begin
    put_pair(Standard_Output,p,pow);
  end put_pair;

  procedure put_terms_line ( p : in Poly ) is
  begin
    put_terms_line(Standard_Output,p,'^');
  end put_terms_line;

  procedure put_terms_pair ( p : in Poly ) is
  begin
    put_terms_pair(Standard_Output,p,'^');
  end put_terms_pair;

  procedure put_line ( p : in Poly ) is
  begin
    put_line(Standard_Output,p,'^');
  end put_line;

  procedure put_pair ( p : in Poly ) is
  begin
    put_pair(Standard_Output,p,'^');
  end put_pair;

  procedure put_terms_line ( file : in file_type; p : in Poly ) is
  begin
    put_terms_line(file,p,'^');
  end put_terms_line;

  procedure put_terms_pair ( file : in file_type; p : in Poly ) is
  begin
    put_terms_pair(file,p,'^');
  end put_terms_pair;

  procedure put_line ( file : in file_type; p : in Poly ) is
  begin
    put_line(file,p,'^');
  end put_line;

  procedure put_pair ( file : in file_type; p : in Poly ) is
  begin
    put_pair(file,p,'^');
  end put_pair;

  procedure put_terms_line ( file : in file_type;
                             p : in Poly; pow : in power ) is

    n : constant natural32 := Number_of_Unknowns(p);
    standard : constant boolean := ( Symbol_Table.Number < n );
    nc : natural32;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      new_line(file);
      if Is_Real(t.cf) then
        if REAL_PART(t.cf) >= 0.0
         then put(file,"+");
        end if;
      else
        if ((REAL_PART(t.cf) /= 0.0) or (IMAG_PART(t.cf) > 0.0))
         then put(file,"+");
        end if;
      end if;
      Init_Line; Write_Number(file,t.cf,nc); Line(file,nc);
      if Sum(t.dg) /= 0 then
        for i in t.dg'range loop
          if t.dg(i) > 0 then
            put(file,'*');
            Write_Factor(file,t.dg(i),natural32(i),standard,pow);
          end if;
        end loop;
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator (process => Write_Term);

  begin
    Write_Terms(p);
  end put_terms_line;

  procedure put_terms_pair ( file : in file_type;
                             p : in Poly; pow : in power ) is

    n : constant natural32 := Number_of_Unknowns(p);
    standard : constant boolean := ( Symbol_Table.Number < n );
    cnt : natural32 := 0;
    nc : natural32;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      if cnt mod 2 = 1
       then new_line(file);
      end if;
      if Is_Real(t.cf) then
        if REAL_PART(t.cf) >= 0.0
         then put(file,"+");
        end if;
      else
        if ((REAL_PART(t.cf) /= 0.0) or (IMAG_PART(t.cf) > 0.0))
         then put(file,"+");
        end if;
      end if;
      Init_Line; Write_Number(file,t.cf,nc); Line(file,nc);
      if Sum(t.dg) /= 0 then
        for i in t.dg'range loop
          if t.dg(i) > 0 then
            put(file,'*');
            Write_Factor(file,t.dg(i),natural32(i),standard,pow);
          end if;
        end loop;
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator (process => Write_Term);

  begin
    Write_Terms(p);
  end put_terms_pair;

  procedure put_line ( file : in file_type; p : in Poly; pow : in power ) is
  begin
    put_terms_line(file,p,pow);
    put(file,delimiter);
    new_line(file);
  end put_line;

  procedure put_pair ( file : in file_type; p : in Poly; pow : in power ) is
  begin
    put_terms_pair(file,p,pow);
    put(file,delimiter);
    new_line(file);
  end put_pair;

  procedure put_terms_line
              ( p : in Poly; s : in Array_of_Symbols; pow : in power ) is
  begin
    put_terms_line(Standard_Output,p,s,pow);
  end put_terms_line;

  procedure put_terms_pair
              ( p : in Poly; s : in Array_of_Symbols; pow : in power ) is
  begin
    put_terms_pair(Standard_Output,p,s,pow);
  end put_terms_pair;

  procedure put_line
              ( p : in Poly; s : in Array_of_Symbols; pow : in power ) is
  begin
    put_line(Standard_Output,p,s,pow);
  end put_line;

  procedure put_pair
              ( p : in Poly; s : in Array_of_Symbols; pow : in power ) is
  begin
    put_pair(Standard_Output,p,s,pow);
  end put_pair;

  procedure put_terms_line ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_terms_line(Standard_Output,p,s,'^');
  end put_terms_line;

  procedure put_terms_pair ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_terms_pair(Standard_Output,p,s,'^');
  end put_terms_pair;

  procedure put_line ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_line(Standard_Output,p,s,'^');
  end put_line;

  procedure put_pair ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_pair(Standard_Output,p,s,'^');
  end put_pair;

  procedure put_terms_line ( file : in file_type; p : in Poly;
                             s : in Array_of_Symbols ) is
  begin
    put_terms_line(file,p,s,'^');
  end put_terms_line;

  procedure put_terms_pair ( file : in file_type; p : in Poly;
                             s : in Array_of_Symbols ) is
  begin
    put_terms_pair(file,p,s,'^');
  end put_terms_pair;

  procedure put_line ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols ) is
  begin
    put_line(file,p,s,'^');
  end put_line;

  procedure put_pair ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols ) is
  begin
    put_pair(file,p,s,'^');
  end put_pair;

  procedure put_terms_line
                ( file : in file_type; p : in Poly;
                  s : in Array_of_Symbols; pow : in power ) is

    procedure Write_Term ( t : in Term; continue : out boolean ) is

      nc : natural32;

    begin
      new_line(file);
      if Is_Real(t.cf) then
        if REAL_PART(t.cf) >= 0.0
         then put(file,"+");
        end if;
      else
        if ((REAL_PART(t.cf) /= 0.0) or (IMAG_PART(t.cf) > 0.0))
         then put(file,"+");
        end if;
      end if;
      Init_Line; Write_Number(file,t.cf,nc); Line(file,nc);
      if Sum(t.dg) /= 0 then
        for i in t.dg'range loop
          if t.dg(i) > 0 then
            put(file,'*');
            Write_Factor(file,t.dg(i),s(i),pow);
          end if;
        end loop;
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
  end put_terms_line;

  procedure put_terms_pair
                ( file : in file_type; p : in Poly;
                  s : in Array_of_Symbols; pow : in power ) is

    cnt : natural32 := 0;
    nc : natural32;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      if cnt mod 2 = 1
       then new_line(file);
      end if;
      if Is_Real(t.cf) then
        if REAL_PART(t.cf) >= 0.0
         then put(file,"+");
        end if;
      else
        if ((REAL_PART(t.cf) /= 0.0) or (IMAG_PART(t.cf) > 0.0))
         then put(file,"+");
        end if;
      end if;
      Init_Line; Write_Number(file,t.cf,nc); Line(file,nc);
      if Sum(t.dg) /= 0 then
        for i in t.dg'range loop
          if t.dg(i) > 0 then
            put(file,'*');
            Write_Factor(file,t.dg(i),s(i),pow);
          end if;
        end loop;
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
  end put_terms_pair;

  procedure put_line ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols; pow : in power ) is
  begin
    put_terms_line(file,p,s,pow);
    put(file,delimiter);
    new_line(file);
  end put_line;

  procedure put_pair ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols; pow : in power ) is
  begin
    put_terms_pair(file,p,s,pow);
    put(file,delimiter);
    new_line(file);
  end put_pair;

  procedure Display_Format is

    s : array(1..24) of string(1..65);

  begin
    s( 1):="  A complex multivariate polynomial is denoted as a  sequence  of";
    s( 2):="terms, separated by `+' and terminated by the semicolon `;'.  The";
    s( 3):="brackets '(' and ')' must be used to isolate a sequence of  terms";
    s( 4):="as a factor in a complex multivariate polynomial.                ";
    s( 5):="  A term can be either a coefficient or a  coefficient,  followed";
    s( 6):="by  '*'  and  a  monomial.  If in the latter case the coefficient";
    s( 7):="equals one, then it may be omitted.                              ";
    s( 8):="  A coefficient may be denoted  as  an  integer,  a  rational,  a";
    s( 9):="floating-point or a complex number.                              ";
    s(10):="  A monomial is a sequence of powers of  unknowns,  separated  by";
    s(11):="'*'.   The power operator is represented by '**' or '^'.  It must";
    s(12):="be followed by a positive natural number.  If  the  power  equals";
    s(13):="one, then it may be omitted.                                     ";
    s(14):="  An unknown can be denoted by at most 80 characters.  The  first";
    s(15):="character  must  be a letter and the other two characters must be";
    s(16):="different from '+', '-', '*', '^', '/', ';', '('  and  ')'.   The";
    s(17):="letter i means sqrt(-1), whence it does not represent an unknown.";
    s(18):="The number of unknowns may not  exceed  the  declared  dimension.";
    s(19):="  Some  examples  of  valid  notations  of  complex  multivariate";
    s(20):="polynomials:                                                     ";
    s(21):="  x**2*y + 1/2*z*y**2 - 2*z + y**3 + x - 1E9/-8.E-6* y + 3;      ";
    s(22):="  x^2*y + z*y^2 - 2*z + y^3 + x - y + 3;                         ";
    s(23):="  (1.01 + 2.8*i)*x1**2*x2 + x3**2*x1 - 3*x1 + 2*x2*x3 - 3;       ";
    s(24):="  (x1^2*x2 + x3^2*x1 - 3*x1 + 2*x2*x3 - 3)*x2**2*(x2-1+i);       ";
    for i in s'range loop
      put_line(s(i));
    end loop;
  end Display_Format;

end Standard_Complex_Polynomials_io;
