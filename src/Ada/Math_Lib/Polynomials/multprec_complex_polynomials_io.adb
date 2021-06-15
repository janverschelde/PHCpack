with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Parse_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Multprec_Parse_Numbers;             use Multprec_Parse_Numbers;
with Multprec_Write_Numbers;             use Multprec_Write_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Symbol_Table_io;
with Write_Factors;                      use Write_Factors;
with Standard_Complex_Polynomials_io;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions;

package body Multprec_Complex_Polynomials_io is

  working_precision : natural32 := 4;

-- SETTING THE WORKING PRECISION :

  procedure Set_Working_Precision ( size : in natural32 ) is
  begin
    working_precision := size;
  end Set_Working_Precision;

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
 
    d,ne : natural32 := 0;
    expo : integer32 := 1;
    sign : character;
    bracket : boolean := false;

  begin
    get(file,char);                           -- skip the '^'
    if char = '('
     then get(file,char); bracket := true;
    end if;
    Standard_Parse_Numbers.Parse(file,char,expo,ne,sign);
    if bracket then
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      if char = ')'
       then get(file,char);      -- skip closing bracket
       else raise BAD_BRACKET;   -- no closing bracket found
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
    term,res,acc : Poly;
    eol : boolean;

  begin
    oper := '+';
    get(file,char);
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if char = '-'
     then oper := '-';
    end if;                         -- the first term can have no sign
    Read_Term(file,bc,char,n,res);  -- therefore read it first
    loop
      case char is
        when ' ' | ASCII.CR | ASCII.LF => get(file,char);    -- skip blanks
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
          Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
          if char = '(' -- or char = ')'
           then raise BAD_BRACKET;
          end if;
          if char = '^' then
            Read_Power_Factor(file,char,term);
          elsif char = '*' then
            look_ahead(file,nextchar,eol);
            if nextchar = '*' then
              get(file,char);
              Read_Power_Factor(file,char,term);
            end if;
          end if;
          case oper is
            when '+' => Add(acc,res); Clear(res);
                        Copy(term,res);
            when '-' => Add(acc,res);Clear(res);
                        Copy(term,res); Min(res);
            when '*' => Mul(res,term);
            when others => raise ILLEGAL_OPERATION;
          end case;
          Clear(term);
        when ')' =>
          bc := bc - 1;
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '*' =>
         -- if res = Null_Poly then     -- input case like (0.0)*(2 + x)
         --   raise ILLEGAL_CHARACTER;  -- should not be an error
         -- else -- the case " ) * " :
            oper := char; get(file,char);  -- skip '*'
            Read_Term(file,bc,char,n,term);
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
    k : integer32 := 0;
    ne : natural32 := 0;
    expo : integer32 := 1;
    sign,nextchar : character;
    eol : boolean;
 
  begin
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if char = '(' then
      bc := bc + 1;
      Read_Polynomial(file,bc,pb);
      get(file,char);
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      if char = '^' then
        Read_Power_Factor(file,char,pb);
      elsif char = '*' then
        look_ahead(file,nextchar,eol);
        if nextchar = '*' then
          get(file,char);
          Read_Power_Factor(file,char,pb);
        end if;
      end if;
      return;
    end if;
    Symbol_Table_io.Parse_Symbol(file,char,sb);
    k := integer32(Symbol_Table.Check_Symbol(n,sb));
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if char = '^' then
      get(file,char);                                    -- skip the '^'
      Standard_Parse_Numbers.Parse(file,char,expo,ne,sign);
      d(k) := d(k) + natural32(expo);
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      if char /= '*'                            -- the case x^2*...
       then return;                             -- end of factor
       else get(file,char);                     -- skip the '*'
      end if; 
    elsif char = '*' then
      get(file,char);
      if char = '*' then
        get(file,char);                 -- the case " x ** expo "
        Standard_Parse_Numbers.Parse(file,char,expo,ne,sign);
        d(k) := d(k) + natural32(expo);
        Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
        if char /= '*'
         then return;                   -- end of factor
         else get(file,char);           -- skip the '*'
        end if;
       else
         d(k) := d(k) + 1;              -- the case " x * ? "
      end if;
    else -- the case " x ?", with ? /= '*' or ' ' or '^' :
      d(k) := d(k) + 1;
      return;
    end if;
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if (char = '-') or (char = '+') 
     then return;
    end if;
    if Convert(char) < 10
     then Read_Term(file,bc,char,n,pb); -- the case " x * c " or " x ** c * c "
     else Read_Factor(file,bc,char,n,d,pb); -- the case " x * y " 
    end if;
  end Read_Factor;

  procedure Read_Term ( file : in file_type; bc : in out integer32;
                        char : in out character;
                        n : in natural32; termp : in out Poly ) is

    c : Complex_Number;
    d : Degrees := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    pb,res : Poly;
    tmp : Term;
    realzero : constant Floating_Number := Create(integer(0));
    realmin1 : constant Floating_Number := Create(integer(-1));
    compzero : constant Complex_Number := Create(realzero);
    compmin1 : constant Complex_Number := Create(realmin1);
    nextchar : character;
    eol : boolean;

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
    Parse(file,working_precision,char,c); -- look for 'i' :
    Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
    if ( Equal(c,compzero) ) and then ((char = 'i') or (char = 'I')) then
      -- the case "+ i"
      c := Create(0.0,1.0); 
      get(file,char);        -- skip 'i'
    elsif ( Equal(c,compmin1) ) and then ((char = 'i') or (char = 'I')) then
      -- the case "- i"
      c := Create(0.0,-1.0);
      get(file,char);        -- skip 'i'
    elsif char = '*' then    -- the case ".. c *.." :
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      get(file,char);  -- skip '*'
      Standard_Parse_Numbers.Skip_Spaces_and_CR(file,char);
      if (char = 'i') or (char = 'I') then -- the case ".. c * i.."
        c := Create(realzero,REAL_PART(c));
        get(file,char);    -- skip 'i'
      else -- the case ".. c * x.." :
        Read_Factor(file,bc,char,n,d,pb);
        if pb /= Null_Poly
         then Clear(res); Copy(pb,res); Clear(pb);
        end if;
      end if;
    else -- the case ".. c ?" will be treated in the loop
      null;
    end if;
    loop
      case char is
        when ' ' | ASCII.CR =>
          get(file,char);
        when '*' => 
          get(file,char); Read_Factor(file,bc,char,n,d,pb);
          Collect_Factor_Polynomial;
        when '+' | '-' =>
          -- if Equal(c,compzero)
          --  then raise ILLEGAL_CHARACTER;
          --  else exit;
          -- end if;
          exit; -- zero coefficient should no be wrong input
        when delimiter =>
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' => 
          if Equal(c,compzero) or else Equal(c,compmin1)
           then c := Create(0.0); exit; -- the case "+ (" or "- (" 
           else raise BAD_BRACKET; -- the case "c  (" 
          end if;
        when ')' =>
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when others =>
          if Equal(c,compzero) then
            c := Create(1.0);
            Read_Factor(file,bc,char,n,d,pb);
          elsif Equal(c,compmin1) then
            Read_Factor(file,bc,char,n,d,pb);
          elsif char = '^' then
            Read_Power_Factor(file,char,res);
          elsif char = '*' then
            look_ahead(file,nextchar,eol);
            if nextchar = '*' then
              get(file,char);
              Read_Power_Factor(file,char,res);
            end if;
          else
            raise ILLEGAL_CHARACTER;
          end if;
          Collect_Factor_Polynomial;
      end case;
    end loop;
    if not Equal(c,compzero) then
      tmp.cf := c;
      tmp.dg := d;
      termp := create(tmp);
      if Number_Of_Unknowns(res) > 0
       then Mul(termp,res); Clear(res);
      end if;
    end if;
  end Read_Term;

-- THE INPUT OPERATIONS :

  procedure get ( n : in out natural32; p : in out Poly ) is
  begin
    get(Standard_Input,n,p);
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

  procedure get ( p : in out Poly ) is
  begin
    get(Standard_Input,p);
  end get;

  procedure get ( file : in file_type; p : in out Poly ) is

    bc : integer32 := 0;

  begin
    Read_Polynomial(file,bc,p);
  end get;

-- THE OUTPUT OPERATIONS FOR TERMS :

  procedure Write ( file : in file_type; d : in Degrees;
                    standard : in boolean; pow : in Power ) is
  begin
    for i in d'range loop
      if d(i) /= 0 then
        put(file,"*");
        Write_Factor(file,d(i),natural32(i),standard,pow);
      end if;
    end loop;
  end Write;

  procedure Write ( file : in file_type; d : in Degrees;
                    s : in Array_of_Symbols; pow : in Power ) is
  begin
    for i in d'range loop
      if d(i) /= 0 then
        put(file,"*");
        Write_Factor(file,d(i),s(i),pow);
      end if;
    end loop;
  end Write;

  procedure put ( file : in file_type; t : in Term;
                  standard : in boolean; pow : in Power ) is
    
    sumtdg : natural32 := 0;

  begin
    Write_Number(file,t.cf);
    for k in t.dg'range loop
      sumtdg := sumtdg + t.dg(k);
    end loop;
    if sumtdg /= 0 -- if Sum(t.dg) /= 0
     then Write(file,t.dg,standard,pow);
    end if;
  end put;

  procedure put ( file : in file_type; t : in Term;
                  s : in Array_of_Symbols; pow : in Power ) is

    sumtdg : natural32 := 0;

  begin
    Write_Number(file,t.cf);
    for k in t.dg'range loop
      sumtdg := sumtdg + t.dg(k);
    end loop;
    if sumtdg /= 0 -- if Sum(t.dg) /= 0
     then Write(file,t.dg,s,pow);
    end if;
  end put;

-- THE OUTPUT OPERATIONS :

  procedure put ( p : in Poly; pow : in Power ) is
  begin
    put(Standard_Output,p,pow);
  end put;

  procedure put ( file : in file_type; p : in Poly; pow : in Power ) is

    n : constant natural32 := Number_of_Unknowns(p);
    standard : constant boolean := (Symbol_Table.Number < n);
    first_time : boolean := true;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      if first_time
       then first_time := false;
       else Write_Plus(file,t.cf);
      end if;
      put(file,t,standard,pow);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put(file,delimiter);
  end put;

  procedure put ( n : in natural32; p : in Poly; pow : in Power ) is
  begin
    put(Standard_Output,n,p,pow);
  end put;

  procedure put ( file : in file_type;
                  n : in natural32; p : in Poly; pow : in Power ) is
  begin
    put(file,n,3); new_line(file);
    put(file,p,pow);
  end put;

  procedure put ( p : in Poly ) is
  begin
    put(Standard_Output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly ) is
  begin
    put(file,p,'*');
  end put;

  procedure put ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put(Standard_Output,p,s);
  end put;

  procedure put ( file : in file_type;
                  p : in Poly; s : in Array_of_Symbols ) is

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      put(file,t,s,'*');
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put(file,delimiter);
  end put;

  procedure put ( n : in natural32; p : in Poly; s : in Array_of_Symbols ) is
  begin
    put(Standard_Output,n,p,s);
  end put;

  procedure put ( file : in file_type;
                  n : in natural32; p : in Poly; s : in Array_of_Symbols ) is
  begin
    put(file,n,3); new_line(file);
    put(file,p,s);
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
    put_terms(file,p,dp);
  end put_terms;

  procedure put ( file : in file_type; p : in Poly; dp : in natural32 ) is
  begin
    put(file,p);
  end put;

  procedure put_line ( p : in Poly ) is
  begin
    put_line(Standard_Output,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly ) is
  begin
    put_line(file,p,'*');
  end put_line;

  procedure put_line ( p : in Poly; pow : in Power ) is
  begin
    put_line(Standard_Output,p,pow);
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly; pow : in Power ) is

    n : constant natural32 := Number_of_Unknowns(p);
    standard : constant boolean := (Symbol_Table.Number < n);

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      put(file,"+ ");
      put(file,t,standard,pow);
      new_line(file);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put(file,delimiter);
  end put_line;

  procedure put_line ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_line(Standard_Output,p,s);
  end put_line;

  procedure put_line ( file : in file_type;
                       p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_line(file,p,s,'*');
  end put_line;

  procedure put_line
              ( p : in Poly; s : in Array_of_Symbols; pow : in Power ) is
  begin
    put_line(Standard_Output,p,s,pow);
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols; pow : in Power ) is

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      put(file,t,s,pow);
      new_line(file);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put(file,delimiter);
  end put_line;

  procedure Display_Format is
  begin
    Standard_Complex_Polynomials_io.Display_Format;
  end Display_Format;

end Multprec_Complex_Polynomials_io;
