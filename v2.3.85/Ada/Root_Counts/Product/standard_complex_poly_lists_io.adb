with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Parse_Numbers;             use Standard_Parse_Numbers;
with Standard_Natural_Vectors;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions;
with Symbol_Table;                       use Symbol_Table;
with Symbol_Table_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;

package body Standard_Complex_Poly_Lists_io is

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
    char,oper : character;
    term,res,acc : Poly := Null_Poly;

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
          if char = '^'
           then Read_Power_Factor(file,char,term);
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
            p := acc + res;
            Clear(acc); Clear(res);
            return;
           -- raise BAD_BRACKET;
          end if;
          exit;
        when '*' =>
          if res = Null_Poly then
            raise ILLEGAL_CHARACTER;
          else -- the case " ) * " :
            oper := char; get(file,char);  -- skip '*'
            Read_Term(file,bc,char,n,term);
            Skip_Spaces_and_CR(file,char);
            if char = '^'
             then Read_Power_Factor(file,char,term);
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
          end if;
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
    k : integer32 := 0;
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
          if c = Create(0.0)
           then raise ILLEGAL_CHARACTER;
           else exit;
          end if;
        when delimiter => 
          if bc /= 0
           then exit; -- raise BAD_BRACKET;
          end if;
          exit;
        when '(' =>  -- no modification of bc by Read_Term
          if c = Create(0.0) or else c = Create(-1.0)
           then c := Create(0.0); exit; -- the case "+ (" or "- (" :
           else raise BAD_BRACKET;      -- the case "c  (" :
          end if;
        when ')' =>  -- no modification of bc by Read_Term
          if bc < 0
           then exit; -- raise BAD_BRACKET;
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
    tmp.cf := c;
    tmp.dg := d;
    termp := create(tmp);
    Clear(tmp);
    if Number_Of_Unknowns(res) > 0
     then Mul(termp,res); Clear(res);
    end if;
  end Read_Term;

-- EXPORTED ROUTINES :

  procedure get ( p : out Prod_Poly ) is
  begin
    get(Standard_Input,p);
  end get;

  procedure get ( file : in file_type; p : out Prod_Poly ) is

    res,res_last : Prod_Poly;
    c : character;
    bc : integer32 := 0;
    q : Poly;

  begin
    get(file,c);
    while c /= ';' loop
      case c is
        when '(' => 
          Read_Polynomial(file,bc,q);
          Append(res,res_last,q);
        when '*' | ' ' => null;
        when others => put("read character "); put(c);
                       put(" unexpected, will skip it...");
      end case;
      get(file,c);
    end loop;
    p := res;
  end get;

  procedure put ( p : in Prod_Poly ) is
  begin
    put(Standard_Output,p);
  end put;

  procedure put ( file : in file_type; p : in Prod_Poly ) is

    tmp : Prod_Poly;
    i : Poly;
    use Lists_of_Standard_Complex_Polynomials;

  begin
    if not Is_Null(p) then
      put(file,"("); put_terms(file,Head_Of(p)); put(file,")");
      tmp := Tail_Of(p);
      while not Is_Null(tmp) loop
        i := Head_Of(tmp);
        put(file,"*("); put_terms(file,i); put(file,")");
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    put(file,";");
  end put;

  procedure put_line ( p : in Prod_Poly ) is
  begin
    put_line(Standard_Output,p);
  end put_line;

  procedure put_line ( file : in file_type; p : in Prod_Poly ) is

    tmp : Prod_Poly;
    i : Poly;
    use Lists_of_Standard_Complex_Polynomials;

  begin
    if not Is_Null(p) then
      new_line(file);
      put(file,"("); put_terms_line(file,Head_Of(p)); put(file,")");
      tmp := Tail_Of(p);
      while not Is_Null(tmp) loop
        i := Head_Of(tmp);
        put(file,"*"); new_line(file); put(file,"("); 
        put_terms_line(file,i); put(file,")");
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    put_line(file,";");
  end put_line;

end Standard_Complex_Poly_Lists_io;
