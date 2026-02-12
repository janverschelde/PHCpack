with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers; 
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Parse_Numbers;             use Standard_Parse_Numbers;
with Standard_Integer_Vectors;
with Symbol_Table;                       use Symbol_Table;
with Symbol_Table_io;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions;

package body Standard_Complex_Laur_Readers is

  procedure Read_Power_Factor
              ( file : in file_type; char : in out character;
                p : in out Poly ) is

    d,ne,ne2 : natural32 := 0;
    expo,expo2 : integer32 := 1;
    sign : character;

  begin
    get(file,char);                           -- skip the '^'
    Parse_also_Brackets(file,char,expo,expo2,ne,ne2,sign);
    if expo < 0
     then raise ILLEGAL_OPERATION;
    end if;
    d := natural32(expo);
    Pow(p,d);
  exception
    when others => put_line("exception raised in Read_Power_Factor"); raise;
  end Read_Power_Factor;

  procedure Read_Polynomial
              ( file : in file_type; bc : in out integer32;
                p : out Poly ) is

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
    end if;                        -- the first term can have no sign
    Read_Term(file,bc,char,n,res); -- therefore read it first
    loop
      case char is
        when ' ' | ASCII.CR =>
          get(file,char);    -- skip blanks
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
            when '+' => Add(acc,res); Clear(res); Copy(term,res);
            when '-' => Add(acc,res); Clear(res); Copy(term,res); Min(res);
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
         -- if res = Null_Poly then    -- an input like (0.0)*(2+x)
         --   raise ILLEGAL_CHARACTER; -- should not raise an error
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
  exception
    when others => put_line("exception raised in Read_Polynomial"); raise;
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
      Parse_also_Brackets(file,char,expo,expo2,ne,ne2,sign);
      d(k) := d(k) + expo;
      Skip_Spaces_and_CR(file,char); 
      if char /= '*'                            -- the case x^2*...
       then return;                             -- end of factor
       else get(file,char);                     -- skip the '*'
      end if; 
    elsif char = '*' then
      get(file,char);
      if char = '*' then
        get(file,char);                 -- the case " x ** expo "
        Parse_also_Brackets(file,char,expo,expo2,ne,ne2,sign);
        d(k) := d(k) + expo;
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
     then Read_Term(file,bc,char,n,pb); -- the case " x * c " or " x ** c * c "
     else Read_Factor(file,bc,char,n,d,pb); -- the case " x * y "
    end if;
  exception
    when others => put_line("exception raised in Read_Factor"); raise;
  end Read_Factor;

  procedure Read_Term ( file : in file_type; bc : in out integer32;
                        char : in out character;
                        n : in natural32; termp : in out Poly ) is

    c : Complex_Number;
    d : Degrees := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
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
  exception
    when others => put_line("exception raised in Read_Term"); raise;
  end Read_Term;

end Standard_Complex_Laur_Readers;
