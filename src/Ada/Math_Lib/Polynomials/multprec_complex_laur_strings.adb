with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Strings_and_Numbers;                use Strings_and_Numbers;
with Standard_Parse_Numbers;
with Multprec_Parse_Numbers;             use Multprec_Parse_Numbers;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Number_Tools;      use Multprec_Complex_Number_Tools;
with Standard_Natural_Vectors;
with Standard_Integer_Vectors;
with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Parse_Polynomial_Exceptions;        use Parse_Polynomial_Exceptions;
with Standard_Complex_Poly_Strings;
--with Standard_Complex_Laur_Strings;

package body Multprec_Complex_Laur_Strings is

-- AUXILIARIES FOR INPUT :

  procedure Parse_Term
              ( s : in string; size : in natural32; bc : in out integer32;
                p : in out integer;
                n : in natural32; termp : in out Poly );

  -- DESCRIPTION :
  --   Parses the string for a term, starting at s(p).

  -- ON ENTRY :
  --   s        string of characters to contain a polynomial;
  --   size     working precision to evaluate fractions;
  --   bc       counts #open - #closed brackets;
  --   p        current position in the string s;
  --   n        number of variables;
  --   termp    accumulating polynomial for the term.

  -- ON RETURN :
  --   bc       updated bracket counter;
  --   p        updated position in the string;
  --   termp    resulting term read from string.

  procedure Parse_Factor
              ( s : in string; size : in natural32;
                bc : in out integer32; p : in out integer; n : in natural32;
                d : in out Degrees; pb : in out Poly );

  -- DESCRIPTION :
  --   Parses the string s for a factor, starting at s(p).

  -- ON ENTRY :
  --   s        string of characters to contain a polynomial;
  --   size     sets the working precision to evaluate fractions;
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
              ( s : in string; size : in natural32;
                k : in out integer; p : in out Poly ) is

  -- DESCRIPTION :
  --   Reads the exponent of a factor stored in s and returns
  --   in p the polynomial p raised to the read exponent.

  -- ON ENTRY :
  --   s        string to contain a polynomial
  --   size     for the working precision to evaluate fractions;
  --   k        s(k) equals the exponentiation symbol '^';
  --   p        current factor.

  -- ON RETURN :
  --   k        updated position in the string;
  --   p        p raised to the read exponent.
 
    nd : natural32 := 0;
    d : integer32;
    sign : character;

  begin
    k := k + 1;                -- skip the '^'
    if k > s'last
     then return;
    end if;
   -- Standard_Complex_Laur_Strings.Read_Exponent(s,k,d);
    Standard_Parse_Numbers.Parse_also_Brackets(s,k,d,nd,sign);
    Pow(p,natural32(d));
   -- skip last digit of exponent 
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,k); 
  exception
    when others =>
      put("Exception raised at character "); put(integer32(k),1);
      put_line(" of " & s & " in Parse_Power_Factor."); raise;
  end Parse_Power_Factor;

  procedure Parse_Polynomial
               ( s : in string; size : in natural32; bc : in out integer32;
                 k : in out integer;
                 n : in natural32; p : in out Poly ) is

  -- DESCRIPTION :
  --   Main parsing routine of a string s into a polynomial in n variables,
  --   starting at the current position k,
  --   while keeping track of the bracket counter bc.

    oper : character;
    term,res,acc : Poly := Null_Poly;

  begin
    oper := '+';
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,k);
    if k > s'last
     then return;
    end if;
    if s(k) = '-'
     then oper := '-';
    end if;                         -- the first term can have no sign
    Parse_Term(s,size,bc,k,n,res);  -- therefore read it first
    loop
      case s(k) is
        when ' ' | ASCII.CR | ASCII.LF => k := k + 1;   -- skip blanks
        when '+' | '-' =>
          oper := s(k);
          Parse_Term(s,size,bc,k,n,term);
          Add(res,term); Clear(term);
        when delimiter =>
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' =>
          bc := bc + 1;
          k := k + 1;
          Parse_Polynomial(s(k..s'last),size,bc,k,n,term);
          Standard_Parse_Numbers.Skip_Spaces_and_CR(s,k);
          if s(k) = '(' -- or s(k) = ')'
           then raise BAD_BRACKET;
          end if;
          if s(k) = '^' then
            Parse_Power_Factor(s,size,k,term);
          elsif s(k) = '*' then
            if s(k+1) = '*' then
              k := k + 1;
              Parse_Power_Factor(s,size,k,term);
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
          -- if res = Null_Poly then     -- zero coefficients should
          --   raise ILLEGAL_CHARACTER;  -- not raise an exception
          -- else -- the case " ) * " :
            oper := s(k); k := k + 1;  -- skip '*'
            Parse_Term(s,size,bc,k,n,term);
            Standard_Parse_Numbers.Skip_Spaces_and_CR(s,k);
            if s(k) = '^' then
              Parse_Power_Factor(s,size,k,term);
            elsif s(k) = '*' then
              if s(k+1) = '*' then
                k := k + 1;
                Parse_Power_Factor(s,size,k,term);
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
           else Parse_Power_Factor(s,size,k,res);
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

  procedure Parse_Factor
              ( s : in string; size : in natural32;
                bc : in out integer32; p : in out integer; n : in natural32;
                d : in out Degrees; pb : in out Poly ) is

    sb : symbol;
    nd : natural32;
    k,expo : integer32;
    sign : character;
 
  begin
    if p > s'last
     then return;
    end if;
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if s(p) = '(' then 
      bc := bc + 1;
      p := p + 1;        -- get a new symbol, skip '('
      if p > s'last
       then return;
      end if;
      Parse_Polynomial(s(p..s'last),size,bc,p,n,pb);
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if s(p) = '^' then
        Parse_Power_Factor(s,size,p,pb);
      elsif s(p) = '*' then
        if s(p+1) = '*' then
          p := p + 1;
          Parse_Power_Factor(s,size,p,pb);
        end if;
      end if;
      return;
    end if;
    Symbol_Table_io.Parse_Symbol(s,p,sb);
    k := integer32(Symbol_Table.Check_Symbol(n,sb));
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if s(p) = '^' then
      p := p + 1;                               -- skip the '^'
      if p > s'last
       then return;
      end if;
     -- Standard_Complex_Laur_Strings.Read_Exponent(s,p,expo);
      Standard_Parse_Numbers.Parse_also_Brackets(s,p,expo,nd,sign);
      d(k) := d(k) + expo;
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if s(p) /= '*'                            -- the case x^2*...
       then return;                             -- end of factor
       else p := p + 1;                         -- skip the '*'
      end if; 
    elsif s(p) = '*' then
      p := p + 1;
      if s(p) = '*' then
        p := p + 1;                             -- the case " x ** expo "
        if p > s'last
         then return;
        end if;
       -- Standard_Complex_Laur_Strings.Read_Exponent(s,p,expo);
        Standard_Parse_Numbers.Parse_also_Brackets(s,p,expo,nd,sign);
        d(k) := d(k) + expo;
        Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
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
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if (s(p) = '-') or (s(p) = '+') 
     then return;
    end if;
    if Convert(s(p)) < 10
     then Parse_Term(s,size,bc,p,n,pb); -- the case " x * c " or " x ** c * c "
     else Parse_Factor(s,size,bc,p,n,d,pb); -- the case " x * y "
    end if;
  exception
    when others =>
      put("Exception raised at character "); put(integer32(p),1);
      put_line(" of " & s & " in Parse_Factor."); raise;
  end Parse_Factor;
 
  procedure Parse_Term
              ( s : in string; size : in natural32; bc : in out integer32;
                p : in out integer;
                n : in natural32; termp : in out Poly ) is

    c : Complex_Number;
    d : Degrees := new Standard_Integer_Vectors.Vector'(1..integer32(n) => 0);
    pb,res : Poly := Null_Poly;
    tmp : Term;
    realzero : constant Floating_Number := Create(integer(0));
    realmin1 : constant Floating_Number := Create(integer(-1));
    compzero : constant Complex_Number := Create(realzero);
    compmin1 : constant Complex_Number := Create(realmin1);

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
    Parse(s,size,p,c); -- look for 'i' :
    Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
    if ( Equal(c,compzero) ) and then ((s(p) = 'i') or (s(p) = 'I')) then
      -- the case "+ i" :
      c := Create(0.0,1.0); 
      p := p + 1;        -- skip 'i'
    elsif ( Equal(c,compmin1) ) and then ((s(p) = 'i') or (s(p) = 'I')) then
      -- the case "- i" :
      c := Create(0.0,-1.0);
      p := p + 1;      -- skip 'i'
    elsif s(p) = '*' then -- the case ".. c *.." :
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      p := p + 1;  -- skip '*'
      Standard_Parse_Numbers.Skip_Spaces_and_CR(s,p);
      if (s(p) = 'i') or (s(p) = 'I') then -- the case ".. c * i.." :
        c := Multprec_Complex_Numbers.Create(realzero,REAL_PART(c));
        p := p + 1;    -- skip 'i'
      else                                 -- the case ".. c * x.." :
        Parse_Factor(s,size,bc,p,n,d,pb);
        if pb /= Null_Poly
         then Clear(res); Copy(pb,res); Clear(pb);
        end if;
      end if;
    end if; -- the case ".. c ?"  will be treated in the loop
    loop
      case s(p) is
        when ' ' | ASCII.CR | ASCII.LF => p := p + 1;
        when '*' =>
          p := p + 1;
          Parse_Factor(s,size,bc,p,n,d,pb);
          Collect_Factor_Polynomial;
        when '+' | '-' => 
         -- if Equal(c,compzero)
         --  then raise ILLEGAL_CHARACTER;
         --  else exit;
         -- end if;
          exit; -- terms with zero coefficient are allowed
        when delimiter => 
          if bc /= 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when '(' => 
          if Equal(c,compzero) or else Equal(c,compmin1)
           then c := Create(0.0); exit; -- the case "+ (" or "- (" :
           else raise BAD_BRACKET;      -- the case "c  (" :
          end if;
        when ')' =>
          if bc < 0
           then raise BAD_BRACKET;
          end if;
          exit;
        when others  =>
          if Equal(c,compzero) then
            c := Create(1.0);
            Parse_Factor(s,size,bc,p,n,d,pb);
          elsif Equal(c,compmin1) then
            Parse_Factor(s,size,bc,p,n,d,pb);
          elsif s(p) = '^' then
            Parse_Power_Factor(s,size,p,res);
          elsif s(p) = '*' then
            if s(p+1) = '*' then
              p := p + 1;
              Parse_Power_Factor(s,size,p,res);
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
      Clear(tmp);
      if Number_Of_Unknowns(res) > 0
       then Mul(termp,res); Clear(res);
      end if;
    end if;
  exception
    when others =>
      put("Exception raised at character "); put(integer32(p),1);
      put_line(" of " & s & " in Parse_Term."); raise;
  end Parse_Term;

-- AUXILIARIES FOR OUTPUT :

  function Write ( t : Term; lead : boolean ) return string is

  -- DESCRIPTION :
  --   Returns the string representation of the term.
  --   When the term is the leading term of the polynomial,
  --   then lead is true.

    function Rec_Write ( k : integer32; accu : string; first : boolean )
                       return string is
    begin
      if k > t.dg'last then
        return accu;
      elsif t.dg(k) = 0 then
        return Rec_Write(k+1,accu,first);
      else
        declare
          sb : constant Symbol := Symbol_Table.Get(natural32(k));
        begin
          if t.dg(k) > 1 or t.dg(k) < 0 then
            if first then
              declare
                new_accu : constant string
                  := Standard_Complex_Poly_Strings.Concat_Symbol0(accu,sb)
                     & '^' & Convert(t.dg(k));
              begin
                return Rec_Write(k+1,new_accu,false);
              end;
            else
              declare
                new_accu : constant string
                  := Standard_Complex_Poly_Strings.Concat_Symbol1(accu,sb)
                     & '^' & Convert(t.dg(k));
              begin
                return Rec_Write(k+1,new_accu,first);
              end;
            end if;
          elsif first then
            declare
              new_accu : constant string
                := Standard_Complex_Poly_Strings.Concat_Symbol0(accu,sb);
            begin
              return Rec_Write(k+1,new_accu,false);
            end;
          else
            declare
              new_accu : constant string
                := Standard_Complex_Poly_Strings.Concat_Symbol1(accu,sb);
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

  function Write ( p : Laur_Sys ) return string is
  begin
    if p'first = p'last
     then return Write(p(p'first));
     else return Write(p(p'first)) & Write(p(p'first+1..p'last));
    end if;
  end Write;

-- TARGET FUNCTIONS :

  procedure Parse ( s : in string; k : in out integer;
                    n,size : in natural32; p : in out Poly ) is

    bc : integer32 := 0;

  begin
    Parse_Polynomial(s,size,bc,k,n,p);
  end Parse;

  function Parse ( n,size : natural32; s : string ) return Poly is

    res : Poly;
    p : integer := s'first;

  begin
    Parse(s,p,n,size,res);
    return res;
  end Parse;

  function Parse ( n,m,size : natural32; s : string ) return Laur_Sys is

    res : Laur_Sys(1..integer32(n));
    ind : constant Standard_Natural_Vectors.Vector(1..integer32(n))
        := Standard_Complex_Poly_Strings.Delimiters(n,s);

  begin
   -- This function may be called when the symbol table has not yet
   -- been initialized and then it should not crash.
    if Symbol_Table.Number < m then
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    res(1) := Parse(m,size,s(s'first..integer(ind(1))));
    for i in 2..integer32(n) loop
      res(i) := Parse(m,size,s(integer(ind(i-1)+1)..integer(ind(i))));
    end loop;
    return res;
  end Parse;

  function Parse ( m,size : natural32;
                   s : Array_of_Strings ) return Laur_Sys is

    res : Laur_Sys(integer32(s'first)..integer32(s'last));
 
  begin
    if Symbol_Table.Number < m then -- same comment as other Parse
      if not Symbol_Table.Empty
       then Symbol_Table.Clear;
      end if;
      Symbol_Table.Init(m);
    end if;
    for i in s'range loop
      res(integer32(i)) := Parse(m,size,s(i).all);
    end loop;
    return res;
  end Parse;

  function Size_Limit ( p : Poly ) return natural32 is

    nbtrm : constant natural64 := natural64(Number_of_Terms(p));
    nbvar : constant natural64 := natural64(Number_of_Unknowns(p));
    symsz : constant natural64 := 5;
    cffsz : constant natural64 := 200;
    bound : constant natural64 := 2**31 - 1;
    res : constant natural64 := nbtrm*nbvar*symsz*cffsz;

  begin
    if res > bound
     then return natural32(bound);
     else return natural32(res);
    end if;
  end Size_Limit;

  function Write ( p : Poly ) return string is

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
         then return Write(terms(tfirst),true);  -- lead term
         else return Write(terms(tfirst),false); -- trailing term
        end if;
      else
        declare
          middle : constant natural32 := (tfirst+tlast)/2;
          midterm : constant string := Write(terms(middle),false);
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

  function Write ( p : Laur_Sys ) return Array_of_Strings is

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

end Multprec_Complex_Laur_Strings;
