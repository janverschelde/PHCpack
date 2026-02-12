with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Write_Numbers;             use Standard_Write_Numbers;
with Symbol_Table_io;
with Line_Breaks;                        use Line_Breaks;
with Standard_Complex_Laur_Readers;

package body Standard_Complex_Laurentials_io is

-- THE INPUT OPERATIONS :

  procedure get ( p : out Poly ) is
  begin
    get(standard_input,p);
  end get;

  procedure get ( file : in file_type; p : out Poly ) is

    bc : integer32 := 0;

  begin
    Standard_Complex_Laur_Readers.Read_Polynomial(file,bc,p);
  end get;
 
  procedure get ( n : out natural32; p : out Poly ) is
  begin
    get(standard_input,n,p);
  end get;

  procedure get ( file : in file_type; n : out natural32; p : out Poly ) is
  begin
    n := 0;
    get(file,n);
    if Symbol_Table.Empty
     then Symbol_Table.Init(n);
    end if;
    get(file,p);
  end get;

-- THE OUTPUT OPERATIONS :

  function All_Zeroes ( d : Degrees ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all d(i) = 0, for all i, otherwise false is returned.

  begin
    for i in d'range loop
      if d(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end All_Zeroes;

  procedure put ( file : in file_type; d,i : in integer32;
                  std : in boolean; pow : in Power ) is
  begin
    if std then
      put(file,'x');
      if i < 10
       then put(file,i,1);
       else put(file,i,2);
      end if;
    else
      Symbol_Table_io.put(file,Symbol_Table.get(natural32(i)));
    end if;
    if d > 1 or d < 0 then
      if pow = '^'
       then put(file,'^');
       else put(file,"**");
      end if;
      if d < 10
       then put(file,d,1);
       else put(file,d,2);
      end if;
    end if;
  end put;

  procedure put ( file : in file_type; d : in integer32;
                  sb : in Symbol; pow : in Power ) is
  begin
    Symbol_Table_io.put(file,sb);
    if d > 1 or d < 0 then
      if pow = '^'
       then put(file,'^');
       else put(file,"**");
      end if;
      if d < 10
       then put(file,d,1);
       else put(file,d,2);
      end if;
    end if;
  end put;

  procedure put ( file : in file_type; d : in Degrees;
                  std : in boolean; pow : in Power ) is

    first : boolean := true;

  begin
    for i in d'range loop
      if d(i) /= 0 then
        if first
         then first := false;
         else put(file,'*');
        end if;
        put(file,d(i),i,std,pow);
      end if;
    end loop;
  end put;

  procedure put ( file : in file_type; d : in Degrees;
                  s : in Array_of_Symbols; pow : in Power ) is

    first : boolean := true;

  begin
    for i in d'range loop
      if d(i) /= 0 then
        if first
         then first := false;
         else put(file,'*');
        end if;
        put(file,d(i),s(i),pow);
      end if;
    end loop;
  end put;

  procedure put ( t : in Term ) is
  begin
    put(standard_output,t);
  end put;

  procedure put ( file : in file_type; t : in Term ) is

    std : constant boolean := (Symbol_Table.Number < natural32(t.dg'last));
    nc : natural32;

  begin
    Write_Plus(file,t.cf,nc);
    Write_Number(file,t.cf,nc);
    if not All_Zeroes(t.dg) then
      put(file,'*');
      put(file,t.dg,std,'^');
    end if;
  end put;

  procedure put ( t : in Term; s : in Array_of_Symbols ) is
  begin
    put(standard_output,t,s);
  end put;

  procedure put ( file : in file_type; t : in Term;
                  s : in Array_of_Symbols ) is

    nc : natural32;

  begin
    Write_Plus(file,t.cf,nc);
    Write_Number(file,t.cf,nc);
    if not All_Zeroes(t.dg) then
      put(file,'*');
      put(file,t.dg,s,'^');
    end if;
  end put;

  procedure put ( p : in Poly ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put(standard_output,p);
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

    begin
      if first_time 
       then first_time := false;
       else Write_Plus(file,t.cf,nc); sc := sc + nc;
      end if;
      if All_Zeroes(t.dg) then
        Write_Number(file,t.cf,nc); sc := sc + nc;
      else
        Write_Coefficient(file,t.cf,nc); sc := sc + nc;
        passed := false;
        for i in t.dg'range loop
          if t.dg(i) /= 0 then
            if passed
             then put(file,'*'); sc := sc + 1;
             else passed := true;
            end if;
            sc := sc + Length(t.dg(i),i,standard,pow);
            put(file,t.dg(i),i,standard,pow);
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

  procedure put ( file : in file_type; p : in Poly ) is
  begin
    put_terms(file,p,'^');
    Line(file,1); put(file,delimiter);
  end put;

  procedure put ( file : in file_type; p : in Poly;
                  s : in Array_of_Symbols ) is
  begin
    null;
  end put;

  procedure put_line ( p : in Poly ) is
  begin
    put_line(standard_output,p);
  end put_line;

  procedure put_line ( p : in Poly; s : in Array_of_Symbols ) is
  begin
    put_line(standard_output,p,s);
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly ) is

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      new_line(file);
      put(file,t);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put_line(file,";");
  end put_line;

  procedure put_line ( file : in file_type; p : in Poly;
                       s : in Array_of_Symbols ) is

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      new_line(file);
      put(file,t,s);
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put_line(file,";");
  end put_line;

  procedure Display_Format is

    s : array(1..24) of string(1..65);

  begin
    s( 1):="  A complex multivariate  Laurent polynomial  is  a  sequence  of";
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
    s(12):="be followed by an integer number, eventually in round brackets.  ";
    s(13):="If the power equals one, then it may be omitted.                 ";
    s(14):="  An unknown can be denoted by at most 80 characters.  The  first";
    s(15):="character  must  be a letter and the other two characters must be";
    s(16):="different from '+', '-', '*', '^', '/', ';', '('  and  ')'.   The";
    s(17):="letter i means sqrt(-1), whence it does not represent an unknown.";
    s(18):="The number of unknowns may not  exceed  the  declared  dimension.";
    s(19):="  Some  examples  of  valid  notations  of  complex  multivariate";
    s(20):="polynomials:                                                     ";
    s(21):="  x**2*y + 1/2*z*y**2 - 2*z + y**3 + x - 1E9/-8.E-6* y + 3;      ";
    s(22):="  x^2*y + z*y^2 - 2*z + y^3 + x^(-1) - y**(-2) + 3;              ";
    s(23):="  (1.01 + 2.8*i)*x1**2*x2 + x3**2*x1 - 3*x1 + 2*x2*x3 - 3;       ";
    s(24):="  (x1^2*x2 + x3^2*x1 - 3*x1 + 2*x2*x3 - 3)*x2**2*(x2-1+i);       ";
    for i in s'range loop
      put_line(s(i));
    end loop;
  end Display_Format;

end Standard_Complex_Laurentials_io;
