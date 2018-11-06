with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with Standard_Complex_Series;
with Standard_Complex_Series_io;        use Standard_Complex_Series_io;
with Standard_Complex_Series_Vectors;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Functions;
with DoblDobl_Complex_Series;
with DoblDobl_Complex_Series_io;        use DoblDobl_Complex_Series_io;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Functions;
with QuadDobl_Complex_Series;
with QuadDobl_Complex_Series_io;        use QuadDobl_Complex_Series_io;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Functions;
with Standard_Random_Series_Vectors;
with DoblDobl_Random_Series_Vectors;
with QuadDobl_Random_Series_Vectors;
with Standard_Polynomial_CSeries;
with DoblDobl_Polynomial_CSeries;
with QuadDobl_Polynomial_CSeries;
with Complex_Series_and_Polynomials;    use Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;

procedure ts_cserpol is

-- DESCRIPTION :
--   Tests the development of polynomials in several variables,
--   with truncated power series as coefficients.
 
  procedure Write ( s : in Standard_CSeries_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Very simple output of a polynomial with series coefficients,
  --   in standard double precision.
 
    cnt : natural32 := 0;

    procedure Visit_Term ( t : in Standard_CSeries_Polynomials.Term;
                           c : out boolean ) is
   
      cf : constant Standard_Complex_Series.Link_to_Series := t.cf;

    begin
      cnt := cnt + 1;
      put("The coefficient of term "); put(cnt); put_line(" :");
      put(cf);
      put("has degree "); put(cf.deg,1);
      put(" and degrees : "); put(t.dg.all); new_line;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new Standard_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
  end Write;
 
  procedure Write ( s : in DoblDobl_CSeries_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Very simple output of a polynomial with series coefficients,
  --   in double double precision.
 
    cnt : natural32 := 0;

    procedure Visit_Term ( t : in DoblDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is
   
      cf : constant DoblDobl_Complex_Series.Link_to_Series := t.cf;

    begin
      cnt := cnt + 1;
      put("The coefficient of term "); put(cnt); put_line(" :");
      put(cf);
      put("has degree "); put(cf.deg,1);
      put(" and degrees : "); put(t.dg.all); new_line;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new DoblDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
  end Write;
 
  procedure Write ( s : in QuadDobl_CSeries_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Very simple output of a polynomial with series coefficients,
  --   in quad double precision.
 
    cnt : natural32 := 0;

    procedure Visit_Term ( t : in QuadDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is
   
      cf : constant QuadDobl_Complex_Series.Link_to_Series := t.cf;

    begin
      cnt := cnt + 1;
      put("The coefficient of term "); put(cnt); put_line(" :");
      put(cf);
      put("has degree "); put(cf.deg,1);
      put(" and degrees : "); put(t.dg.all); new_line;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new QuadDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
  end Write;

  procedure Standard_Test_Conversion is 

  -- DESCRIPTION :
  --   Prompts the user for the number of variables.
  --   and reads in a regular polynomial in several variables,
  --   for conversion into a series polynomial,
  --   in standard double precision.

    n : natural32 := 0;
    i : integer32 := 0;
    p,q : Standard_Complex_Polynomials.Poly;
    s : Standard_CSeries_Polynomials.Poly;
    ans : character;
 
  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.init(n);
    put("Give a polynomial : "); get(p);
    put("> your polynomial : "); put(p); new_line;
    new_line;
    put("Give the index of the series parameter : "); get(i);
    new_line;
    put("Extra output during the conversion ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then s := Polynomial_to_Series_Polynomial(p,i,true);
     else s := Polynomial_to_Series_Polynomial(p,i);
    end if;
    new_line;
    put_line("The series polynomial s :");
    Write(s);
    q := Series_Polynomial_to_Polynomial(s,i,ans = 'y');
    put("s as poly : "); put(q); new_line;
  end Standard_Test_Conversion;

  procedure DoblDobl_Test_Conversion is 

  -- DESCRIPTION :
  --   Prompts the user for the number of variables.
  --   and reads in a regular polynomial in several variables,
  --   for conversion into a series polynomial,
  --   in double double precision.

    n : natural32 := 0;
    i : integer32 := 0;
    p,q : DoblDobl_Complex_Polynomials.Poly;
    s : DoblDobl_CSeries_Polynomials.Poly;
    ans : character;
 
  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.init(n);
    put("Give a polynomial : "); get(p);
    put("> your polynomial : "); put(p); new_line;
    new_line;
    put("Give the index of the series parameter : "); get(i);
    new_line;
    put("Extra output during the conversion ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then s := Polynomial_to_Series_Polynomial(p,i,true);
     else s := Polynomial_to_Series_Polynomial(p,i);
    end if;
    new_line;
    put_line("The series polynomial s :");
    Write(s);
    q := Series_Polynomial_to_Polynomial(s,i,ans = 'y');
    put("s as poly : "); put(q); new_line;
  end DoblDobl_Test_Conversion;

  procedure QuadDobl_Test_Conversion is 

  -- DESCRIPTION :
  --   Prompts the user for the number of variables.
  --   and reads in a regular polynomial in several variables,
  --   for conversion into a series polynomial,
  --   in quad double precision.

    n : natural32 := 0;
    i : integer32 := 0;
    p,q : QuadDobl_Complex_Polynomials.Poly;
    s : QuadDobl_CSeries_Polynomials.Poly;
    ans : character;
 
  begin
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.init(n);
    put("Give a polynomial : "); get(p);
    put("> your polynomial : "); put(p); new_line;
    new_line;
    put("Give the index of the series parameter : "); get(i);
    new_line;
    put("Extra output during the conversion ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then s := Polynomial_to_Series_Polynomial(p,i,true);
     else s := Polynomial_to_Series_Polynomial(p,i);
    end if;
    new_line;
    put_line("The series polynomial s :");
    Write(s);
    q := Series_Polynomial_to_Polynomial(s,i,ans = 'y');
    put("s as poly : "); put(q); new_line;
  end QuadDobl_Test_Conversion;

  function Factor ( n,k : integer32;
                    s : Standard_Complex_Series.Series )
                  return Standard_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns x[k] - s as a polynomial in n variables.
  --   All coefficients of the polynomial on return have the same degree,
  --   the same as s.deg, and are in standard double precision.

  -- REQUIRED : k is in range 1..n.

    res : Standard_CSeries_Polynomials.Poly;
    one : constant Standard_Complex_Series.Link_to_Series
        := Standard_Complex_Series.Create(1);
    trm : Standard_CSeries_Polynomials.Term;

  begin
    trm.cf := one;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := Standard_CSeries_Polynomials.Create(trm);
    trm.cf := new Standard_Complex_Series.Series'(s);
    trm.dg(k) := 0;
    Standard_CSeries_Polynomials.Sub(res,trm);
    return res;
  end Factor;

  function Factor ( n,k : integer32;
                    s : DoblDobl_Complex_Series.Series )
                  return DoblDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns x[k] - s as a polynomial in n variables.
  --   All coefficients of the polynomial on return have the same degree,
  --   the same as s.deg, and are in double double precision.

  -- REQUIRED : k is in range 1..n.

    res : DoblDobl_CSeries_Polynomials.Poly;
    one : constant DoblDobl_Complex_Series.Link_to_Series
        := DoblDobl_Complex_Series.Create(1);
    trm : DoblDobl_CSeries_Polynomials.Term;

  begin
    trm.cf := one;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := DoblDobl_CSeries_Polynomials.Create(trm);
    trm.cf := new DoblDobl_Complex_Series.Series'(s);
    trm.dg(k) := 0;
    DoblDobl_CSeries_Polynomials.Sub(res,trm);
    return res;
  end Factor;

  function Factor ( n,k : integer32;
                    s : QuadDobl_Complex_Series.Series )
                  return QuadDobl_CSeries_Polynomials.Poly is

  -- DESCRIPTION :
  --   Returns x[k] - s as a polynomial in n variables.
  --   All coefficients of the polynomial on return have the same degree,
  --   the same as s.deg, and are in quad double precision.

  -- REQUIRED : k is in range 1..n.

    res : QuadDobl_CSeries_Polynomials.Poly;
    one : constant QuadDobl_Complex_Series.Link_to_Series
        := QuadDobl_Complex_Series.Create(1);
    trm : QuadDobl_CSeries_Polynomials.Term;

  begin
    trm.cf := one;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := QuadDobl_CSeries_Polynomials.Create(trm);
    trm.cf := new QuadDobl_Complex_Series.Series'(s);
    trm.dg(k) := 0;
    QuadDobl_CSeries_Polynomials.Sub(res,trm);
    return res;
  end Factor;

  function Product ( s : Standard_Complex_Series_Vectors.Vector )
                   return Standard_CSeries_Polynomials.Poly is

  -- DESCRIPION :
  --   Returns the product of the factors x[k] - s[k], for k in s'range,
  --   where s'first = 1, in standard double precision.

    dim : constant integer32 := s'last;
    res : Standard_CSeries_Polynomials.Poly := Factor(dim,1,s(1).all);
    fac : Standard_CSeries_Polynomials.Poly;

  begin
    for k in 2..s'last loop
      fac := Factor(dim,k,s(k).all);
      Standard_CSeries_Polynomials.Mul(res,fac);
      Standard_CSeries_Polynomials.Clear(fac);
    end loop;
    return res;
  end Product;

  function Product ( s : DoblDobl_Complex_Series_Vectors.Vector )
                   return DoblDobl_CSeries_Polynomials.Poly is

  -- DESCRIPION :
  --   Returns the product of the factors x[k] - s[k], for k in s'range,
  --   where s'first = 1, in double double precision.

    dim : constant integer32 := s'last;
    res : DoblDobl_CSeries_Polynomials.Poly := Factor(dim,1,s(1).all);
    fac : DoblDobl_CSeries_Polynomials.Poly;

  begin
    for k in 2..s'last loop
      fac := Factor(dim,k,s(k).all);
      DoblDobl_CSeries_Polynomials.Mul(res,fac);
      DoblDobl_CSeries_Polynomials.Clear(fac);
    end loop;
    return res;
  end Product;

  function Product ( s : QuadDobl_Complex_Series_Vectors.Vector )
                   return QuadDobl_CSeries_Polynomials.Poly is

  -- DESCRIPION :
  --   Returns the product of the factors x[k] - s[k], for k in s'range,
  --   where s'first = 1, in double double precision.

    dim : constant integer32 := s'last;
    res : QuadDobl_CSeries_Polynomials.Poly := Factor(dim,1,s(1).all);
    fac : QuadDobl_CSeries_Polynomials.Poly;

  begin
    for k in 2..s'last loop
      fac := Factor(dim,k,s(k).all);
      QuadDobl_CSeries_Polynomials.Mul(res,fac);
      QuadDobl_CSeries_Polynomials.Clear(fac);
    end loop;
    return res;
  end Product;

  procedure Standard_Test_Evaluation is

  -- DESCRIPTION :
  --   Prompts for the number of variables and the degree of the series.
  --   Then as many random series as the number of variables are generated.
  --   The polynomial is of the product of x[k] - s[k], where k ranges
  --   over the number of variables 'x' and series 's'.
  --   So the evaluation at the series should produce zero.

    degree,dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables in the polynomial : "); get(dim);
    put("Give the degree of the power series : "); get(degree);
    declare
      rns : constant Standard_Complex_Series_Vectors.Vector(1..dim)
          := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
      pol : constant Standard_CSeries_Polynomials.Poly := Product(rns);
      eva : constant Standard_Complex_Series.Link_to_Series
          := Standard_CSeries_Poly_Functions.Eval(pol,rns);
    begin
      for i in 1..dim loop
        put("random series "); put(i,1); put_line(" :");
        put(rns(i));
      end loop;
      put_line("The polynomial :"); Write(pol);
      put_line("The value at the polynomial :"); put(eva);
    end;
  end Standard_Test_Evaluation;

  procedure DoblDobl_Test_Evaluation is

  -- DESCRIPTION :
  --   Prompts for the number of variables and the degree of the series.
  --   Then as many random series as the number of variables are generated.
  --   The polynomial is of the product of x[k] - s[k], where k ranges
  --   over the number of variables 'x' and series 's'.
  --   So the evaluation at the series should produce zero.

    degree,dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables in the polynomial : "); get(dim);
    put("Give the degree of the power series : "); get(degree);
    declare
      rns : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
          := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
      pol : constant DoblDobl_CSeries_Polynomials.Poly := Product(rns);
      eva : constant DoblDobl_Complex_Series.Link_to_Series
          := DoblDobl_CSeries_Poly_Functions.Eval(pol,rns);
    begin
      for i in 1..dim loop
        put("random series "); put(i,1); put_line(" :");
        put(rns(i));
      end loop;
      put_line("The polynomial :"); Write(pol);
      put_line("The value at the polynomial :"); put(eva);
    end;
  end DoblDobl_Test_Evaluation;

  procedure QuadDobl_Test_Evaluation is

  -- DESCRIPTION :
  --   Prompts for the number of variables and the degree  of the series.
  --   Then as many random series as the number of variables are generated.
  --   The polynomial is of the product of x[k] - s[k], where k ranges
  --   over the number of variables 'x' and series 's'.
  --   So the evaluation at the series should produce zero.

    degree,dim : integer32 := 0;

  begin
    new_line;
    put("Give the number of variables in the polynomial : "); get(dim);
    put("Give the degree of the power series : "); get(degree);
    declare
      rns : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
          := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
      pol : constant QuadDobl_CSeries_Polynomials.Poly := Product(rns);
      eva : constant QuadDobl_Complex_Series.Link_to_Series
          := QuadDobl_CSeries_Poly_Functions.Eval(pol,rns);
    begin
      for i in 1..dim loop
        put("random series "); put(i,1); put_line(" :");
        put(rns(i));
      end loop;
      put_line("The polynomial :"); Write(pol);
      put_line("The value at the polynomial :"); put(eva);
    end;
  end QuadDobl_Test_Evaluation;

  procedure Standard_Test_Input_Output is

  -- DESCRIPTION :
  --   Reads a series in symbolic format and writes the series back,
  --   in standard double precision.

    maxdeg : constant integer32 := 16;
    s : Standard_Complex_Series.Series(maxdeg);

  begin
    new_line;
    put_line("Give a series, terminate with ;");
    Complex_Series_and_Polynomials_io.get(s);
    new_line;
    put_line("The coefficients of the series : ");
    put(s);
    new_line;
    put_line("The series : ");
    Complex_Series_and_Polynomials_io.put(s);
    new_line;
  end Standard_Test_Input_Output;

  procedure DoblDobl_Test_Input_Output is

  -- DESCRIPTION :
  --   Reads a series in symbolic format and writes the series back,
  --   in double double precision.

    maxdeg : constant integer32 := 16;
    s : DoblDobl_Complex_Series.Series(maxdeg);

  begin
    new_line;
    put_line("Give a series, terminate with ;");
    Complex_Series_and_Polynomials_io.get(s);
    new_line;
    put_line("The coefficients of the series : ");
    put(s);
    new_line;
    put_line("The series : ");
    Complex_Series_and_Polynomials_io.put(s);
    new_line;
  end DoblDobl_Test_Input_Output;

  procedure QuadDobl_Test_Input_Output is

  -- DESCRIPTION :
  --   Reads a series in symbolic format and writes the series back,
  --   in quad double precision.

    maxdeg : constant integer32 := 16;
    s : QuadDobl_Complex_Series.Series(maxdeg);

  begin
    new_line;
    put_line("Give a series, terminate with ;");
    Complex_Series_and_Polynomials_io.get(s);
    new_line;
    put_line("The coefficients of the series : ");
    put(s);
    new_line;
    put_line("The series : ");
    Complex_Series_and_Polynomials_io.put(s);
    new_line;
  end QuadDobl_Test_Input_Output;

  procedure Standard_Test_Polynomial_Series is

  -- DESCRIPTION :
  --   Reads a polynomial in several variables and converts
  --   the polynomial into a polynomial series, with complex
  --   coefficients in standard double precision.

    maxdeg : constant integer32 := 16;
    sp : Standard_CSeries_Polynomials.Poly;
    ps : Standard_Polynomial_CSeries.Poly(maxdeg);
    n : natural32 := 0;

  begin
    new_line;
    put("Give the total number of symbols : "); get(n);
    Symbol_Table.init(n);
    new_line;
    put_line("Reading a series polynomial, series parameter comes first.");
    Complex_Series_and_Polynomials_io.get(sp,1);
    new_line;
    put_line("Your polynomial :");
    Complex_Series_and_Polynomials_io.put(sp,1);
    ps := Standard_Polynomial_CSeries.Create(sp);
    Standard_CSeries_Polynomials.Clear(sp);
    sp := Standard_Polynomial_CSeries.Create(ps);
    new_line;
    put_line("The polynomial from the created polynomial series :");
    Complex_Series_and_Polynomials_io.put(sp,1);
  end Standard_Test_Polynomial_Series;

  procedure DoblDobl_Test_Polynomial_Series is

  -- DESCRIPTION :
  --   Reads a polynomial in several variables and converts
  --   the polynomial into a polynomial series, with complex
  --   coefficients in double double precision.

    maxdeg : constant integer32 := 16;
    sp : DoblDobl_CSeries_Polynomials.Poly;
    ps : DoblDobl_Polynomial_CSeries.Poly(maxdeg);
    n : natural32 := 0;

  begin
    new_line;
    put("Give the total number of symbols : "); get(n);
    Symbol_Table.init(n);
    new_line;
    put_line("Reading a series polynomial, series parameter comes first.");
    Complex_Series_and_Polynomials_io.get(sp,1);
    new_line;
    put_line("Your polynomial :");
    Complex_Series_and_Polynomials_io.put(sp,1);
    ps := DoblDobl_Polynomial_CSeries.Create(sp);
    DoblDobl_CSeries_Polynomials.Clear(sp);
    sp := DoblDobl_Polynomial_CSeries.Create(ps);
    new_line;
    put_line("The polynomial from the created polynomial series :");
    Complex_Series_and_Polynomials_io.put(sp,1);
  end DoblDobl_Test_Polynomial_Series;

  procedure QuadDobl_Test_Polynomial_Series is

  -- DESCRIPTION :
  --   Reads a polynomial in several variables and converts
  --   the polynomial into a polynomial series, with complex
  --   coefficients in double double precision.

    maxdeg : constant integer32 := 16;
    sp : QuadDobl_CSeries_Polynomials.Poly;
    ps : QuadDobl_Polynomial_CSeries.Poly(maxdeg);
    n : natural32 := 0;

  begin
    new_line;
    put("Give the total number of symbols : "); get(n);
    Symbol_Table.init(n);
    new_line;
    put_line("Reading a series polynomial, series parameter comes first.");
    Complex_Series_and_Polynomials_io.get(sp,1);
    new_line;
    put_line("Your polynomial :");
    Complex_Series_and_Polynomials_io.put(sp,1);
    ps := QuadDobl_Polynomial_CSeries.Create(sp);
    QuadDobl_CSeries_Polynomials.Clear(sp);
    sp := QuadDobl_Polynomial_Cseries.Create(ps);
    new_line;
    put_line("The polynomial from the created polynomial series :");
    Complex_Series_and_Polynomials_io.put(sp,1);
  end QuadDobl_Test_Polynomial_Series;
  
  procedure Main is

  -- DESCRIPTION :
  --   Displays the menu of tests
  --   and prompts the user to make a choice.

    ans,prc : character;

  begin
    new_line;
    put_line("MENU to test polynomials with series coefficients :");
    put_line("  0. test conversion from/to ordinary polynomials;");
    put_line("  1. test evaluation at power series;");
    put_line("  2. test symbolic input and output;");
    put_line("  3. test conversion to polynomial series.");
    put("Type 0, 1, 2, or 3 to select a test : ");
    Ask_Alternative(ans,"0123");
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision;");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(prc,"012");
    case ans is
      when '0' =>
        case prc is
          when '0' => Standard_Test_Conversion;
          when '1' => DoblDobl_Test_Conversion;
          when '2' => QuadDobl_Test_Conversion;
          when others => null;
        end case;
      when '1' =>
        case prc is
          when '0' => Standard_Test_Evaluation;
          when '1' => DoblDobl_Test_Evaluation;
          when '2' => QuadDobl_Test_Evaluation;
          when others => null;
        end case;
      when '2' =>
        case prc is
          when '0' => Standard_Test_Input_Output;
          when '1' => DoblDobl_Test_Input_Output;
          when '2' => QuadDobl_Test_Input_Output;
          when others => null;
        end case;
      when '3' =>
        case prc is
          when '0' => Standard_Test_Polynomial_Series;
          when '1' => DoblDobl_Test_Polynomial_Series;
          when '2' => QuadDobl_Test_Polynomial_Series;
          when others => null;
        end case;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cserpol;
