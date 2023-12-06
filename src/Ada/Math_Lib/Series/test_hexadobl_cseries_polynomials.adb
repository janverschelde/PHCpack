with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Symbol_Table;
with HexaDobl_Complex_Polynomials;
with HexaDobl_Complex_Polynomials_io;   use HexaDobl_Complex_Polynomials_io;
with HexaDobl_Complex_Series_io;        use HexaDobl_Complex_Series_io;
with HexaDobl_CSeries_Poly_Functions;
with HexaDobl_Random_Series_Vectors;
with HexaDobl_Polynomial_CSeries;
with Complex_Series_and_Polynomials;    use Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;

package body Test_HexaDobl_CSeries_Polynomials is

  procedure Write ( s : in HexaDobl_CSeries_Polynomials.Poly ) is

    cnt : natural32 := 0;

    procedure Visit_Term ( t : in HexaDobl_CSeries_Polynomials.Term;
                           c : out boolean ) is
   
      cf : constant HexaDobl_Complex_Series.Link_to_Series := t.cf;

    begin
      cnt := cnt + 1;
      put("The coefficient of term "); put(cnt); put_line(" :");
      put(cf);
      put("has degree "); put(cf.deg,1);
      put(" and degrees : "); put(t.dg.all); new_line;
      c := true;
    end Visit_Term;
    procedure Visit_Terms is
      new HexaDobl_CSeries_Polynomials.Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(s);
  end Write;
 
  procedure HexaDobl_Test_Conversion is 

    n : natural32 := 0;
    i : integer32 := 0;
    p,q : HexaDobl_Complex_Polynomials.Poly;
    s : HexaDobl_CSeries_Polynomials.Poly;
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
  end HexaDobl_Test_Conversion;

  function Factor ( n,k : integer32;
                    s : HexaDobl_Complex_Series.Series )
                  return HexaDobl_CSeries_Polynomials.Poly is

    res : HexaDobl_CSeries_Polynomials.Poly;
    one : constant HexaDobl_Complex_Series.Link_to_Series
        := HexaDobl_Complex_Series.Create(1);
    trm : HexaDobl_CSeries_Polynomials.Term;

  begin
    trm.cf := one;
    trm.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    trm.dg(k) := 1;
    res := HexaDobl_CSeries_Polynomials.Create(trm);
    trm.cf := new HexaDobl_Complex_Series.Series'(s);
    trm.dg(k) := 0;
    HexaDobl_CSeries_Polynomials.Sub(res,trm);
    return res;
  end Factor;

  function Product ( s : HexaDobl_Complex_Series_Vectors.Vector )
                   return HexaDobl_CSeries_Polynomials.Poly is

    dim : constant integer32 := s'last;
    res : HexaDobl_CSeries_Polynomials.Poly := Factor(dim,1,s(1).all);
    fac : HexaDobl_CSeries_Polynomials.Poly;

  begin
    for k in 2..s'last loop
      fac := Factor(dim,k,s(k).all);
      HexaDobl_CSeries_Polynomials.Mul(res,fac);
      HexaDobl_CSeries_Polynomials.Clear(fac);
    end loop;
    return res;
  end Product;

  procedure HexaDobl_Frequent_Evaluation 
              ( dim,deg : in integer32;
                pol : in HexaDobl_CSeries_Polynomials.Poly ) is

    freq : natural32 := 0;
    rns : HexaDobl_Complex_Series_Vectors.Vector(1..dim);
    eva : HexaDobl_Complex_Series.Link_to_Series;

  begin
    new_line;
    put("Give the number of evaluations : "); get(freq);
    new_line;
    put("Evaluating "); put(freq,1); put_line(" times ...");
    for i in 1..freq loop
      rns := HexaDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
      eva := HexaDobl_CSeries_Poly_Functions.Eval(pol,rns);
      HexaDobl_Complex_Series_Vectors.Clear(rns);
      HexaDobl_Complex_Series.Clear(eva);
    end loop;
  end HexaDobl_Frequent_Evaluation;

  procedure HexaDobl_Test_Evaluation is

    degree,dim : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the number of variables in the polynomial : "); get(dim);
    put("Give the degree of the power series : "); get(degree);
    declare
      rns : constant HexaDobl_Complex_Series_Vectors.Vector(1..dim)
          := HexaDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,degree);
      pol : constant HexaDobl_CSeries_Polynomials.Poly := Product(rns);
      eva : constant HexaDobl_Complex_Series.Link_to_Series
          := HexaDobl_CSeries_Poly_Functions.Eval(pol,rns);
    begin
      for i in 1..dim loop
        put("random series "); put(i,1); put_line(" :");
        put(rns(i));
      end loop;
      put_line("The polynomial :"); Write(pol);
      put_line("The value at the polynomial :"); put(eva);
      new_line; 
      put("Run a frequency test on memory management ? ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then HexaDobl_Frequent_Evaluation(dim,degree,pol);
      end if;
    end;
  end HexaDobl_Test_Evaluation;

  procedure HexaDobl_Test_Input_Output is

    maxdeg : constant integer32 := 16;
    s : HexaDobl_Complex_Series.Series(maxdeg);

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
  end HexaDobl_Test_Input_Output;

  procedure HexaDobl_Test_Polynomial_Series is

    maxdeg : constant integer32 := 16;
    sp : HexaDobl_CSeries_Polynomials.Poly;
    ps : HexaDobl_Polynomial_CSeries.Poly(maxdeg);
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
    ps := HexaDobl_Polynomial_CSeries.Create(sp);
    HexaDobl_CSeries_Polynomials.Clear(sp);
    sp := HexaDobl_Polynomial_CSeries.Create(ps);
    new_line;
    put_line("The polynomial from the created polynomial series :");
    Complex_Series_and_Polynomials_io.put(sp,1);
  end HexaDobl_Test_Polynomial_Series;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to test polynomials with series coefficients :");
    put_line("  0. test conversion from/to ordinary polynomials;");
    put_line("  1. test evaluation at power series;");
    put_line("  2. test symbolic input and output;");
    put_line("  3. test conversion to polynomial series.");
    put("Type 0, 1, 2, or 3 to select a test : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => HexaDobl_Test_Conversion;
      when '1' => HexaDobl_Test_Evaluation;
      when '2' => HexaDobl_Test_Input_Output;
      when '3' => HexaDobl_Test_Polynomial_Series;
      when others => null;
    end case;
  end Main;

end Test_HexaDobl_CSeries_Polynomials;
