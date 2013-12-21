with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Lists_of_Floating_Vectors_io;       use Lists_of_Floating_Vectors_io;
with Arrays_of_Floating_Vector_Lists;    use Arrays_of_Floating_Vector_Lists;
with Symbol_Table;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
with Standard_Tableau_Formats;           use Standard_Tableau_Formats;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;

procedure ts_supports is

-- DESCRIPTION :
--   This procedure allows interactive testing of supports
--   of a polynomial system.

  procedure Select_Supported_Subpolynomial ( p : in Poly ) is

    n : constant natural32 := Number_of_Unknowns(p);
    k : natural32 := 0;
    s : List;
    sub_p : Poly;

  begin
    new_line;
    put("Give #selected monomials from polynomial : ");
    get(k);
    put("Give "); put(k,1); put(" vectors of length ");
    put(n,1); put_line(" :"); get(n,k,s);
    put_line("The selected support : "); put(s);
    sub_p := Select_Terms(p,s);
    put_line("The selected polynomial : "); put(sub_p); new_line;
    declare
      dim : constant natural32
          := Standard_Complex_Polynomials.Number_of_Unknowns(p);
      deg : Standard_Complex_Polynomials.Degrees
          := new Standard_Natural_Vectors.Vector(1..integer32(dim));
      len : constant natural32 := Length_Of(s);
      cff : Standard_Complex_Vectors.Vector(1..integer32(len));
    begin
      Select_Coefficients(p,s,dim,deg,cff);
      put_line("The coefficients : "); put_line(cff);
      Standard_Complex_Polynomials.Clear(deg);
    end;
  end Select_Supported_Subpolynomial;

  procedure Select_Supported_Subsystem
              ( p : in Poly_Sys; r : in integer32 ) is

    n : constant natural32 := natural32(p'last);
    k : natural32 := 0;
    m : Standard_Integer_Vectors.Vector(1..r);
    s : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    sub_p : Poly_Sys(p'range);

  begin
    put("Give the type of mixture : "); get(m);
    for i in 1..r loop
      new_line;
      put("Give #selected monomials from support "); put(i,1);
      put(" : "); get(k);
      put("Give "); put(k,1); put(" vectors of length ");
      put(n,1); put_line(" :"); get(n,k,s(i));
      put_line("The selected support : "); put(s(i));
    end loop;
    sub_p := Select_Terms(p,m,s);
    put_line("The selected subsystem is "); put(sub_p); 
    declare
      dim : constant natural32
          := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      deg : Standard_Complex_Polynomials.Degrees
          := new Standard_Natural_Vectors.Vector(1..integer32(dim));
      cff : Standard_Complex_VecVecs.VecVec(p'range);
      ind : integer32 := 0;
    begin
      for i in m'range loop
        for j in 1..m(i) loop
          ind := ind + 1;
          cff(ind) := new Standard_Complex_Vectors.Vector
                            (1..integer32(Length_Of(s(i))));
        end loop;
      end loop;
      Select_Coefficients(p,m,s,dim,deg,cff);
      put_line("The coefficients : "); put_line(cff);
      Standard_Complex_Polynomials.Clear(deg);
    end;
  end Select_Supported_Subsystem;

  procedure Select_Supported_Subsystem ( p : in Poly_Sys ) is

    n : constant natural32 := natural32(p'last);
    k : natural32 := 0;
    s : Arrays_of_Floating_Vector_Lists.Array_of_Lists(p'range);
    sub_p : Poly_Sys(p'range);

  begin
    for i in p'range loop
      new_line;
      put("Give #selected monomials from polynomial "); put(i,1);
      put(" : "); get(k);
      put("Give "); put(k,1); put(" vectors of length ");
      put(n,1); put_line(" :"); get(n,k,s(i));
      put_line("The selected support : "); put(s(i));
    end loop;
    sub_p := Select_Terms(p,s);
    put_line("The selected subsystem is "); put(sub_p); 
    declare
      dim : constant natural32
          := Standard_Complex_Polynomials.Number_of_Unknowns(p(p'first));
      deg : Standard_Complex_Polynomials.Degrees
          := new Standard_Natural_Vectors.Vector(1..integer32(dim));
      cff : Standard_Complex_VecVecs.VecVec(p'range);
    begin
      for i in cff'range loop
        cff(i) := new Standard_Complex_Vectors.Vector
                        (1..integer32(Length_Of(s(i))));
      end loop;
      Select_Coefficients(p,s,dim,deg,cff);
      put_line("The coefficients : "); put_line(cff);
      Standard_Complex_Polynomials.Clear(deg);
    end;
  end Select_Supported_Subsystem;

  procedure Compute_Tableau ( p : in Poly_Sys ) is

    m : integer32;

  begin
    put(p'last,1); new_line;
    for i in p'range loop
      m := integer32(Standard_Complex_Polynomials.Number_of_Terms(p(i)));
      declare
        cff : Standard_Complex_Vectors.Vector(1..m);
        exp : Standard_Natural_VecVecs.VecVec(1..m);
      begin
        Extract_Coefficients_and_Exponents(p(i),cff,exp);
        put(m,1); new_line;
        for j in 1..m loop
          put(cff(j));
          put("  ");
          put(exp(j));
          new_line;
        end loop;
      end;
    end loop;
  end Compute_Tableau;

  procedure Main is

    ans : character;
    lp : Link_to_Poly_Sys;
    n,r : natural32 := 0;
    p : Poly;

  begin
    new_line;
    put_line("Testing operations on supports of polynomial systems...");
    new_line;
    put_line("Choose one of the following operations : ");
    put_line("  0. Show the supports of a polynomial system;");
    put_line("  1. Select terms from a given polynomial;");
    put_line("  2. Select a supported subsystem of a fully-mixed system.");
    put_line("  3. Select a supported subsystem of a semi-mixed system.");
    put_line("  4. Compute tableau format of a polynomial system.");
    put("Type 1, 2, 3, or 4 to make your choice : ");
    Ask_Alternative(ans,"01234");
    new_line;
    case ans is 
      when '1' =>
        put("Give the number of unknowns : "); get(n);
        Symbol_Table.Init(n);
        put("Give a polynomial in "); put(n,1); put(" variables,");
        put_line(" terminate with semi-colon:");
        get(p);
        Select_Supported_Subpolynomial(p);
      when '2' =>
        get(lp);
        Select_Supported_Subsystem(lp.all);
      when '3' =>
        get(lp);
        new_line;
        put("Give the number of different supports : "); get(r);
        Select_Supported_Subsystem(lp.all,integer32(r));
      when '4' =>
        get(lp);
        Compute_Tableau(lp.all);
      when others =>
        get(lp);
        declare
          supports
            : constant Arrays_of_Integer_Vector_Lists.Array_of_Lists(lp'range)
            := Create(lp.all);
        begin
          put_line("The supports of the system : "); put(supports);
        end;
    end case;
  end Main;

begin
  Main;
end ts_supports;
