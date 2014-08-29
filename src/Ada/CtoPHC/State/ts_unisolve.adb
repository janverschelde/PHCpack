with text_io;                           use text_io;
with Interfaces.C;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Random_Polynomials;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Random_Polynomials;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Random_Polynomials;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Random_Polynomials;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;     use Multprec_Complex_Solutions_io;
with Standard_PolySys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_Solutions_Container;
with Multprec_PolySys_Container;
with Multprec_Solutions_Container;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with unisolve;

procedure ts_unisolve is

-- DESCRIPTION :
--   Test on computing the roots of a polynomial in one variable.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a choice between a random or a user given
  --   polynomials, asks for the number of iterations and then makes
  --   the call to unisolve, in standard double precision.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    d : natural32 := 0;
    p : Poly;
    s : Poly_Sys(1..1);
    ar,br : C_Integer_Array(0..0);
    ap,bp : C_IntArrs.Pointer;
    cr : C_Double_Array(0..0);
    eps : constant double_float := 1.0E-12;
    cp : C_DblArrs.Pointer;
    r : integer32;
    max : natural32 := 0;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put("Random or given polynomial ? (r/g) ");
    Ask_Alternative(ans,"rg");
    if ans = 'g' then
      Symbol_Table.Init(1);
      put("Give a polynomial : "); get(p);
    else
      new_line;
      put("Give the degree : "); get(d);
      p := Standard_Random_Polynomials.Random_Dense_Poly(1,d,0);
    end if;
    put_line("the polynomial :"); put(p); new_line;
    new_line;
    put("Give the maximum number of iterations : "); get(max);
    ar(0) := Interfaces.C.int(max);
    br(0) := Interfaces.C.int(0);
    cr(0) := Interfaces.C.double(eps);
    ap := ar(0)'unchecked_access;
    bp := br(0)'unchecked_access;
    cp := cr(0)'unchecked_access;
    s(1) := p;
    Standard_PolySys_Container.Initialize(s);
    r := unisolve(1,ap,bp,cp);
    sols := Standard_Solutions_Container.Retrieve;
    if not Is_Null(sols)
     then put(Length_Of(sols),1,sols);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a choice between a random or a user given
  --   polynomials, asks for the number of iterations and then makes
  --   the call to unisolve, in double double precision.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    d : natural32 := 0;
    p : Poly;
    s : Poly_Sys(1..1);
    ar,br : C_Integer_Array(0..0);
    ap,bp : C_IntArrs.Pointer;
    cr : C_Double_Array(0..0);
    eps : constant double_float := 1.0E-28;
    cp : C_DblArrs.Pointer;
    r : integer32;
    max : natural32 := 0;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put("Random or given polynomial ? (r/g) ");
    Ask_Alternative(ans,"rg");
    if ans = 'g' then
      Symbol_Table.Init(1);
      new_line;
      put("Give a polynomial : "); get(p);
    else
      new_line;
      put("Give the degree : "); get(d);
      p := DoblDobl_Random_Polynomials.Random_Dense_Poly(1,d,0);
    end if;
    put_line("the polynomial :"); put(p);
    new_line;
    put("Give the maximum number of iterations : "); get(max);
    ar(0) := Interfaces.C.int(max);
    br(0) := Interfaces.C.int(0);
    cr(0) := Interfaces.C.double(eps);
    ap := ar(0)'unchecked_access;
    bp := br(0)'unchecked_access;
    cp := cr(0)'unchecked_access;
    s(1) := p;
    DoblDobl_PolySys_Container.Initialize(s);
    r := unisolve(2,ap,bp,cp);
    sols := DoblDobl_Solutions_Container.Retrieve;
    if not Is_Null(sols)
     then put(Length_Of(sols),1,sols);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a choice between a random or a user given
  --   polynomials, asks for the number of iterations and then makes
  --   the call to unisolve, in quad double precision.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    d : natural32 := 0;
    p : Poly;
    s : Poly_Sys(1..1);
    ar,br : C_Integer_Array(0..0);
    ap,bp : C_IntArrs.Pointer;
    cr : C_Double_Array(0..0);
    eps : constant double_float := 1.0E-28;
    cp : C_DblArrs.Pointer;
    r : integer32;
    max : natural32 := 0;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put("Random or given polynomial ? (r/g) ");
    Ask_Alternative(ans,"rg");
    if ans = 'g' then
      Symbol_Table.Init(1);
      new_line;
      put("Give a polynomial : "); get(p);
    else
      new_line;
      put("Give the degree : "); get(d);
      p := QuadDobl_Random_Polynomials.Random_Dense_Poly(1,d,0);
    end if;
    put_line("the polynomial :"); put(p);
    new_line;
    put("Give the maximum number of iterations : "); get(max);
    ar(0) := Interfaces.C.int(max);
    br(0) := Interfaces.C.int(0);
    cr(0) := Interfaces.C.double(eps);
    ap := ar(0)'unchecked_access;
    bp := br(0)'unchecked_access;
    cp := cr(0)'unchecked_access;
    s(1) := p;
    QuadDobl_PolySys_Container.Initialize(s);
    r := unisolve(3,ap,bp,cp);
    sols := QuadDobl_Solutions_Container.Retrieve;
    if not Is_Null(sols)
     then put(Length_Of(sols),1,sols);
    end if;
  end QuadDobl_Test;

  procedure Multprec_Test is

  -- DESCRIPTION :
  --   Prompts the user for a choice between a random or a user given
  --   polynomials, asks for the number of iterations and then makes
  --   the call to unisolve, in arbitrary multiprecision.

    use Multprec_Floating_Numbers;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Solutions;

    d : natural32 := 0;
    p : Poly;
    s : Poly_Sys(1..1);
    ar : C_Integer_Array(0..1);
    br : C_Integer_Array(0..0);
    ap,bp : C_IntArrs.Pointer;
    cr : C_Double_Array(0..0);
    eps : double_float;
    cp : C_DblArrs.Pointer;
    r : integer32;
    deci,max : natural32 := 0;
    size : natural32;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put("Give the number of decimal places in the working precision : ");
    get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Multprec_Complex_Polynomials_io.Set_Working_Precision(size);
    eps := 10.0**(-natural((deci*3)/4));
    new_line;
    put("Random or given polynomial ? (r/g) ");
    Ask_Alternative(ans,"rg");
    if ans = 'g' then
      Symbol_Table.Init(1);
      new_line;
      put("Give a polynomial : "); get(p);
    else
      new_line;
      put("Give the degree : "); get(d);
      p := Multprec_Random_Polynomials.Random_Dense_Poly(1,d,0);
    end if;
    put_line("the polynomial :"); put(p);
    new_line;
    put("Give the maximum number of iterations : "); get(max);
    ar(0) := Interfaces.C.int(deci);
    ar(1) := Interfaces.C.int(max);
    br(0) := Interfaces.C.int(0);
    cr(0) := Interfaces.C.double(eps);
    ap := ar(0)'unchecked_access;
    bp := br(0)'unchecked_access;
    cp := cr(0)'unchecked_access;
    s(1) := p;
    Multprec_PolySys_Container.Initialize(s);
    r := unisolve(4,ap,bp,cp);
    sols := Multprec_Solutions_Container.Retrieve;
    if not Is_Null(sols)
     then put(Length_Of(sols),1,sols);
    end if;
  end Multprec_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision level
  --   and calls the corresponding test procedure.

    ans : character;

  begin
    new_line;
    put_line("Calling the gateway to the univariate root finders ...");
    new_line;
    put_line("MENU to set the precision level :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put_line("  3. arbitrary multiprecision.");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when '3' => Multprec_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_unisolve;
