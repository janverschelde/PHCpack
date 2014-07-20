with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Standard_PolySys_Container;
with DoblDobl_PolySys_Container;
with QuadDobl_PolySys_Container;
with Multprec_PolySys_Container;
with Laurent_Systems_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_LaurSys_Container;

procedure ts_syscon is

-- DESCRIPTION :
--   This procedure is an interactive testing facility on the
--   operations in the package systems container.

  procedure Standard_Test_Retrievals is

  -- DESCRIPTION :
  --   Tests the selectors in the package.

    use Standard_Complex_Polynomials;

    n : constant natural32 := Standard_PolySys_Container.Dimension;
    m : natural32;

  begin
    put("Dimension of the system : "); put(n,1); new_line;
    for i in 1..integer32(n) loop
      m := Standard_PolySys_Container.Number_of_Terms(i);
      put("Number of terms in polynomial "); put(i,1);
      put(" : "); put(m,1); new_line;
      declare
        p : Poly := Null_Poly;
      begin
        for j in 1..m loop
          declare
            t : Term := Standard_PolySys_Container.Retrieve_Term(i,j);
          begin
            Add(p,t);
            Clear(t);
          end;
        end loop;
        put(p); new_line;
        put(Standard_PolySys_Container.Retrieve(i)); new_line;
      end;
    end loop;
  end Standard_Test_Retrievals;

  procedure DoblDobl_Test_Retrievals is

  -- DESCRIPTION :
  --   Tests the selectors in the package.

    use DoblDobl_Complex_Polynomials;

    n : constant natural32 := DoblDobl_PolySys_Container.Dimension;
    m : natural32;

  begin
    put("Dimension of the system : "); put(n,1); new_line;
    for i in 1..integer32(n) loop
      m := DoblDobl_PolySys_Container.Number_of_Terms(i);
      put("Number of terms in polynomial "); put(i,1);
      put(" : "); put(m,1); new_line;
      declare
        p : Poly := Null_Poly;
      begin
        for j in 1..m loop
          declare
            t : Term := DoblDobl_PolySys_Container.Retrieve_Term(i,j);
          begin
            Add(p,t);
            Clear(t);
          end;
        end loop;
        put(p); new_line;
        put(DoblDobl_PolySys_Container.Retrieve(i)); new_line;
      end;
    end loop;
  end DoblDobl_Test_Retrievals;

  procedure QuadDobl_Test_Retrievals is

  -- DESCRIPTION :
  --   Tests the selectors in the package.

    use QuadDobl_Complex_Polynomials;

    n : constant natural32 := QuadDobl_PolySys_Container.Dimension;
    m : natural32;

  begin
    put("Dimension of the system : "); put(n,1); new_line;
    for i in 1..integer32(n) loop
      m := QuadDobl_PolySys_Container.Number_of_Terms(i);
      put("Number of terms in polynomial "); put(i,1);
      put(" : "); put(m,1); new_line;
      declare
        p : Poly := Null_Poly;
      begin
        for j in 1..m loop
          declare
            t : Term := QuadDobl_PolySys_Container.Retrieve_Term(i,j);
          begin
            Add(p,t);
            Clear(t);
          end;
        end loop;
        put(p); new_line;
        put(QuadDobl_PolySys_Container.Retrieve(i)); new_line;
      end;
    end loop;
  end QuadDobl_Test_Retrievals;

  procedure Multprec_Test_Retrievals is

  -- DESCRIPTION :
  --   Tests the selectors in the package.

    use Multprec_Complex_Polynomials;

    n : constant natural32 := Multprec_PolySys_Container.Dimension;
    m : natural32;

  begin
    put("Dimension of the system : "); put(n,1); new_line;
    for i in 1..integer32(n) loop
      m := Multprec_PolySys_Container.Number_of_Terms(i);
      put("Number of terms in polynomial "); put(i,1);
      put(" : "); put(m,1); new_line;
      declare
        p : Poly := Null_Poly;
      begin
        for j in 1..m loop
          declare
            t : Term := Multprec_PolySys_Container.Retrieve_Term(i,j);
          begin
            Add(p,t);
            Clear(t);
          end;
        end loop;
        put(p); new_line;
        put(Multprec_PolySys_Container.Retrieve(i)); new_line;
      end;
    end loop;
  end Multprec_Test_Retrievals;

  procedure Standard_Test_Laurent_Retrievals is

  -- DESCRIPTION :
  --   Test on retrieving data from the Laurent system container,
  --   for coefficients in standard double precision.

    n : constant natural32 := Laurent_Systems_Container.Dimension;
    lp : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    put("Dimension of the system : "); put(n,1); new_line;
    lp := Laurent_Systems_Container.Retrieve;
    put_line("The Laurent polynomial system : "); put(lp.all);
  end Standard_Test_Laurent_Retrievals;

  procedure DoblDobl_Test_Laurent_Retrievals is

  -- DESCRIPTION :
  --   Test on retrieving data from the Laurent system container,
  --   for coefficients in double double precision.

    n : constant natural32 := DoblDobl_LaurSys_Container.Dimension;
    lp : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    put("Dimension of the system : "); put(n,1); new_line;
    lp := DoblDobl_LaurSys_Container.Retrieve;
    put_line("The Laurent polynomial system : "); put(lp.all);
  end DoblDobl_Test_Laurent_Retrievals;

  procedure QuadDobl_Test_Laurent_Retrievals is

  -- DESCRIPTION :
  --   Test on retrieving data from the Laurent system container,
  --   for coefficients in double double precision.

    n : constant natural32 := QuadDobl_LaurSys_Container.Dimension;
    lp : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;

  begin
    put("Dimension of the system : "); put(n,1); new_line;
    lp := QuadDobl_LaurSys_Container.Retrieve;
    put_line("The Laurent polynomial system : "); put(lp.all);
  end QuadDobl_Test_Laurent_Retrievals;

  procedure Standard_Test_Additions
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the constructors of the package.

    use Standard_Complex_Polynomials;

    ind : integer32;

    procedure Add_Term ( t : in Term; continue : out boolean ) is
    begin
      Standard_PolySys_Container.Add_Term(ind,t);
      continue := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    Standard_PolySys_Container.Initialize(p'last);
    for i in p'range loop
      ind := i;
      Add_Terms(p(i));
    end loop;
    put_line("The polynomial system : ");
    put(Standard_PolySys_Container.Retrieve.all);
  end Standard_Test_Additions;

  procedure DoblDobl_Test_Additions
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the constructors of the package.

    use DoblDobl_Complex_Polynomials;

    ind : integer32;

    procedure Add_Term ( t : in Term; continue : out boolean ) is
    begin
      DoblDobl_PolySys_Container.Add_Term(ind,t);
      continue := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    DoblDobl_PolySys_Container.Initialize(p'last);
    for i in p'range loop
      ind := i;
      Add_Terms(p(i));
    end loop;
    put_line("The polynomial system : ");
    put(DoblDobl_PolySys_Container.Retrieve.all);
  end DoblDobl_Test_Additions;

  procedure QuadDobl_Test_Additions
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the constructors of the package.

    use QuadDobl_Complex_Polynomials;

    ind : integer32;

    procedure Add_Term ( t : in Term; continue : out boolean ) is
    begin
      QuadDobl_PolySys_Container.Add_Term(ind,t);
      continue := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    QuadDobl_PolySys_Container.Initialize(p'last);
    for i in p'range loop
      ind := i;
      Add_Terms(p(i));
    end loop;
    put_line("The polynomial system : ");
    put(QuadDobl_PolySys_Container.Retrieve.all);
  end QuadDobl_Test_Additions;

  procedure Multprec_Test_Additions
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the constructors of the package.

    use Multprec_Complex_Polynomials;

    ind : integer32;

    procedure Add_Term ( t : in Term; continue : out boolean ) is
    begin
      Multprec_PolySys_Container.Add_Term(ind,t);
      continue := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    Multprec_PolySys_Container.Initialize(p'last);
    for i in p'range loop
      ind := i;
      Add_Terms(p(i));
    end loop;
    put_line("The polynomial system : ");
    put(Multprec_PolySys_Container.Retrieve.all);
  end Multprec_Test_Additions;

  procedure Main is

    st_lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    dd_lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    qd_lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp_lp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    st_lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    dd_lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    qd_lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    ans : character;

  begin
    new_line;
    put_line("MENU to test the operations in the systems container :");
    put_line("  1. test with standard double complex numbers;");
    put_line("  2. test with double double complex numbers;");
    put_line("  3. test with quad double complex numbers;");
    put_line("  4. test with arbitrary precision complex numbers;");
    put_line("  5. test Laurent systems in standard double precision;");
    put_line("  6. test Laurent systems in double double precision;");
    put_line("  7. test Laurent systems in quad double precision;");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to select the precision : ");
    Ask_Alternative(ans,"1234567");
    new_line;
    case ans is
      when '1' =>
        get(st_lp);
        Standard_PolySys_Container.Initialize(st_lp.all);
        Standard_Test_Retrievals;
        Standard_PolySys_Container.Clear;
        Standard_Test_Additions(st_lp.all);
      when '2' =>
        get(dd_lp);
        DoblDobl_PolySys_Container.Initialize(dd_lp.all);
        DoblDobl_Test_Retrievals;
        DoblDobl_PolySys_Container.Clear;
        DoblDobl_Test_Additions(dd_lp.all);
      when '3' =>
        get(qd_lp);
        QuadDobl_PolySys_Container.Initialize(qd_lp.all);
        QuadDobl_Test_Retrievals;
        QuadDobl_PolySys_Container.Clear;
        QuadDobl_Test_Additions(qd_lp.all);
      when '4' =>
        get(mp_lp);
        Multprec_PolySys_Container.Initialize(mp_lp.all);
        Multprec_Test_Retrievals;
        Multprec_PolySys_Container.Clear;
        Multprec_Test_Additions(mp_lp.all);
      when '5' =>
        get(st_lq);
        Laurent_Systems_Container.Initialize(st_lq.all);
        Standard_Test_Laurent_Retrievals;
      when '6' =>
        get(dd_lq);
        DoblDobl_LaurSys_Container.Initialize(dd_lq.all);
        DoblDobl_Test_Laurent_Retrievals;
      when '7' =>
        get(qd_lq);
        QuadDobl_LaurSys_Container.Initialize(qd_lq.all);
        QuadDobl_Test_Laurent_Retrievals;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_syscon;
