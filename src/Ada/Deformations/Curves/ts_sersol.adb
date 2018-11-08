with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_System_and_Solutions_io;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with Standard_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_Systems;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with Series_and_Solutions;
with Power_Series_Methods;              use Power_Series_Methods;

procedure ts_sersol is

-- DESCRIPTION :
--   Test on the development of series as solutions of polynomial systems.

  procedure Run_Newton
             ( nq,idx,dim : in integer32; echelon : in boolean;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Applies Newton's method in standard double precision,
  --   at the polynomial system p, starting at the series in s.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   dim     number of coordinates in the series;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a sequence of series to start Newton's method at.

    len : constant integer32 := s'last;
    nbrit : integer32 := 0;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    if echelon then
      put_line("Echelon Newton will be applied.");
      Run_Echelon_Newton(nbrit,p,s,true,true);
    else
      if nq = dim then
        put_line("LU Newton will be applied.");
        Run_LU_Newton(nbrit,p,s,true,true);
      else
        put_line("QR Newton will be applied.");
        Run_QR_Newton(nbrit,p,s,true,true);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx,dim : in integer32; echelon : in boolean;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Applies Newton's method in double double precision,
  --   at the polynomial system p, starting at the series in s.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   dim     number of coordinates in the series;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a sequence of series to start Newton's method at.

    len : constant integer32 := s'last;
    nbrit : integer32 := 0;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    if echelon then
      put_line("Echelon Newton will be applied.");
      Run_Echelon_Newton(nbrit,p,s,true,true);
    else
      if nq = dim then
        put_line("LU Newton will be applied.");
        Run_LU_Newton(nbrit,p,s,true,true);
      else
        put_line("QR Newton will be applied.");
        Run_QR_Newton(nbrit,p,s,true,true);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx,dim : in integer32; echelon : in boolean;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Series_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Applies Newton's method in quad double precision,
  --   at the polynomial system p, starting at the series in s.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   dim     number of coordinates in the series;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a sequence of series to start Newton's method at.

    len : constant integer32 := s'last;
    nbrit : integer32 := 0;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    if echelon then
      put_line("Echelon Newton will be applied.");
      Run_Echelon_Newton(nbrit,p,s,true,true);
    else
      if nq = dim then
        put_line("LU Newton will be applied.");
        Run_LU_Newton(nbrit,p,s,true,true);
      else
        put_line("QR Newton will be applied.");
        Run_QR_Newton(nbrit,p,s,true,true);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in standard double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant Standard_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,idx,dim,echelon,srp,srv);
    Standard_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in double double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant DoblDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,idx,dim,echelon,srp,srv);
    DoblDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.
  --   Newton's method is performed in quad double precision.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   echelon is the flag for the echelon Newton's method to be used;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant QuadDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,idx,dim,echelon,srp,srv);
    QuadDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Standard_Test_at_Zero_Order ( echelon : in boolean ) is

  -- DESCRIPTION :
  --   Performs a test in standard double precision,
  --   prompting the user for a system and its solutions,
  --   starting Newton's method at zero-th order series.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the file name for a system and solutions ...");
    Standard_System_and_Solutions_io.get(lp,sols);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton(nq,idx,echelon,lp.all,sols);
  end Standard_Test_at_Zero_Order;

  procedure DoblDobl_Test_at_Zero_Order ( echelon : in boolean ) is

  -- DESCRIPTION :
  --   Performs a test in double double precision,
  --   prompting the user for a system and its solutions,
  --   starting Newton's method at zero-th order series.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the file name for a system and solutions ...");
    DoblDobl_System_and_Solutions_io.get(lp,sols);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton(nq,idx,echelon,lp.all,sols);
  end DoblDobl_Test_at_Zero_Order;

  procedure QuadDobl_Test_at_Zero_Order ( echelon : in boolean ) is

  -- DESCRIPTION :
  --   Performs a test in quad double precision,
  --   prompting the user for a system and its solutions,
  --   starting Newton's method at zero-th order series.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the file name for a system and solutions ...");
    QuadDobl_System_and_Solutions_io.get(lp,sols);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton(nq,idx,echelon,lp.all,sols);
  end QuadDobl_Test_at_Zero_Order;

  procedure Standard_Test_at_Series ( echelon : in boolean ) is

  -- DESCRIPTION :
  --   Performs a test in standard double precision,
  --   prompting the user for a system and its solutions,
  --   starting Newton's method at a series.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : Standard_Complex_Series_Vectors.Link_to_Vector;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    declare
      s : Standard_Complex_Series_VecVecs.VecVec(1..1);
      srp : Standard_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(nq,idx,nv,echelon,srp,s);
    end;
  end Standard_Test_at_Series;

  procedure DoblDobl_Test_at_Series ( echelon : in boolean ) is

  -- DESCRIPTION :
  --   Performs a test in double double precision,
  --   prompting the user for a system and its solutions,
  --   starting Newton's method at a series.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : DoblDobl_Complex_Series_Vectors.Link_to_Vector;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    declare
      s : DoblDobl_Complex_Series_VecVecs.VecVec(1..1);
      srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(nq,idx,nv,echelon,srp,s);
    end;
  end DoblDobl_Test_at_Series;

  procedure QuadDobl_Test_at_Series ( echelon : in boolean ) is

  -- DESCRIPTION :
  --   Performs a test in quad double precision,
  --   prompting the user for a system and its solutions,
  --   starting Newton's method at a series.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : QuadDobl_Complex_Series_Vectors.Link_to_Vector;

  begin
    new_line;
    put_line("Reading a polynomial system ...");
    get(lp);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    new_line;
    put_line("Reading a series to start Newton's method at ...");
    Complex_Series_and_Polynomials_io.get(srv);
    declare
      s : QuadDobl_Complex_Series_VecVecs.VecVec(1..1);
      srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(nq,idx,nv,echelon,srp,s);
    end;
  end QuadDobl_Test_at_Series;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user to select the working precision
  --   and then launches the corresponding test.

    prc,ans : character;
    echelon : boolean;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the working precision : ");
    Ask_Alternative(prc,"012");
    new_line;
    put("Use the lower triangular echelon form ? (y/n) ");
    Ask_Yes_or_No(ans);
    echelon := (ans = 'y');
    new_line;
    put("Start Newton's method at zero order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      case prc is
        when '0' => Standard_Test_at_Zero_Order(echelon);
        when '1' => DoblDobl_Test_at_Zero_Order(echelon);
        when '2' => QuadDobl_Test_at_Zero_Order(echelon);
        when others => null;
      end case;
    else
      case prc is
        when '0' => Standard_Test_at_Series(echelon);
        when '1' => DoblDobl_Test_at_Series(echelon);
        when '2' => QuadDobl_Test_at_Series(echelon);
        when others => null;
      end case;
    end if;
  end Main;

begin
  Main;
end ts_sersol;
