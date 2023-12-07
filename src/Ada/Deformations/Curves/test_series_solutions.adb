with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_System_and_Solutions_io;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Systems_io;   use TripDobl_Complex_Poly_Systems_io;
with TripDobl_System_and_Solutions_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with QuadDobl_System_and_Solutions_io;
with Standard_Complex_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with TripDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with Complex_Series_and_Polynomials;
with Complex_Series_and_Polynomials_io;
with Series_and_Solutions;
with Power_Series_Methods;              use Power_Series_Methods;

package body Test_Series_Solutions is

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in Standard_CSeries_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Series_VecVecs.VecVec ) is

    maxdeg,nbrit : integer32 := 0;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    put("Give the maximal degree of the series : "); get(maxdeg); new_line;
    if echelon then
      put_line("Echelon Newton will be applied.");
      Run_Echelon_Newton(maxdeg,nbrit,p,s,true,true);
    else
      if nq = dim then
        put_line("LU Newton will be applied.");
        Run_LU_Newton(maxdeg,nbrit,p,s,true,true);
      else
        put_line("QR Newton will be applied.");
        Run_QR_Newton(maxdeg,nbrit,p,s,true,true);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in DoblDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Series_VecVecs.VecVec ) is

    maxdeg,nbrit : integer32 := 0;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    put("Give the maximal degree of the series : "); get(maxdeg); new_line;
    if echelon then
      put_line("Echelon Newton will be applied.");
      Run_Echelon_Newton(maxdeg,nbrit,p,s,true,true);
    else
      if nq = dim then
        put_line("LU Newton will be applied.");
        Run_LU_Newton(maxdeg,nbrit,p,s,true,true);
      else
        put_line("QR Newton will be applied.");
        Run_QR_Newton(maxdeg,nbrit,p,s,true,true);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in TripDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Series_VecVecs.VecVec ) is

    maxdeg,nbrit : integer32 := 0;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    put("Give the maximal degree of the series : "); get(maxdeg); new_line;
    if echelon then
      put_line("Echelon Newton will be applied.");
      Run_Echelon_Newton(maxdeg,nbrit,p,s,true,true);
    else
      if nq = dim then
        put_line("LU Newton will be applied.");
        Run_LU_Newton(maxdeg,nbrit,p,s,true,true);
      else
        put_line("QR Newton will be applied.");
        Run_QR_Newton(maxdeg,nbrit,p,s,true,true);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
             ( nq,dim : in integer32; echelon : in boolean;
               p : in QuadDobl_CSeries_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Series_VecVecs.VecVec ) is

    maxdeg,nbrit : integer32 := 0;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    put("Give the maximal degree of the series : ");
    get(maxdeg); new_line;
    if echelon then
      put_line("Echelon Newton will be applied.");
      Run_Echelon_Newton(maxdeg,nbrit,p,s,true,true);
    else
      if nq = dim then
        put_line("LU Newton will be applied.");
        Run_LU_Newton(maxdeg,nbrit,p,s,true,true);
      else
        put_line("QR Newton will be applied.");
        Run_QR_Newton(maxdeg,nbrit,p,s,true,true);
      end if;
    end if;
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant Standard_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,dim,echelon,srp,srv);
    Standard_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               s : in DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant DoblDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : DoblDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,dim,echelon,srp,srv);
    DoblDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in TripDobl_Complex_Poly_Systems.Poly_Sys;
               s : in TripDobl_Complex_Solutions.Solution_List ) is

    use TripDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant TripDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : TripDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,dim,echelon,srp,srv);
    TripDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32; echelon : in boolean;
               p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               s : in QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant QuadDobl_Complex_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : QuadDobl_CSeries_Poly_Systems.Poly_Sys(p'range)
        := Complex_Series_and_Polynomials.System_to_Series_System(p,idx);

  begin
    Run_Newton(nq,dim,echelon,srp,srv);
    QuadDobl_CSeries_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Standard_Test_at_Zero_Order ( echelon : in boolean ) is

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

  procedure TripDobl_Test_at_Zero_Order ( echelon : in boolean ) is

    use TripDobl_Complex_Polynomials;
    use TripDobl_Complex_Poly_Systems;
    use TripDobl_Complex_Solutions;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the file name for a system and solutions ...");
    TripDobl_System_and_Solutions_io.get(lp,sols);
    new_line;
    nq := lp'last;
    nv := integer32(Number_of_Unknowns(lp(lp'first)));
    put("Read a system of "); put(nq,1);
    put(" equations in "); put(nv,1); put_line(" unknowns.");
    put("Read "); put(integer32(Length_Of(sols)),1); put_line(" solutions.");
    new_line;
    put("Give the index of the parameter : "); get(idx);
    Run_Newton(nq,idx,echelon,lp.all,sols);
  end TripDobl_Test_at_Zero_Order;

  procedure QuadDobl_Test_at_Zero_Order ( echelon : in boolean ) is

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
      srp : constant Standard_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(nq,nv,echelon,srp,s);
    end;
  end Standard_Test_at_Series;

  procedure DoblDobl_Test_at_Series ( echelon : in boolean ) is

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
      srp : constant DoblDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(nq,nv,echelon,srp,s);
    end;
  end DoblDobl_Test_at_Series;

  procedure TripDobl_Test_at_Series ( echelon : in boolean ) is

    use TripDobl_Complex_Polynomials;
    use TripDobl_Complex_Poly_Systems;

    lp : Link_to_Poly_Sys;
    nq,nv,idx : integer32 := 0;
    srv : TripDobl_Complex_Series_Vectors.Link_to_Vector;

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
      s : TripDobl_Complex_Series_VecVecs.VecVec(1..1);
      srp : constant TripDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(nq,nv,echelon,srp,s);
    end;
  end TripDobl_Test_at_Series;

  procedure QuadDobl_Test_at_Series ( echelon : in boolean ) is

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
      srp : constant QuadDobl_CSeries_Poly_Systems.Poly_Sys(lp'range)
          := Complex_Series_and_Polynomials.System_to_Series_System(lp.all,idx);
    begin
      s(1) := srv;
      Run_Newton(nq,nv,echelon,srp,s);
    end;
  end QuadDobl_Test_at_Series;

  procedure Main is

    ans,prc : character;
    echelon : boolean;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. triple double precision");
    put_line("  3. quad double precision");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(prc,"0123");
    new_line;
    put("Use the lower triangular echelon form ? (y/n) ");
    Ask_Yes_or_No(ans); echelon := (ans = 'y');
    new_line;
    put("Start Newton's method at zero order ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      case prc is
        when '0' => Standard_Test_at_Zero_Order(echelon);
        when '1' => DoblDobl_Test_at_Zero_Order(echelon);
        when '2' => TripDobl_Test_at_Zero_Order(echelon);
        when '3' => QuadDobl_Test_at_Zero_Order(echelon);
        when others => null;
      end case;
    else
      case prc is
        when '0' => Standard_Test_at_Series(echelon);
        when '1' => DoblDobl_Test_at_Series(echelon);
        when '2' => TripDobl_Test_at_Series(echelon);
        when '3' => QuadDobl_Test_at_Series(echelon);
        when others => null;
      end case;
    end if;
  end Main;

end Test_Series_Solutions;
