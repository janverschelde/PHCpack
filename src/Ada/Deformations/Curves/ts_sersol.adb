with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_VecVecs;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with Series_and_Polynomials;
with Series_and_Polynomials_io;
with Series_and_Solutions;
with Standard_Newton_Series;

procedure ts_sersol is

-- DESCRIPTION :
--   Test on the development of series as solutions of polynomial systems.

  procedure Run_LU_Newton
             ( nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               verbose : in boolean := false ) is

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p. 

    use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-20;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);

  begin
    if verbose
     then LU_Newton_Steps(standard_output,p,order,nbrit,s,info);
     else LU_Newton_Steps(p,order,nbrit,s,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(s,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(s);
      eva := Standard_Series_Poly_SysFun.Eval(p,s);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Run_LU_Newton;

  procedure Run_QR_Newton
             ( nbrit : in integer32;
               p : in Standard_Series_Poly_Systems.Poly_Sys;
               s : in out Standard_Dense_Series_Vectors.Vector;
               verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Applies as many steps with Newton's method as the value of nbrit,
  --   starting at the solution in s to the system p. 

    use Standard_Newton_Series;

    order : integer32 := 1;
    info : integer32;
    tol : constant double_float := 1.0E-12;
    eva : Standard_Dense_Series_Vectors.Vector(p'range);

  begin
    if verbose
     then QR_Newton_Steps(standard_output,p,order,nbrit,s,info);
     else QR_Newton_Steps(p,order,nbrit,s,info);
    end if;
    if info /= 0 then
      put("info = "); put(info,1); new_line;
    else
      Series_and_Polynomials.Filter(s,tol);
      put_line("The updated power series solution :");
      Series_and_Polynomials_io.put(s);
      eva := Standard_Series_Poly_SysFun.Eval(p,s);
      Series_and_Polynomials.Filter(eva,tol);
      put_line("The evaluated solution :");
      Series_and_Polynomials_io.put(eva);
    end if;
  end Run_QR_Newton;

  procedure Run_Newton
             ( nq,idx : in integer32;
               p : in Standard_Complex_Poly_Systems.Poly_Sys;
               s : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   The coordinates of the solution vectors in s are the leading
  --   coefficients in a power series solution to p.

  -- ON ENTRY :
  --   nq      number of equations in p;
  --   idx     index to the series parameter;
  --   p       a polynomial of nq equations in nv unknowns;
  --   s       a list of solutions.

    use Standard_Complex_Solutions;

    len : constant integer32 := integer32(Length_Of(s));
    dim : constant integer32 := Head_Of(s).n - 1;
    srv : constant Standard_Dense_Series_VecVecs.VecVec(1..len)
        := Series_and_Solutions.Create(s,idx);
    srp : Standard_Series_Poly_Systems.Poly_Sys(p'range)
        := Series_and_Polynomials.System_to_Series_System(p,idx);
    nbrit : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Number of coordinates in the series : "); put(dim,1); put_line(".");
    new_line;
    put("Give the number of steps in Newton's method : ");
    get(nbrit); new_line;
    if nq = dim then
      put_line("LU Newton will be applied.");
      for i in srv'range loop
        put("Running on solution "); put(i,1); put_line(" ...");
        Run_LU_Newton(nbrit,srp,srv(i).all);
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end loop;
    else
      put_line("QR Newton will be applied.");
      for i in srv'range loop
        put("Running on solution "); put(i,1); put_line(" ...");
        Run_QR_Newton(nbrit,srp,srv(i).all);
        put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end loop;
    end if;
    Standard_Series_Poly_Systems.Clear(srp);
  end Run_Newton;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Performs a test in standard double precision.

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
    Run_Newton(nq,idx,lp.all,sols);
  end Standard_Test;

  procedure Main is
  begin
    Standard_Test;
  end Main;

begin
  Main;
end ts_sersol;
