with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Random_Numbers;
with DoblDobl_Random_Numbers;
with QuadDobl_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;        use QuadDobl_Complex_VecVecs_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Vector_Splitters;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_Parameter_Systems;
with Solution_Drops;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Coefficient_Convolutions;
with Standard_Convolution_Splitters;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy_Convolutions_io;
with Standard_Newton_Convolutions;
with DoblDobl_Newton_Convolutions;
with QuadDobl_Newton_Convolutions;
with Newton_Power_Convolutions;
with Staggered_Newton_Convolutions;
with Convergence_Radius_Estimates;

procedure ts_fabryhom is

-- DESCRIPTION :
--   Development of a test on the Newton-Fabry convergence radius
--   for artificial or natural-parameter homotopies.

  procedure Standard_Newton_Fabry
              ( cfs : in Standard_Coefficient_Convolutions.Link_to_System;
                sol : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs,
  --   in double precision, on coefficient convolution systems,
  --   with a staggered, inlined, and indexed Newton's method.

    dim : constant integer32 := sol'last;
    deg : constant integer32 := cfs.deg;
    scf : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Standard_Newton_Convolutions.Series_Coefficients(sol,deg);
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rc,ic,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    tol : double_float := 1.0E-12;
    ans : character;
    rcond,absdx,rad,err : double_float;
    maxit,nbrit,idxtoldx : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    fail : boolean;
    scale : constant boolean := false;
    zpt : Standard_Complex_Numbers.Complex_Number;

  begin
    Standard_Floating_VecVecVecs.Allocate(rv,1,deg,1,dim,1,dim);
    Standard_Floating_VecVecVecs.Allocate(iv,1,deg,1,dim,1,dim);
    rc := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ic := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rb := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ib := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Tolerance for the accuracy : "); put(tol,3); new_line;
      put("Change the tolerance ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      if ans = 'y' then
        put("Give the new tolerance for the accuracy : "); get(tol);
      end if;
    end loop;
    Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
      (standard_output,cfs,scf,rx,ix,maxit,nbrit,tol,idxtoldx,absdx,
       fail,rcond,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
    put_line("The coefficients of the series : "); put_line(scf);
    Convergence_Radius_Estimates.Fabry
      (standard_output,scf,zpt,rad,err,fail,1,true);
    put("the convergence radius :"); put(rad,3);
    put("   error estimate :"); put(err,3); new_line;
    put(zpt); put_line("  estimates nearest singularity");
    if fail
     then put_line("Reported failure.");
     else put_line("Reported success.");
    end if;
    Standard_Floating_Vectors.Clear(ry);
    Standard_Floating_Vectors.Clear(iy);
    Standard_Floating_VecVecs.Deep_Clear(rc);
    Standard_Floating_VecVecs.Deep_Clear(ic);
    Standard_Floating_VecVecs.Deep_Clear(rb);
    Standard_Floating_VecVecs.Deep_Clear(ib);
  end Standard_Newton_Fabry;

  procedure DoblDobl_Newton_Fabry
              ( cfs : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs,
  --   in double double precision.

    dim : constant integer32 := sol'last;
    deg : constant integer32 := cfs.deg;
    scf : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
        := DoblDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
    maxit,nbrit : integer32 := 0;
    ans : character;
    tol : double_float := 1.0E-20;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..dim); -- dim = #equations
    fail : boolean;
    absdx,rcond,rad,err : double_double;
    scale : constant boolean := false;
    zpt : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Tolerance for the accuracy : "); put(tol,3); new_line;
      put("Change the tolerance ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      if ans = 'y' then
        put("Give the new tolerance for the accuracy : "); get(tol);
      end if;
    end loop;
    Newton_Power_Convolutions.LU_Newton_Steps
      (standard_output,cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,
       ipvt,wrk,scale);
    put_line("The coefficients of the series : "); put_line(scf);
    Convergence_Radius_Estimates.Fabry
      (standard_output,scf,zpt,rad,err,fail,1,true);
    put("the convergence radius : "); put(rad,3);
    put("   error estimate : "); put(err,3); new_line;
    put(zpt); put_line("  estimates nearest singularity");
    if fail
     then put_line("Reported failure.");
     else put_line("Reported success.");
    end if;
    DoblDobl_Complex_Vectors.Clear(wrk);
  end DoblDobl_Newton_Fabry;

  procedure QuadDobl_Newton_Fabry
              ( cfs : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Runs Newton's method and applies Fabry's theorem
  --   starting at the solution for the homotopy in cfs,
  --   in quad double precision.

    dim : constant integer32 := sol'last;
    deg : constant integer32 := cfs.deg;
    scf : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
        := QuadDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
    maxit,nbrit : integer32 := 0;
    ans : character;
    tol : double_float := 1.0E-32;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : QuadDobl_Complex_Vectors.Link_to_Vector
        := new QuadDobl_Complex_Vectors.Vector(1..dim); -- dim = #equations
    fail : boolean;
    absdx,rcond,rad,err : quad_double;
    scale : constant boolean := false;
    zpt : QuadDobl_Complex_Numbers.Complex_Number;

  begin
    new_line;
    put("Give the maximum number of iterations : "); get(maxit);
    loop
      put("Tolerance for the accuracy : "); put(tol,3); new_line;
      put("Change the tolerance ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      if ans = 'y' then
        put("Give the new tolerance for the accuracy : "); get(tol);
      end if;
    end loop;
    Newton_Power_Convolutions.LU_Newton_Steps
      (standard_output,cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,
       ipvt,wrk,scale);
    put_line("The coefficients of the series : "); put_line(scf);
    Convergence_Radius_Estimates.Fabry
      (standard_output,scf,zpt,rad,err,fail,1,true);
    put("the convergence radius : "); put(rad,3);
    put("   error estimate : "); put(err,3); new_line;
    put(zpt); put_line("  estimates nearest singularity");
    if fail
     then put_line("Reported failure.");
     else put_line("Reported success.");
    end if;
    QuadDobl_Complex_Vectors.Clear(wrk);
  end QuadDobl_Newton_Fabry;

  procedure Standard_Run
              ( nbequ,idxpar,deg : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With the homotopy defined in double precision,
  --   starting at the solutions in sols, runs Newton's method
  --   on power series and applies Fabry's theorem.

  -- ON ENTRY :
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

    cvh : Standard_Speelpenning_Convolutions.Link_to_System;
    cfh : Standard_Coefficient_Convolutions.Link_to_System;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    ans : character;

  begin
    cvh := Standard_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    cfh := Standard_Convolution_Splitters.Split(cvh);
    while not Standard_Complex_Solutions.Is_Null(tmp) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      Standard_Newton_Fabry(cfh,ls.v);
      put("Continue with the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Standard_Run;

  procedure DoblDobl_Run
              ( nbequ,idxpar,deg : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With the homotopy defined in double double precision,
  --   starting at the solutions in sols, runs Newton's method
  --   on power series and applies Fabry's theorem.

  -- ON ENTRY :
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

    cvh : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    ans : character;

  begin
    cvh := DoblDobl_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      DoblDobl_Newton_Fabry(cvh,ls.v);
      put("Continue with the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end DoblDobl_Run;

  procedure QuadDobl_Run
              ( nbequ,idxpar,deg : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   With the homotopy defined in quad double precision,
  --   starting at the solutions in sols, runs Newton's method
  --   on power series and applies Fabry's theorem.

  -- ON ENTRY :
  --   nbequ    number of equations in the homotopy;
  --   idxpar   index of the continuation parameter in the homotopy;
  --   deg      degree of the power series;
  --   sols     start solutions in the homotopy.

    cvh : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : QuadDobl_Complex_Solutions.Solution_List := sols;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    ans : character;

  begin
    cvh := QuadDobl_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    while not QuadDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      QuadDobl_Newton_Fabry(cvh,ls.v);
      put("Continue with the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end QuadDobl_Run;

  procedure Standard_Artificial_Test is

  -- DESCRIPTION :
  --   Promps the user for an artifical-parameter homotopy.
  --   If the number of start solutions is positive,
  --   then the homotopy is defined, in double precision.

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;
    nbequ,nbvar,nbsols,deg : integer32 := 0;
    ans : character;

    use Standard_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for the target system ...");
    get(target);
    nbequ := target'last;
    nbvar := integer32(Number_of_Unknowns(target(target'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(nbvar,1); put_line(" variables.");
    new_line;
    put_line
      ("Reading the name of a file for the start system and solutions ...");
    Standard_System_and_Solutions_io.get(start,sols);
    nbsols := integer32(Standard_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      nbvar := Standard_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(nbvar,1); put_line(".");
      new_line;
      put("Random gamma ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then gamma := Standard_Random_Numbers.Random1;
       else gamma := Standard_Complex_Numbers.Create(1.0);
      end if;
      Standard_Homotopy.Create(target.all,start.all,1,gamma);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      Standard_Run(nbequ,nbvar+1,deg,sols);
    end if;
  end Standard_Artificial_Test;

  procedure DoblDobl_Artificial_Test is

  -- DESCRIPTION :
  --   Promps the user for an artifical-parameter homotopy.
  --   If the number of start solutions is positive,
  --   then the homotopy is defined, in double double precision.

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;
    nbequ,nbvar,nbsols,deg : integer32 := 0;
    ans : character;

    use DoblDobl_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for the target system ...");
    get(target);
    nbequ := target'last;
    nbvar := integer32(Number_of_Unknowns(target(target'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(nbvar,1); put_line(" variables.");
    new_line;
    put_line
      ("Reading the name of a file for the start system and solutions ...");
    DoblDobl_System_and_Solutions_io.get(start,sols);
    nbsols := integer32(DoblDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      nbvar := DoblDobl_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(nbvar,1); put_line(".");
      new_line;
      put("Random gamma ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then gamma := DoblDobl_Random_Numbers.Random1;
       else gamma := DoblDobl_Complex_Numbers.Create(integer(1));
      end if;
      DoblDobl_Homotopy.Create(target.all,start.all,1,gamma);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      DoblDobl_Run(nbequ,nbvar+1,deg,sols);
    end if;
  end DoblDobl_Artificial_Test;

  procedure QuadDobl_Artificial_Test is

  -- DESCRIPTION :
  --   Promps the user for an artifical-parameter homotopy.
  --   If the number of start solutions is positive,
  --   then the homotopy is defined, in quad double precision.

    target,start : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;
    nbequ,nbvar,nbsols,deg : integer32 := 0;
    ans : character;

    use QuadDobl_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for the target system ...");
    get(target);
    nbequ := target'last;
    nbvar := integer32(Number_of_Unknowns(target(target'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(nbvar,1); put_line(" variables.");
    new_line;
    put_line
      ("Reading the name of a file for the start system and solutions ...");
    QuadDobl_System_and_Solutions_io.get(start,sols);
    nbsols := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      nbvar := QuadDobl_Complex_Solutions.Head_Of(sols).n;
      new_line;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(nbvar,1); put_line(".");
      new_line;
      put("Random gamma ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then gamma := QuadDobl_Random_Numbers.Random1;
       else gamma := QuadDobl_Complex_Numbers.Create(integer(1));
      end if;
      QuadDobl_Homotopy.Create(target.all,start.all,1,gamma);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      QuadDobl_Run(nbequ,nbvar+1,deg,sols);
    end if;
  end QuadDobl_Artificial_Test;

  procedure Standard_Natural_Test is

  -- DESCRIPTION :
  --   Promps the user for a natural-parameter homotopy,
  --   with start solutions in double precision.

    hom : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : Standard_Complex_Solutions.Solution_List;
    nbequ,sysnbvar,solnbvar,nbsols,idxpar,deg : integer32 := 0;
    par : Standard_Integer_Vectors.Vector(1..1);

    use Standard_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for a homotopy ...");
    Standard_System_and_Solutions_io.get(hom,sols);
    nbequ := hom'last;
    sysnbvar := integer32(Number_of_Unknowns(hom(hom'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(sysnbvar,1); put_line(" variables.");
    nbsols := integer32(Standard_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      solnbvar := Standard_Complex_Solutions.Head_Of(sols).n;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(solnbvar,1); put_line(".");
      new_line;
      par := Standard_Parameter_Systems.Define_Parameters(nbequ,sysnbvar,1);
      idxpar := par(1);
      put("The index to the continuation parameter : ");
      put(idxpar,1); new_line;
      if solnbvar = nbequ then
        put_line("Solution dimension is okay.");
      else
        put_line("Need to drop one coordinate of each solution.");
        dropsols := Solution_Drops.Drop(sols,natural32(idxpar));
      end if;
      Standard_Homotopy.Create(hom.all,idxpar);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      if solnbvar = nbequ
       then Standard_Run(nbequ,idxpar,deg,sols);
       else Standard_Run(nbequ,idxpar,deg,dropsols);
      end if;
    end if;
  end Standard_Natural_Test;

  procedure DoblDobl_Natural_Test is

  -- DESCRIPTION :
  --   Promps the user for a natural-parameter homotopy,
  --   with start solutions in double double precision.

    hom : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : DoblDobl_Complex_Solutions.Solution_List;
    nbequ,sysnbvar,solnbvar,nbsols,idxpar,deg : integer32 := 0;
    par : Standard_Integer_Vectors.Vector(1..1);

    use DoblDobl_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for a homotopy ...");
    DoblDobl_System_and_Solutions_io.get(hom,sols);
    nbequ := hom'last;
    sysnbvar := integer32(Number_of_Unknowns(hom(hom'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(sysnbvar,1); put_line(" variables.");
    nbsols := integer32(DoblDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      solnbvar := DoblDobl_Complex_Solutions.Head_Of(sols).n;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(solnbvar,1); put_line(".");
      new_line;
      par := Standard_Parameter_Systems.Define_Parameters(nbequ,sysnbvar,1);
      idxpar := par(1);
      put("The index to the continuation parameter : ");
      put(idxpar,1); new_line;
      if solnbvar = nbequ then
        put_line("Solution dimension is okay.");
      else
        put_line("Need to drop one coordinate of each solution.");
        dropsols := Solution_Drops.Drop(sols,natural32(idxpar));
      end if;
      DoblDobl_Homotopy.Create(hom.all,idxpar);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      if solnbvar = nbequ
       then DoblDobl_Run(nbequ,idxpar,deg,sols);
       else DoblDobl_Run(nbequ,idxpar,deg,dropsols);
      end if;
    end if;
  end DoblDobl_Natural_Test;

  procedure QuadDobl_Natural_Test is

  -- DESCRIPTION :
  --   Promps the user for a natural-parameter homotopy,
  --   with start solutions in quad double precision.

    hom : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : QuadDobl_Complex_Solutions.Solution_List;
    nbequ,sysnbvar,solnbvar,nbsols,idxpar,deg : integer32 := 0;
    par : Standard_Integer_Vectors.Vector(1..1);

    use QuadDobl_Complex_Polynomials;

  begin
    new_line;
    put_line("Reading the name of a file for a homotopy ...");
    QuadDobl_System_and_Solutions_io.get(hom,sols);
    nbequ := hom'last;
    sysnbvar := integer32(Number_of_Unknowns(hom(hom'first)));
    new_line;
    put("Read "); put(nbequ,1); put(" polynomials in ");
    put(sysnbvar,1); put_line(" variables.");
    nbsols := integer32(QuadDobl_Complex_Solutions.Length_Of(sols));
    if nbsols = 0 then
      put_line("No solutions read.");
    else
      solnbvar := QuadDobl_Complex_Solutions.Head_Of(sols).n;
      put("Read "); put(nbsols,1); put(" solutions in dimension ");
      put(solnbvar,1); put_line(".");
      new_line;
      par := Standard_Parameter_Systems.Define_Parameters(nbequ,sysnbvar,1);
      idxpar := par(1);
      put("The index to the continuation parameter : ");
      put(idxpar,1); new_line;
      if solnbvar = nbequ then
        put_line("Solution dimension is okay.");
      else
        put_line("Need to drop one coordinate of each solution.");
        dropsols := Solution_Drops.Drop(sols,natural32(idxpar));
      end if;
      QuadDobl_Homotopy.Create(hom.all,idxpar);
      new_line;
      put("Give the degree of the power series : "); get(deg);
      if solnbvar = nbequ
       then QuadDobl_Run(nbequ,idxpar,deg,sols);
       else QuadDobl_Run(nbequ,idxpar,deg,dropsols);
      end if;
    end if;
  end QuadDobl_Natural_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision
  --   and the type of homotopy.

    prc,ans : character;

  begin
    prc := Communications_with_User.Prompt_for_Precision;
    new_line;
    put("Artificial-parameter homotopy ? (y/n) "); Ask_Yes_or_No(ans);
    case prc is
      when '0' =>
        if ans = 'y'
         then Standard_Artificial_Test;
         else Standard_Natural_Test;
        end if;
      when '1' =>
        if ans = 'y'
         then DoblDobl_Artificial_Test;
         else DoblDobl_Natural_Test;
        end if;
      when '2' =>
        if ans = 'y'
         then QuadDobl_Artificial_Test;
         else QuadDobl_Natural_Test;
        end if;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_fabryhom;
