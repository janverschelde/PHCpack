with Ada.Calendar;
with Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecVecs;
with Standard_Vector_Splitters;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_System_and_Solutions_io;
with Standard_Homotopy;
with Standard_Parameter_Systems;
with Solution_Drops;
with Standard_Convolution_Splitters;
with Standard_Homotopy_Convolutions_io;
with Standard_Newton_Convolutions;
with Staggered_Newton_Convolutions;
with Convergence_Radius_Estimates;
with Fabry_on_Homotopy_Helpers;
with Multitasked_Power_Newton;

package body Standard_Fabry_on_Homotopy is

  procedure Newton_Fabry
              ( nbt : in natural32;
                cvh : in Standard_Speelpenning_Convolutions.Link_to_System;
                cfs : in Standard_Coefficient_Convolutions.Link_to_System;
                sol : in Standard_Complex_Vectors.Vector ) is

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
    rcond,absdx,rad,err : double_float;
    maxit : integer32 := deg/2;
    natnbt : natural32;
    nbtasks : integer32 := integer32(nbt);
    nbrit,idxtoldx,info : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    fail,verbose : boolean;
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
    Fabry_on_Homotopy_Helpers.Prompt_for_Parameters(maxit,tol,verbose);
    if nbtasks = 0 then
      new_line;
      put("Give the number of tasks (0 for no multitasking) : ");
      Numbers_io.Read_Natural(natnbt); nbtasks := integer32(natnbt);
    end if;
    if verbose then
      if nbtasks = 0 then
        Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
          (standard_output,cfs,scf,rx,ix,maxit,nbrit,tol,idxtoldx,absdx,
           fail,rcond,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
      else
        Multitasked_Power_Newton.Standard_Run
          (standard_output,nbtasks,dim,maxit,cvh,scf,tol,true,fail,
           info,nbrit,rcond,absdx);
      end if;
    else
      if nbtasks = 0 then
        Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
          (cfs,scf,rx,ix,maxit,nbrit,tol,idxtoldx,absdx,
           fail,rcond,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale,false);
      else
        Multitasked_Power_Newton.Standard_Run
          (nbtasks,dim,maxit,cvh,scf,tol,true,fail,
           info,nbrit,rcond,absdx,false,false);
      end if;
    end if;
    put_line("The coefficients of the series : "); put_line(scf);
    Convergence_Radius_Estimates.Fabry
      (standard_output,scf,zpt,rad,err,fail,1,true);
    Fabry_on_Homotopy_Helpers.Write_Report(standard_output,rad,err,zpt,fail);
    Standard_Floating_Vectors.Clear(ry);
    Standard_Floating_Vectors.Clear(iy);
    Standard_Floating_VecVecs.Deep_Clear(rc);
    Standard_Floating_VecVecs.Deep_Clear(ic);
    Standard_Floating_VecVecs.Deep_Clear(rb);
    Standard_Floating_VecVecs.Deep_Clear(ib);
  end Newton_Fabry;

  procedure Run ( nbt : in natural32; nbequ,idxpar,deg : in integer32;
                  sols : in out Standard_Complex_Solutions.Solution_List ) is

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
      Newton_Fabry(nbt,cvh,cfh,ls.v);
      put("Continue with the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Run;

  procedure Run ( file : in file_type;
                  nbt : in natural32; nbequ,idxpar,deg : in integer32;
                  maxit : in integer32; tol : in double_float;
                  sols : in out Standard_Complex_Solutions.Solution_List;
                  verbose : in boolean ) is

    cvh : Standard_Speelpenning_Convolutions.Link_to_System;
    cfh : Standard_Coefficient_Convolutions.Link_to_System;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution
       := Standard_Complex_Solutions.Head_Of(sols);
    dim : constant integer32 := ls.v'last;
    scf : Standard_Complex_VecVecs.VecVec(1..dim)
        := Standard_Newton_Convolutions.Series_Coefficients(ls.v,deg);
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    rc,ic,rb,ib : Standard_Floating_VecVecs.Link_to_VecVec;
    ry,iy : Standard_Floating_Vectors.Link_to_Vector;
    rv,iv : Standard_Floating_VecVecVecs.Link_to_VecVecVec;
    nbtasks : constant integer32 := integer32(nbt);
    nbrit,idxtoldx,info : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : Standard_Complex_Vectors.Link_to_Vector
        := new Standard_Complex_Vectors.Vector(1..dim); -- dim = #equations
    fail : boolean;
    absdx,rcond,rad,err : double_float;
    scale : constant boolean := false;
    zpt : Standard_Complex_Numbers.Complex_Number;
    tstart,tstop : Ada.Calendar.Time;
    cnt : integer32 := 0;

  begin
    Standard_Floating_VecVecVecs.Allocate(rv,1,deg,1,dim,1,dim);
    Standard_Floating_VecVecVecs.Allocate(iv,1,deg,1,dim,1,dim);
    rc := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    ic := Standard_Vector_Splitters.Allocate(dim,dim,1,1);
    rb := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ib := Standard_Vector_Splitters.Allocate(deg,dim,0,1);
    ry := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    iy := new Standard_Floating_Vectors.Vector'(1..dim => 0.0);
    cvh := Standard_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    cfh := Standard_Convolution_Splitters.Split(cvh);
    tstart := Ada.Calendar.Clock;
    loop
      if verbose then
        if nbtasks = 0 then
          Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
            (file,cfh,scf,rx,ix,maxit,nbrit,tol,idxtoldx,absdx,
             fail,rcond,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale);
        else
          Multitasked_Power_Newton.Standard_Run
            (file,nbtasks,dim,maxit,cvh,scf,tol,true,fail,
             info,nbrit,rcond,absdx);
        end if;
      else
        if nbtasks = 0 then
          Staggered_Newton_Convolutions.Indexed_LU_Newton_Steps
            (cfh,scf,rx,ix,maxit,nbrit,tol,idxtoldx,absdx,
             fail,rcond,ipvt,rc,ic,rv,iv,rb,ib,ry,iy,scale,false);
        else
          Multitasked_Power_Newton.Standard_Run
            (nbtasks,dim,maxit,cvh,scf,tol,true,fail,
             info,nbrit,rcond,absdx,false,false);
        end if;
      end if;
      put(file,"The coefficients of the series for solution ");
      cnt := cnt + 1; put(file,cnt,1); put_line(file," :");
      put_line(file,scf);
      new_line(file);
      put(file,"rcond : "); put(file,rcond,3); new_line(file);
      put(file,"absdx : "); put(file,absdx,3); new_line(file);
      put(file,"nbrit : "); put(file,nbrit,1);
      if fail
       then put(file,"  failed to reach ");
       else put(file,"  reached tolerance ");
      end if;
      put(file,tol,3); new_line(file); flush(file);
      new_line(file);
      Convergence_Radius_Estimates.Fabry
        (file,scf,zpt,rad,err,fail,1,true);
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
      exit when Standard_Complex_Solutions.Is_Null(tmp);
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      scf := Standard_Newton_Convolutions.Series_Coefficients(ls.v,deg);
    end loop;
    tstop := Ada.Calendar.Clock;
    new_line(file);
    Time_Stamps.Write_Elapsed_Time(file,tstart,tstop);
    Standard_Floating_Vectors.Clear(ry);
    Standard_Floating_Vectors.Clear(iy);
    Standard_Floating_VecVecs.Deep_Clear(rc);
    Standard_Floating_VecVecs.Deep_Clear(ic);
    Standard_Floating_VecVecs.Deep_Clear(rb);
    Standard_Floating_VecVecs.Deep_Clear(ib);
    Standard_Complex_Vectors.Clear(wrk);
    Standard_Speelpenning_Convolutions.Clear(cvh);
  end Run;

  procedure Artificial_Setup
              ( nbtasks : in natural32; vrblvl : in integer32 := 0 ) is

    target,start : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    gamma : Standard_Complex_Numbers.Complex_Number;
    posdeg : positive;
    nbequ,nbvar,nbsols,deg,mxt : integer32 := 0;
    ans : character;
    tofile : boolean;
    outfile : file_type;
    nbt : natural32 := nbtasks;
    tol : double_float := 1.0E-12;
    vrb : boolean;

    use Standard_Complex_Polynomials;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_fabry_on_homotopy.Artificial_Setup ...");
    end if;
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
      put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
      tofile := (ans = 'y');
      if tofile then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(outfile);
      end if;
      new_line;
      put("Random gamma ? (y/n) "); Ask_Yes_or_No(ans);
      if ans = 'y'
       then gamma := Standard_Random_Numbers.Random1;
       else gamma := Standard_Complex_Numbers.Create(1.0);
      end if;
      new_line;
      put("Give the degree of the power series : ");
      Numbers_io.Read_Positive(posdeg); deg := integer32(posdeg);
      if tofile then
        mxt := deg/2;
        Fabry_on_Homotopy_Helpers.Prompt_and_Write(outfile,nbt,mxt,tol,vrb);
        put(outfile,"gamma : "); put(outfile,gamma); new_line(outfile);
        put(outfile,"degree : "); put(outfile,deg,1); new_line(outfile);
        new_line(outfile);
        put(outfile,target'last); new_line(outfile);
        put(outfile,target.all);
        new_line(outfile);
        put_line(outfile,"THE START SYSTEM :");
        Standard_System_and_Solutions_io.put
          (outfile,start.all,sols,"THE START SOLUTIONS :");
      end if;
      Standard_Homotopy.Create(target.all,start.all,1,gamma);
      if not tofile then
        Run(nbtasks,nbequ,nbvar+1,deg,sols);
      else
        new_line(outfile);
        Run(outfile,nbt,nbequ,nbvar+1,deg,mxt,tol,sols,vrb);
      end if;
    end if;
  end Artificial_Setup;

  procedure Natural_Setup
              ( nbtasks : in natural32; vrblvl : in integer32 := 0 ) is

    hom : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : Standard_Complex_Solutions.Solution_List;
    posdeg : positive;
    nbequ,sysnbvar,solnbvar,nbsols,idxpar,deg,mxt : integer32 := 0;
    par : Standard_Integer_Vectors.Vector(1..1);
    ans : character;
    tofile : boolean;
    outfile : file_type;
    nbt : natural32 := nbtasks;
    tol : double_float := 1.0E-12;
    vrb : boolean;

    use Standard_Complex_Polynomials;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_fabry_on_homotopy.Natural_Setup ...");
    end if;
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
        put_line("Dropping one coordinate of each solution ...");
        dropsols := Solution_Drops.Drop(sols,natural32(idxpar));
      end if;
      new_line;
      put("Output to file ? (y/n) "); Ask_Yes_or_No(ans);
      tofile := (ans = 'y');
      if tofile then
        new_line;
        put_line("Reading the name of the output file ...");
        Read_Name_and_Create_File(outfile);
      end if;
      new_line;
      put("Give the degree of the power series : ");
      Numbers_io.Read_Positive(posdeg); deg := integer32(posdeg);
      if tofile then
        mxt := deg/2;
        Fabry_on_Homotopy_Helpers.Prompt_and_Write(outfile,nbt,mxt,tol,vrb);
        put(outfile,"degree : "); put(outfile,deg,1); new_line(outfile);
        new_line(outfile);
        if solnbvar = nbequ then
          Standard_System_and_Solutions_io.put
            (outfile,hom.all,sols,"THE START SOLUTIONS :");
        else
          Standard_System_and_Solutions_io.put
            (outfile,hom.all,dropsols,"THE START SOLUTIONS :");
        end if;
      end if;
      Standard_Homotopy.Create(hom.all,idxpar);
      if not tofile then
        if solnbvar = nbequ
         then Run(nbtasks,nbequ,idxpar,deg,sols);
         else Run(nbtasks,nbequ,idxpar,deg,dropsols);
        end if;
      else
        new_line(outfile);
        if solnbvar = nbequ
         then Run(outfile,nbt,nbequ,idxpar,deg,mxt,tol,sols,vrb);
         else Run(outfile,nbt,nbequ,idxpar,deg,mxt,tol,dropsols,vrb);
        end if;
      end if;
    end if;
  end Natural_Setup;

  procedure Main ( nbtasks : in natural32; vrblvl : in integer32 := 0 ) is

    ans : character;

  begin
    if vrblvl > 0 then
      put_line("-> in standard_fabry_on_homotopy.Main ...");
    end if;
    new_line;
    put("Artificial-parameter homotopy ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Artificial_Setup(nbtasks,vrblvl-1);
     else Natural_Setup(nbtasks,vrblvl-1);
    end if;
  end Main;

end Standard_Fabry_on_Homotopy;
