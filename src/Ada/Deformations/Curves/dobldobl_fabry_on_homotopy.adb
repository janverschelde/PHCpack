with Ada.Calendar;
with Time_Stamps;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with DoblDobl_Random_Numbers;
with Standard_Integer_Vectors;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;        use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_System_and_Solutions_io;
with DoblDobl_Homotopy;
with Standard_Parameter_Systems;
with Solution_Drops;
with DoblDobl_Homotopy_Convolutions_io;
with DoblDobl_Newton_Convolutions;
with DoblDobl_Newton_Convolution_Steps;
with Convergence_Radius_Estimates;
with Fabry_on_Homotopy_Helpers;
with Multitasked_Power_Newton;

package body DoblDobl_Fabry_on_Homotopy is

  procedure Newton_Fabry
              ( nbt : in natural32;
                cfs : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sol : in DoblDobl_Complex_Vectors.Vector ) is

    dim : constant integer32 := sol'last;
    deg : constant integer32 := cfs.deg;
    scf : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
        := DoblDobl_Newton_Convolutions.Series_Coefficients(sol,deg);
    maxit : integer32 := deg/2;
    natnbt : natural32;
    nbtasks : integer32 := integer32(nbt);
    nbrit,info : integer32 := 0;
    tol : double_float := 1.0E-20;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..dim); -- dim = #equations
    fail,verbose : boolean;
    absdx,rcond,rad,err : double_double;
    scale : constant boolean := false;
    zpt : DoblDobl_Complex_Numbers.Complex_Number;

  begin
    Fabry_on_Homotopy_Helpers.Prompt_for_Parameters(maxit,tol,verbose);
    if nbtasks = 0 then
      new_line;
      put("Give the number of tasks (0 for no multitasking) : ");
      Numbers_io.Read_Natural(natnbt); nbtasks := integer32(natnbt);
    end if;
    if verbose then
      if nbtasks = 0 then
        DoblDobl_Newton_Convolution_Steps.LU_Newton_Steps
          (standard_output,cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,
           ipvt,wrk,scale);
      else
        Multitasked_Power_Newton.DoblDobl_Run
          (standard_output,nbtasks,dim,maxit,cfs,scf,tol,true,fail,
           info,nbrit,rcond,absdx);
      end if;
    else
      if nbtasks = 0 then
        DoblDobl_Newton_Convolution_Steps.LU_Newton_Steps
          (cfs,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale,false);
      else
        Multitasked_Power_Newton.DoblDobl_Run
          (nbtasks,dim,maxit,cfs,scf,tol,true,fail,
           info,nbrit,rcond,absdx,false,false);
      end if;
    end if;
    put_line("The coefficients of the series : "); put_line(scf);
    Convergence_Radius_Estimates.Fabry
      (standard_output,scf,zpt,rad,err,fail,1,true);
    Fabry_on_Homotopy_Helpers.Write_Report(standard_output,rad,err,zpt,fail);
    DoblDobl_Complex_Vectors.Clear(wrk);
  end Newton_Fabry;

  procedure Run ( nbt : in natural32; nbequ,idxpar,deg : in integer32;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    cvh : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    ans : character;

  begin
    cvh := DoblDobl_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    while not DoblDobl_Complex_Solutions.Is_Null(tmp) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      Newton_Fabry(nbt,cvh,ls.v);
      put("Continue with the next solution ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
  end Run;

  procedure Run ( file : in file_type;
                  nbt : in natural32; nbequ,idxpar,deg : in integer32;
                  maxit : in integer32; tol : in double_float;
                  sols : in out DoblDobl_Complex_Solutions.Solution_List;
                  verbose : in boolean ) is

    cvh : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    tmp : DoblDobl_Complex_Solutions.Solution_List := sols;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution
       := DoblDobl_Complex_Solutions.Head_Of(sols);
    dim : constant integer32 := ls.v'last;
    scf : DoblDobl_Complex_VecVecs.VecVec(1..dim)
        := DoblDobl_Newton_Convolutions.Series_Coefficients(ls.v,deg);
    nbtasks : constant integer32 := integer32(nbt);
    nbrit,info : integer32 := 0;
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    wrk : DoblDobl_Complex_Vectors.Link_to_Vector
        := new DoblDobl_Complex_Vectors.Vector(1..dim); -- dim = #equations
    fail : boolean;
    absdx,rcond,rad,err : double_double;
    scale : constant boolean := false;
    zpt : DoblDobl_Complex_Numbers.Complex_Number;
    tstart,tstop : Ada.Calendar.Time;
    cnt : integer32 := 0;

  begin
    cvh := DoblDobl_Homotopy_Convolutions_io.Make_Homotopy(nbequ,idxpar,deg);
    tstart := Ada.Calendar.Clock;
    loop
      if verbose then
        if nbtasks = 0 then
          DoblDobl_Newton_Convolution_Steps.LU_Newton_Steps
            (file,cvh,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale);
        else
          Multitasked_Power_Newton.DoblDobl_Run
            (file,nbtasks,dim,maxit,cvh,scf,tol,true,fail,
             info,nbrit,rcond,absdx);
        end if;
      else
        if nbtasks = 0 then
          DoblDobl_Newton_Convolution_Steps.LU_Newton_Steps
            (cvh,scf,maxit,nbrit,tol,absdx,fail,rcond,ipvt,wrk,scale,false);
        else
          Multitasked_Power_Newton.DoblDobl_Run
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
      Convergence_Radius_Estimates.Fabry(file,scf,zpt,rad,err,fail,1,true);
      Fabry_on_Homotopy_Helpers.Write_Report(file,rad,err,zpt,fail);
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
      exit when DoblDobl_Complex_Solutions.Is_Null(tmp);
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      scf := DoblDobl_Newton_Convolutions.Series_Coefficients(ls.v,deg);
    end loop;
    tstop := Ada.Calendar.Clock;
    new_line(file);
    Time_Stamps.Write_Elapsed_Time(file,tstart,tstop);
    DoblDobl_Complex_Vectors.Clear(wrk);
    DoblDobl_Speelpenning_Convolutions.Clear(cvh);
  end Run;

  procedure Artificial_Setup
              ( nbtasks : in natural32; vrblvl : in integer32 := 0 ) is

    target,start : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;
    posdeg : positive;
    nbequ,nbvar,nbsols,deg,mxt : integer32 := 0;
    ans : character;
    tofile : boolean;
    outfile : file_type;
    nbt : natural32 := nbtasks;
    tol : double_float := 1.0E-12;
    vrb : boolean;

    use DoblDobl_Complex_Polynomials;

  begin
    if vrblvl > 0 then
      put_line("-> in dobldobl_fabry_on_homotopy.Artificial_Setup ...");
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
       then gamma := DoblDobl_Random_Numbers.Random1;
       else gamma := DoblDobl_Complex_Numbers.Create(integer(1));
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
        DoblDobl_System_and_Solutions_io.put
          (outfile,start.all,sols,"THE START SOLUTIONS :");
      end if;
      DoblDobl_Homotopy.Create(target.all,start.all,1,gamma);
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

    hom : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,dropsols : DoblDobl_Complex_Solutions.Solution_List;
    posdeg : positive;
    nbequ,sysnbvar,solnbvar,nbsols,idxpar,deg,mxt : integer32 := 0;
    par : Standard_Integer_Vectors.Vector(1..1);
    ans : character;
    tofile : boolean;
    outfile : file_type;
    nbt : natural32 := nbtasks;
    tol : double_float := 1.0E-12;
    vrb : boolean;

    use DoblDobl_Complex_Polynomials;

  begin
    if vrblvl > 0 then
      put_line("-> in dobldobl_fabry_on_homotopy.Natural_Setup ...");
    end if;
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
          DoblDobl_System_and_Solutions_io.put
            (outfile,hom.all,sols,"THE START SOLUTIONS :");
        else
          DoblDobl_System_and_Solutions_io.put
            (outfile,hom.all,dropsols,"THE START SOLUTIONS :");
        end if;
      end if;
      DoblDobl_Homotopy.Create(hom.all,idxpar);
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
    if vrblvl > 0
     then put_line("-> in dobldobl_fabry_on_homotopy.Main ...");
    end if;
    new_line;
    put("Artificial-parameter homotopy ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Artificial_Setup(nbtasks,vrblvl-1);
     else Natural_Setup(nbtasks,vrblvl-1);
    end if;
  end Main;

end DoblDobl_Fabry_on_Homotopy;
