with Ada.Calendar;
with Timing_Package,Time_Stamps;         use Timing_Package,Time_Stamps;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Cascade_Homotopies;
with Cascade_Homotopy_Filters;
with Monodromy_Homotopies;
with Monodromy_Homotopies_io;
with Extrinsic_Diagonal_Homotopies_io;
with Path_Counts_Table;

package body Running_Cascades is

  deftol : constant double_float := 1.0E-8;     -- tolerance for zero
  defrcotol : constant double_float := 1.0E-6;  -- tolerance for regular
  defrestol : constant double_float := 1.0E-12; -- tolerance on residual
  defhomtol : constant double_float := 1.0E-6;  -- htpy membership test
  defmaxloops : constant natural32 := 20;       -- max #monodromy loops

  procedure Standard_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.Standard_Run_Cascade 1 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition
          (standard_output,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
    end if;
    Path_Counts_Table.Write_Path_Counts(standard_output,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts(standard_output,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts
          (standard_output,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(standard_output,deco);
        if Standard_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put("number of isolated solutions : ");
          put(Standard_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line;
        end if;
      end if;     
    end if;
    new_line;
    put("The CPU time for the cascade filters : ");
    if filter then
      if factor 
       then print_hms(standard_output,alltime+totfac);
       else print_hms(standard_output,alltime);
      end if;
    else
      print_hms(standard_output,totcas);
    end if;
    new_line;
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
  end Standard_Run_Cascade;

  procedure Standard_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : Standard_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.Standard_Run_Cascade 2 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition
          (standard_output,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
    end if;
    Path_Counts_Table.Write_Path_Counts(standard_output,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts
        (standard_output,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts
          (standard_output,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(standard_output,deco);
        if Standard_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put("number of isolated solutions : ");
          put(Standard_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line;
        end if;
      end if;
    end if;
    new_line;
    put("The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(standard_output,alltime+totfac);
       else print_hms(standard_output,alltime);
      end if;
    else
      print_hms(standard_output,totcas);
    end if;
    new_line;
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
  end Standard_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.DoblDobl_Run_Cascade 1 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition
          (standard_output,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
    end if;
    Path_Counts_Table.Write_Path_Counts(standard_output,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts
        (standard_output,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts
          (standard_output,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(standard_output,deco);
        if DoblDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put("number of isolated solutions : ");
          put(DoblDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line;
        end if;
      end if;
    end if;
    new_line;
    put("The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(standard_output,alltime+totfac);
       else print_hms(standard_output,alltime);
      end if;
    else
      print_hms(standard_output,totcas);
    end if;
    new_line;
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
  end DoblDobl_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.DoblDobl_Run_Cascade 2 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition
          (standard_output,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
    end if;
    Path_Counts_Table.Write_Path_Counts(standard_output,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts
        (standard_output,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts
          (standard_output,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(standard_output,deco);
        if DoblDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put("number of isolated solutions : ");
          put(DoblDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line;
        end if;
      end if;
    end if;
    new_line;
    put("The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(standard_output,alltime+totfac);
       else print_hms(standard_output,alltime);
      end if;
    else
      print_hms(standard_output,totcas);
    end if;
    new_line;
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
  end DoblDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.QuadDobl_Run_Cascade 1 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition
          (standard_output,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
    end if;
    Path_Counts_Table.Write_Path_Counts(standard_output,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts
        (standard_output,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts
          (standard_output,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(standard_output,deco);
        if QuadDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put("number of isolated solutions : ");
          put(QuadDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line;
        end if;
      end if;
    end if;
    new_line;
    put("The CPU time for the cascade filters : ");
    if filter then
      if factor 
       then print_hms(standard_output,alltime+totfac);
       else print_hms(standard_output,alltime);
      end if;
    else
      print_hms(standard_output,totcas);
    end if;
    new_line;
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
  end QuadDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.QuadDobl_Run_Cascade 2 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition
          (standard_output,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(standard_output,ep,gpts);
    end if;
    Path_Counts_Table.Write_Path_Counts(standard_output,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts
        (standard_output,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts
          (standard_output,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(standard_output,deco);
        if QuadDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put("number of isolated solutions : ");
          put(QuadDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line;
        end if;
      end if;
    end if;
    new_line;
    put("The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(standard_output,alltime+totfac);
       else print_hms(standard_output,alltime);
      end if;
    else
      print_hms(standard_output,totcas);
    end if;
    new_line;
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(standard_output,start_moment,ended_moment);
  end QuadDobl_Run_Cascade;

  procedure Standard_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.Standard_Run_Cascade 3 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (file,nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition
          (file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
    for k in gpts'range loop  -- write the witness sets to separate files
      if not Is_Null(gpts(k)) then
        Extrinsic_Diagonal_Homotopies_io.Write_Witness_Set
          (file,name,natural32(k),ep(k).all,gpts(k));
      end if;
    end loop;
    Path_Counts_Table.Write_Path_Counts(file,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts(file,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts(file,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(file,deco);
        if Standard_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put(file,"number of isolated solutions : ");
          put(file,Standard_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line(file);
        end if;
      end if;     
    end if;
    new_line(file);
    put(file,"The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(file,alltime+totfac);
       else print_hms(file,alltime);
      end if;
    else
      print_hms(file,totcas);
    end if;
    new_line(file);
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(file,start_moment,ended_moment);
  end Standard_Run_Cascade;

  procedure Standard_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : Standard_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.Standard_Run_Cascade 4 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (file,nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
    for k in gpts'range loop  -- write the witness sets to separate files
      if not Is_Null(gpts(k)) then
        Extrinsic_Diagonal_Homotopies_io.Write_Witness_Set
          (file,name,natural32(k),ep(k).all,gpts(k));
      end if;
    end loop;
    Path_Counts_Table.Write_Path_Counts(file,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts(file,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts(file,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(file,deco);
        if Standard_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put(file,"number of isolated solutions : ");
          put(file,Standard_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line(file);
        end if;
      end if;
    end if;
    new_line(file);
    put(file,"The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(file,alltime+totfac);
       else print_hms(file,alltime);
      end if;
    else
      print_hms(file,totcas);
    end if;
    new_line(file);
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(file,start_moment,ended_moment);
  end Standard_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.DoblDobl_Run_Cascade 3 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (file,nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
    for k in gpts'range loop  -- write the witness sets to separate files
      if not Is_Null(gpts(k)) then
        Extrinsic_Diagonal_Homotopies_io.Write_Witness_Set
          (file,name,natural32(k),ep(k).all,gpts(k));
      end if;
    end loop;
    Path_Counts_Table.Write_Path_Counts(file,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts(file,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts(file,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(file,deco);
        if DoblDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put(file,"number of isolated solutions : ");
          put(file,DoblDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line(file);
        end if;
      end if;
    end if;
    new_line(file);
    put(file,"The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(file,alltime+totfac);
       else print_hms(file,alltime);
      end if;
    else
      print_hms(file,totcas);
    end if;
    new_line(file);
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(file,start_moment,ended_moment);
  end DoblDobl_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.DoblDobl_Run_Cascade 4 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (file,nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
    for k in gpts'range loop  -- write the witness sets to separate files
      if not Is_Null(gpts(k)) then
        Extrinsic_Diagonal_Homotopies_io.Write_Witness_Set
          (file,name,natural32(k),ep(k).all,gpts(k));
      end if;
    end loop;
    Path_Counts_Table.Write_Path_Counts(file,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts(file,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts(file,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(file,deco);
        if DoblDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put(file,"number of isolated solutions : ");
          put(file,DoblDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line(file);
        end if;
      end if;
    end if;
    new_line(file);
    put(file,"The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(file,alltime+totfac);
       else print_hms(file,alltime);
      end if;
    else
      print_hms(file,totcas);
    end if;
    new_line(file);
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(file,start_moment,ended_moment);
  end DoblDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.QuadDobl_Run_Cascade 3 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (file,nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
    for k in gpts'range loop  -- write the witness sets to separate files
      if not Is_Null(gpts(k)) then
        Extrinsic_Diagonal_Homotopies_io.Write_Witness_Set
          (file,name,natural32(k),ep(k).all,gpts(k));
      end if;
    end loop;
    Path_Counts_Table.Write_Path_Counts(file,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts(file,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts(file,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(file,deco);
        if QuadDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put(file,"number of isolated solutions : ");
          put(file,QuadDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line(file);
        end if;
      end if;
    end if;
    new_line(file);
    put(file,"The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(file,alltime+totfac);
       else print_hms(file,alltime);
      end if;
    else
      print_hms(file,totcas);
    end if;
    new_line(file);
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(file,start_moment,ended_moment);
  end QuadDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    start_moment : constant Ada.Calendar.Time := Ada.Calendar.Clock;
    ended_moment : Ada.Calendar.Time;
    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if verbose > 0 then
      put("-> in running_cascades.QuadDobl_Run_Cascade 4 ...");
    end if;
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (file,nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
    for k in gpts'range loop  -- write the witness sets to separate files
      if not Is_Null(gpts(k)) then
        Extrinsic_Diagonal_Homotopies_io.Write_Witness_Set
          (file,name,natural32(k),ep(k).all,gpts(k));
      end if;
    end loop;
    Path_Counts_Table.Write_Path_Counts(file,pc,castm,totcas);
    if filter then
      Path_Counts_Table.Write_Filter_Counts(file,fc,filtm,totfil);
      if factor then
        Path_Counts_Table.Write_Factor_Counts(file,deco,factm,totfac);
        Path_Counts_Table.Write_Decomposition(file,deco);
        if QuadDobl_Complex_Solutions.Length_Of(gpts(0)) > 0 then
          put(file,"number of isolated solutions : ");
          put(file,QuadDobl_Complex_Solutions.Length_Of(gpts(0)),1);
          new_line(file);
        end if;
      end if;
    end if;
    new_line(file);
    put(file,"The CPU time for the cascade filters : ");
    if filter then
      if factor
       then print_hms(file,alltime+totfac);
       else print_hms(file,alltime);
      end if;
    else
      print_hms(file,totcas);
    end if;
    new_line(file);
    ended_moment := Ada.Calendar.Clock;
    Write_Elapsed_Time(file,start_moment,ended_moment);
  end QuadDobl_Run_Cascade;

  procedure Standard_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use Standard_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter_Callback
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime,Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
      filtcnt := new Standard_Natural_VecVecs.VecVec'(fc);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        idxfac := new Standard_Natural_VecVecs.Array_of_VecVecs'(deco);
      end if;
    else
      Cascade_Homotopies.Witness_Generate_Callback
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas,
         Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
    end if;
  end Standard_Cascade_Callback;

  procedure Standard_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in Standard_Complex_Laur_Systems.Laur_Sys;
                    ws : in Standard_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use Standard_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : Standard_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter_Callback
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime,Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
      filtcnt := new Standard_Natural_VecVecs.VecVec'(fc);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        idxfac := new Standard_Natural_VecVecs.Array_of_VecVecs'(deco);
      end if;
    else
      Cascade_Homotopies.Witness_Generate_Callback
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas,
         Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
    end if;
  end Standard_Cascade_Callback;

  procedure DoblDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use DoblDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter_Callback
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime,Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
      filtcnt := new Standard_Natural_VecVecs.VecVec'(fc);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        idxfac := new Standard_Natural_VecVecs.Array_of_VecVecs'(deco);
      end if;
    else
      Cascade_Homotopies.Witness_Generate_Callback
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas,
         Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
    end if;
  end DoblDobl_Cascade_Callback;

  procedure DoblDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in DoblDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use DoblDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter_Callback
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime,Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
      filtcnt := new Standard_Natural_VecVecs.VecVec'(fc);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        idxfac := new Standard_Natural_VecVecs.Array_of_VecVecs'(deco);
      end if;
    else
      Cascade_Homotopies.Witness_Generate_Callback
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas,
         Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
    end if;
  end DoblDobl_Cascade_Callback;

  procedure QuadDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use QuadDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter_Callback
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime,Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
      filtcnt := new Standard_Natural_VecVecs.VecVec'(fc);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        idxfac := new Standard_Natural_VecVecs.Array_of_VecVecs'(deco);
      end if;
    else
      Cascade_Homotopies.Witness_Generate_Callback
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas,
         Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
    end if;
  end QuadDobl_Cascade_Callback;

  procedure QuadDobl_Cascade_Callback
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean;
                pathcnt,filtcnt : out Standard_Natural_VecVecs.Link_to_VecVec;
                idxfac : out Standard_Natural_VecVecs.Link_to_Array_of_VecVecs;
                Report_Witness_Set : access procedure
                  ( ep : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    ws : in QuadDobl_Complex_Solutions.Solution_List;
                    dim : in natural32 ) ) is

    use QuadDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := deftol;
    rcotol : constant double_float := defrcotol;
    restol : constant double_float := defrestol;
    homtol : constant double_float := defhomtol;
    ep : QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := defmaxloops;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter_Callback
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,
         ep,gpts,pc,fc,castm,filtm,totcas,totfil,alltime,Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
      filtcnt := new Standard_Natural_VecVecs.VecVec'(fc);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        idxfac := new Standard_Natural_VecVecs.Array_of_VecVecs'(deco);
      end if;
    else
      Cascade_Homotopies.Witness_Generate_Callback
        (nt,embsys,sols,topdim,lowdim,restol,rcotol,ep,gpts,pc,castm,totcas,
         Report_Witness_Set);
      pathcnt := new Standard_Natural_VecVecs.VecVec'(pc);
    end if;
  end QuadDobl_Cascade_Callback;

end Running_Cascades;
