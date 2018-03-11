with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_VecVecs;
with Cascade_Homotopies;
with Cascade_Homotopy_Filters;
with Monodromy_Homotopies;
with Monodromy_Homotopies_io;
with Path_Counts_Table;

package body Running_Cascades is

  procedure Standard_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use Standard_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
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
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
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
  end Standard_Run_Cascade;

  procedure Standard_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use Standard_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : Standard_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
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
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
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
  end Standard_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
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
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
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
  end DoblDobl_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
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
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
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
  end DoblDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
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
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
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
  end QuadDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
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
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
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
  end QuadDobl_Run_Cascade;

  procedure Standard_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use Standard_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
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
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
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
  end Standard_Run_Cascade;

  procedure Standard_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use Standard_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : Standard_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
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
  end Standard_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
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
  end DoblDobl_Run_Cascade;

  procedure DoblDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use DoblDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
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
  end DoblDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
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
  end QuadDobl_Run_Cascade;

  procedure QuadDobl_Run_Cascade
              ( file : in file_type; name : in string;
                nt,topdim,lowdim : in natural32;
                embsys : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                filter,factor : in boolean ) is

    use QuadDobl_Complex_Solutions;

    ns : constant integer32 := integer32(topdim);
    tol : constant double_float := 1.0E-8;
    rcotol : constant double_float := 1.0E-6;
    restol : constant double_float := 1.0E-6;
    homtol : constant double_float := 1.0E-6;
    ep : QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys(0..ns);
    gpts : Array_of_Solution_Lists(0..ns);
    pc,fc : Standard_Natural_VecVecs.VecVec(0..ns);
    nbl : constant natural32 := 20;
    deco : Standard_Natural_VecVecs.Array_of_VecVecs(1..ns);
    castm : Array_of_Duration(0..integer(ns));
    filtm : Array_of_Duration(0..integer(ns));
    factm : Array_of_Duration(0..integer(ns));
    totcas,totfil,totfac,alltime : duration;

  begin
    if filter then
      Cascade_Homotopy_Filters.Witness_Filter
        (nt,embsys,sols,topdim,lowdim,tol,rcotol,restol,homtol,ep,gpts,pc,fc,
         castm,filtm,totcas,totfil,alltime);
      if factor then
        Monodromy_Homotopies.Witness_Factor
          (false,ep,gpts,topdim,nbl,tol,deco,factm,totfac);
        Monodromy_Homotopies_io.Write_Decomposition(file,ep,gpts,deco);
      else
        Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
      end if;
    else
      Cascade_Homotopies.Witness_Generate
        (nt,embsys,sols,topdim,lowdim,tol,ep,gpts,pc,castm,totcas);
      Monodromy_Homotopies_io.Write_Components(file,ep,gpts);
    end if;
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
  end QuadDobl_Run_Cascade;

end Running_Cascades;
