with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Standard_Complex_VecVecs;
with Multprec_Complex_Poly_Systems;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Sampling_Machine;
with Sampling_Laurent_Machine;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with DoblDobl_Sample_Lists;              use DoblDobl_Sample_Lists;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with QuadDobl_Sample_Lists;              use QuadDobl_Sample_Lists;
with Irreducible_Components;             use Irreducible_Components;
with Irreducible_Component_Lists;        use Irreducible_Component_Lists;
with Irreducible_Component_Lists_io;     use Irreducible_Component_Lists_io;
with Drivers_to_Component_Creators;      use Drivers_to_Component_Creators;
with Monodromy_Component_Breakup;        use Monodromy_Component_Breakup;

package body Drivers_to_Factor_Components is

  procedure Call_Standard_Interpolate
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim : in natural32; itp : in natural32 ) is

    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    deco,deco_last : Standard_Irreducible_Component_List;
    hyp : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    full_output : boolean := false;
    stoptol,membtol : double_float := 1.0E-8;

  begin
    Standard_Tuner(file,full_output,stoptol,membtol);
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    new_line(file);
    put_line(file,"Diagnostics of filtering process :");
    Sampling_Machine.Initialize(p);
    case itp is
      when 0 => Standard_Massive_Interpolate
                  (file,full_output,sps,stoptol,membtol,deco,deco_last);
      when 1 => Standard_Incremental_Interpolate
                  (file,full_output,sps,stoptol,membtol,deco,deco_last);
      when 2 => Standard_Incremental_Interpolate_with_Span
                  (file,full_output,sps,stoptol,membtol,deco,deco_last);
      when 3 => Standard_Incremental_Central_Interpolate
                  (file,full_output,sps,stoptol,membtol,deco,deco_last);
      when others => null;
    end case;
    Sampling_Machine.Clear;
    tstop(timer);
    put_line(file,"The list of irreducible components : ");
    put(file,deco);
    new_line(file);
    Write_Summary(file,dim,deco);
    new_line(file);
    print_times(file,timer,"Numerical Irreducible Decomposition");
  end Call_Standard_Interpolate;

  procedure Call_Multprec_Interpolate
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim : in natural32; itp : in natural32 ) is

    use Standard_Complex_Solutions;

    timer : Timing_Widget;
    deco,deco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    size,deci : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    full_output : boolean := false;
    stoptol,membtol : double_float;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,dim);
    stoptol := 10.0**integer(-deci/2);
    membtol := stoptol;
    Multprec_Tuner(file,full_output,size,stoptol,membtol);
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    new_line(file);
    put_line(file,"Diagnostics of filtering process :");
    Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
    case itp is
      when 0 => Multprec_Massive_Interpolate
                  (file,full_output,sps,size,stoptol,membtol,deco,deco_last);
      when 1 => Multprec_Incremental_Interpolate
                  (file,full_output,sps,size,stoptol,membtol,deco,deco_last);
      when 2 => Multprec_Incremental_Interpolate_with_Span
                  (file,full_output,sps,size,stoptol,membtol,deco,deco_last);
      when 3 => Multprec_Incremental_Central_Interpolate
                  (file,full_output,sps,size,stoptol,membtol,deco,deco_last);
      when others => null;
    end case;
    Sampling_Machine.Clear;
    tstop(timer);
    put_line(file,"The list of irreducible components : ");
    put(file,deco);
    new_line(file);
    Write_Summary(file,dim,deco);
    new_line(file);
    print_times(file,timer,"Numerical Irreducible Decomposition");
  end Call_Multprec_Interpolate;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    Call_Monodromy_Breakup(file,p,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    Call_Monodromy_Breakup(file,p,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    Call_Monodromy_Breakup(file,p,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    Call_Monodromy_Breakup(file,p,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    Call_Monodromy_Breakup(file,p,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    f : Standard_Natural_VecVecs.Link_to_VecVec;

  begin
    Call_Monodromy_Breakup(file,p,sols,dim,f);
    Standard_Natural_VecVecs.Deep_Clear(f);
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    nbl : constant natural32 := 20;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    Sampling_Machine.Initialize(p);
   -- Sampling_Machine.Default_Tune_Sampler(0);
    new_line;
    Sampling_Machine.Interactive_Tune_Sampler(file);
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    grid := Create(file,p,sols,dim);
    Factor(file,dim,nbl,grid,f);
    Sampling_Machine.Clear;
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    nbl : constant natural32 := 20;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    DoblDobl_Sampling_Machine.Initialize(p);
   -- DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    new_line;
    DoblDobl_Sampling_Machine.Interactive_Tune_Sampler(file);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    grid := Create(file,p,sols,dim);
    Factor(file,dim,nbl,grid,f);
    DoblDobl_Sampling_Machine.Clear;
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    nbl : constant natural32 := 20;
    grid : Array_of_QuadDobl_Sample_Lists(0..2);

  begin
    QuadDobl_Sampling_Machine.Initialize(p);
   -- DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    new_line;
    QuadDobl_Sampling_Machine.Interactive_Tune_Sampler(file);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    grid := Create(file,p,sols,dim);
    Factor(file,dim,nbl,grid,f);
    QuadDobl_Sampling_Machine.Clear;
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in Standard_Complex_Laur_Systems.Laur_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    nbl : constant natural32 := 20;
    grid : Array_of_Standard_Sample_Lists(0..2);

  begin
    Sampling_Laurent_Machine.Initialize(p);
   -- Sampling_Machine.Default_Tune_Sampler(0);
    new_line;
    Sampling_Laurent_Machine.Interactive_Tune_Sampler(file);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    grid := Create(file,p,sols,dim);
    Laurent_Factor(file,dim,nbl,grid,f);
    Sampling_Laurent_Machine.Clear;
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    nbl : constant natural32 := 20;
    grid : Array_of_DoblDobl_Sample_Lists(0..2);

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(p);
   -- DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    new_line;
    DoblDobl_Sampling_Laurent_Machine.Interactive_Tune_Sampler(file);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    grid := Create(file,p,sols,dim);
    Laurent_Factor(file,dim,nbl,grid,f);
    DoblDobl_Sampling_Laurent_Machine.Clear;
  end Call_Monodromy_Breakup;

  procedure Call_Monodromy_Breakup
              ( file : in file_type;
                p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32;
                f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    grid : Array_of_QuadDobl_Sample_Lists(0..2);
    nbl : constant natural32 := 20;

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(p);
   -- DoblDobl_Sampling_Machine.Default_Tune_Sampler(0);
    new_line;
    QuadDobl_Sampling_Laurent_Machine.Interactive_Tune_Sampler(file);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    grid := Create(file,p,sols,dim);
    Laurent_Factor(file,dim,nbl,grid,f);
    QuadDobl_Sampling_Laurent_Machine.Clear;
  end Call_Monodromy_Breakup;

  function Last_Element ( sps : Standard_Sample_List )
                        return Standard_Sample_List is

  -- DESCRIPTION :
  --   Returns a list whose first element is the last element in sps.

    res,tmp : Standard_Sample_List := sps;

  begin
    if not Is_Null(sps) then
      tmp := Tail_Of(sps);
      while not Is_Null(tmp) loop
        res := Tail_Of(res);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
    return res;
  end Last_Element;

  procedure Split
              ( deco : in Standard_Irreducible_Component_List;
                d : in natural32;
                ld,ld_last : in out Standard_Irreducible_Component_List;
                md,md_last : in out Multprec_Irreducible_Component_List ) is

  -- DESCRIPTION :
  --   The list deco is split into two lists: those components of degree
  --   less or equal than d go into ld, the others into md.

    tmp : Standard_Irreducible_Component_List := deco;

  begin
    while not Is_Null(tmp) loop
      declare
        stcp : constant Standard_Irreducible_Component := Head_Of(tmp); 
        mpcp : Multprec_Irreducible_Component;
        sps,sps_last : Standard_Sample_List;
      begin
        if Degree(stcp) <= d then
          Append(ld,ld_last,stcp);
        else
          Initialize_Labels(mpcp,Labels(stcp).all);
          sps := Points(stcp);
          sps_last := Last_Element(sps);
          Add_Points(mpcp,sps,sps_last);
          Append(md,md_last,mpcp);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Split;

  procedure Write_Numerical_Results
              ( file : in file_type; numres : in Matrix ) is

  -- DESCRIPTION :
  --   Writes the table with numerical results from the Newton interpolation.

  begin
    put_line(file,"-------------------------------------------------------");
    put(file,"|   d |");
    put(file,"    eps    |");
    put(file," distance  |");
    put(file,"  gridres  |");
    put(file,"  testres  |"); new_line(file);
    put_line(file,"-------------------------------------------------------");
    for i in numres'range(1) loop
      put(file,"| ");
      put(file,integer32(numres(i,1)),3); put(file," |");
      for j in 2..integer32(5) loop
        put(file,numres(i,j),2,3,3); put(file," |");
      end loop;
      new_line(file);
    end loop;
    put_line(file,"-------------------------------------------------------");
  end Write_Numerical_Results;

  procedure Call_Newton_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    stdeco,stdeco_last : Standard_Irreducible_Component_List;
    mpdeco,mpdeco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size,thresd : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    timer : Timing_Widget;

  begin
    new_line;
    put_line
      ("Reading the name of the file for the labels in the decomposition");
    Read_Name_and_Open_File(decofile);
    get_labels(decofile,deco,deco_last);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(p,mp,size,dim);
      new_line;
      Write_Summary(Standard_Output,dim,deco);
      new_line;
      put("What is the threshold degree for standard arithmetic ? ");
      get(thresd);
      Distribute_Points(deco,sps);
      Split(deco,thresd,stdeco,stdeco_last,mpdeco,mpdeco_last);
      Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
      Sampling_Machine.Default_Tune_Sampler(2);
      Sampling_Machine.Default_Tune_Refiner;
      Sampling_Machine.Interactive_Tune_Refiner(size);
      new_line;
      put_line("See the output file for results...");
      new_line;
      tstart(timer);
      if not Is_Null(stdeco) then
        put(file,"Newton Interpolation with Standard Arithmetic ");
        put(file,"for degrees <= "); put(file,thresd,1);
        put_line(file," :");
        declare
          numres : Matrix(1..integer32(Length_Of(stdeco)),1..5);
        begin
          Standard_Newton_Interpolate(file,stdeco,numres);
          Write_Numerical_Results(file,numres);
        end;
      end if;
      if not Is_Null(mpdeco) then
        put(file,"Newton Interpolation with ");
        put(file,deci,1); put(file," Decimal Places ");
        put(file,"for Degrees > "); put(file,thresd,1);
        put_line(file," :");
        declare
          numres : Matrix(1..integer32(Length_Of(mpdeco)),1..5);
        begin
          Multprec_Newton_Interpolate(file,size,mpdeco,numres);
          Write_Numerical_Results(file,numres);
        end;
      end if;
    else
      new_line;
      Write_Summary(Standard_Output,dim,deco);
      new_line;
      put_line("See the output file for results...");
      new_line;
      Distribute_Points(deco,sps);
      Sampling_Machine.Initialize(p);
      Sampling_Machine.Default_Tune_Sampler(2);
      Sampling_Machine.Default_Tune_Refiner;
      put_line(file,"Newton Interpolation with Standard Arithmetic :");
      tstart(timer);
      declare
        numres : Matrix(1..integer32(Length_Of(deco)),1..5);
      begin
        Standard_Newton_Interpolate(file,deco,numres);
        Write_Numerical_Results(file,numres);
      end;
    end if;
    tstop(timer);
    Sampling_Machine.Clear;
    new_line(file);
    print_times(file,timer,"Newton Interpolation with Divided Differences");
  end Call_Newton_Interpolate;

  procedure Call_Trace_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    stdeco,stdeco_last : Standard_Irreducible_Component_List;
    mpdeco,mpdeco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size,thresd : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    timer : Timing_Widget;

  begin
    new_line;
    put_line
      ("Reading the name of the file for the labels in the decomposition");
    Read_Name_and_Open_File(decofile);
    get_labels(decofile,deco,deco_last);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(p,mp,size,dim);
      new_line;
      Write_Summary(Standard_Output,dim,deco);
      new_line;
      put("What is the threshold degree for standard arithmetic ? ");
      get(thresd);
      Distribute_Points(deco,sps);
      Split(deco,thresd,stdeco,stdeco_last,mpdeco,mpdeco_last);
      Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
      Sampling_Machine.Default_Tune_Sampler(2);
      Sampling_Machine.Default_Tune_Refiner;
      Sampling_Machine.Interactive_Tune_Refiner(size);
      new_line;
      put_line("See the output file for results...");
      new_line;
      tstart(timer);
      if not Is_Null(stdeco) then
        put(file,"Traces Interpolation with Standard Arithmetic ");
        put(file,"for Degrees <= "); put(file,thresd,1);
        put_line(file," :");
        declare
          numres : Matrix(1..integer32(Length_Of(stdeco)),1..5);
        begin
          Standard_Trace_Interpolate(file,stdeco,numres);
          Write_Numerical_Results(file,numres);
        end;
      end if;
      if not Is_Null(mpdeco) then
        put(file,"Traces Interpolation with ");
        put(file,deci,1); put(file," Decimal Places ");
        put(file,"for degrees > "); put(file,thresd,1);
        put_line(file," :");
        declare
          numres : Matrix(1..integer32(Length_Of(mpdeco)),1..5);
        begin
          Multprec_Trace_Interpolate(file,size,mpdeco,numres);
          Write_Numerical_Results(file,numres);
        end;
      end if;
    else
      new_line;
      Write_Summary(Standard_Output,dim,deco);
      new_line;
      put_line("See the output file for results...");
      new_line;
      Distribute_Points(deco,sps);
      Sampling_Machine.Initialize(p);
      Sampling_Machine.Default_Tune_Sampler(2);
      Sampling_Machine.Default_Tune_Refiner;
      put_line(file,"Traces Interpolation with Standard Arithmetic :");
      tstart(timer);
      declare
        numres : Matrix(1..integer32(Length_Of(deco)),1..5);
      begin
        Standard_Trace_Interpolate(file,deco,numres);
        Write_Numerical_Results(file,numres);
      end;
    end if;
    tstop(timer);
    Sampling_Machine.Clear;
    new_line(file);
    print_times(file,timer,"Interpolation with Traces");
  end Call_Trace_Interpolate;

  procedure Call_Power_Trace_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    stdeco,stdeco_last : Standard_Irreducible_Component_List;
    mpdeco,mpdeco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size,thresd : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    timer : Timing_Widget;

  begin
    new_line;
    put_line
      ("Reading the name of the file for the labels in the decomposition");
    Read_Name_and_Open_File(decofile);
    get_labels(decofile,deco,deco_last);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(p,mp,size,dim);
      new_line;
      Write_Summary(Standard_Output,dim,deco);
      new_line;
      put("What is the threshold degree for standard arithmetic ? ");
      get(thresd);
      Distribute_Points(deco,sps);
      Split(deco,thresd,stdeco,stdeco_last,mpdeco,mpdeco_last);
      Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
      Sampling_Machine.Default_Tune_Sampler(2);
      Sampling_Machine.Default_Tune_Refiner;
      Sampling_Machine.Interactive_Tune_Refiner(size);
      new_line;
      put_line("See the output file for results...");
      new_line;
      tstart(timer);
      if not Is_Null(stdeco) then
        put(file,"Power Trace Interpolation with Standard Arithmetic ");
        put(file,"for Degrees <= "); put(file,thresd,1);
        put_line(file," :");
        declare
          numres : Matrix(1..integer32(Length_Of(stdeco)),1..5);
        begin
          Standard_Power_Trace_Interpolate(file,stdeco,numres);
          Write_Numerical_Results(file,numres);
        end;
      end if;
      if not Is_Null(mpdeco) then
        put(file,"Power Trace Interpolation with ");
        put(file,deci,1); put(file," Decimal Places ");
        put(file,"for Degrees > "); put(file,thresd,1);
        put_line(file," :");
        declare
          numres : Matrix(1..integer32(Length_Of(mpdeco)),1..5);
        begin
          Multprec_Power_Trace_Interpolate(file,size,mpdeco,numres);
          Write_Numerical_Results(file,numres);
        end;
      end if;
    else
      new_line;
      Write_Summary(Standard_Output,dim,deco);
      new_line;
      put_line("See the output file for results...");
      new_line;
      Distribute_Points(deco,sps);
      Sampling_Machine.Initialize(p);
      Sampling_Machine.Default_Tune_Sampler(2);
      Sampling_Machine.Default_Tune_Refiner;
      put_line(file,"Power Trace Interpolation with Standard Arithmetic :");
      tstart(timer);
      declare
        numres : Matrix(1..integer32(Length_Of(deco)),1..5);
      begin
        Standard_Power_Trace_Interpolate(file,deco,numres);
        Write_Numerical_Results(file,numres);
      end;
    end if;
    tstop(timer);
    Sampling_Machine.Clear;
    new_line(file);
    print_times(file,timer,"Interpolation with Power Traces");
  end Call_Power_Trace_Interpolate;

  procedure Call_Linear_Trace_Interpolate
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 ) is

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    hyp : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    timer : Timing_Widget;

  begin
    new_line;
    put_line
      ("Reading the name of the file for the labels in the decomposition");
    Read_Name_and_Open_File(decofile);
    get_labels(decofile,deco,deco_last);
    new_line;
    Write_Summary(Standard_Output,dim,deco);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Distribute_Points(deco,sps);
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    put_line(file,"Linear Traces Interpolation with Standard Arithmetic :");
    tstart(timer);
    declare
      numres : Matrix(1..integer32(Length_Of(deco)),1..5);
    begin
      Standard_Linear_Trace_Interpolate(file,deco,numres);
      Write_Numerical_Results(file,numres);
    end;
    tstop(timer);
    Sampling_Machine.Clear;
    new_line(file);
    print_times(file,timer,"Interpolation with Linear Traces");
  end Call_Linear_Trace_Interpolate;

end Drivers_to_Factor_Components;
