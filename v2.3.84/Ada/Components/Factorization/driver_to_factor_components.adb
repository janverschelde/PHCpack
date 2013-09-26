with text_io;                            use text_io;
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
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Sampling_Machine;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Irreducible_Components;             use Irreducible_Components;
with Irreducible_Component_Lists;        use Irreducible_Component_Lists;
with Irreducible_Component_Lists_io;     use Irreducible_Component_Lists_io;
with Drivers_to_Component_Creators;      use Drivers_to_Component_Creators;

procedure Driver_to_Factor_Components
            ( file : in file_type;
              ep : in Standard_Complex_Poly_Systems.Poly_Sys;
              sols : in Standard_Complex_Solutions.Solution_List;
              dim : in natural32 ) is

-- DESCRIPTION :
--   This driver displays two menus (incremental and monodromy) and then
--   calls the selected interpolation method.

  procedure Call_Standard_Interpolate
               ( file : in file_type;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List;
                 dim : in natural32; itp : in natural32 ) is

  -- DESCRIPTION :
  --   The incremental interpolation scheme is a flexible way to determine
  --   the decomposition of an equidimensional solution set into irreducible
  --   components of solutions.  However, the use standard arithmetic limits
  --   the application range.

    timer : Timing_Widget;
    deco,deco_last : Standard_Irreducible_Component_List;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
   -- spt : Standard_Sample := Create(Head_Of(sols).all,hyp);
    sps : constant Standard_Sample_List := Create(sols,hyp);
   -- maxdeg : constant natural32 := Length_Of(sols);
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

    timer : Timing_Widget;
    deco,deco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    size,deci : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
   -- spt : Standard_Sample := Create(Head_Of(sols).all,hyp);
    sps : constant Standard_Sample_List := Create(sols,hyp);
   -- maxdeg : constant natural32 := Length_Of(sols);
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

    timer : Timing_Widget;
    tol : constant double_float := 1.0E-12;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    threshold,nit : natural32 := 0;
    deco,deco_last : Standard_Irreducible_Component_List;

  begin
    new_line;
    put("Give threshold for number of stable iterations : ");
    get(threshold);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    tstart(timer);
    Monodromy_Breakup(file,sps,threshold,tol,deco,deco_last,nit);
    tstop(timer);
    Sampling_Machine.Clear;
    new_line(file);
    Write_Summary(file,dim,deco);
    new_line(file);
    put_line(file,"The labels of the list of components : ");
    put_labels(file,deco);
    new_line(file);
    put(file,"Number of iterations : "); put(file,nit,1); new_line(file);
    new_line(file);
    print_times(file,timer,"Monodromy Group Action Breakup");
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
      for j in 2..5 loop
        put(file,numres(i,integer32(j)),2,3,3); put(file," |");
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

  -- DESCRIPTION :
  --   The application of Newton interpolation requires the samples to
  --   lie on a structured grid.  Therefore we need a decomposition of
  --   the set of generic points.

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    stdeco,stdeco_last : Standard_Irreducible_Component_List;
    mpdeco,mpdeco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size,thresd : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
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

  -- DESCRIPTION :
  --   The use of traces leads to another form of interpolating polynomials
  --   through generic points of a solution component.  Also here we need
  --   a predicted decomposition to set up a structured grid of samples.

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    stdeco,stdeco_last : Standard_Irreducible_Component_List;
    mpdeco,mpdeco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size,thresd : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
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

  -- DESCRIPTION :
  --   The name "power traces" refers to the Newton identities for
  --   transforming values of the elementary symmetric functions into power
  --   sums.  With these identities we need fewer samples to set up the
  --   trace form of the interpolating polynomial.

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    stdeco,stdeco_last : Standard_Irreducible_Component_List;
    mpdeco,mpdeco_last : Multprec_Irreducible_Component_List;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size,thresd : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
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

  -- DESCRIPTION :
  --   The major advantage of working with traces is that linear ones
  --   suffice to perform the validation process.

    decofile : file_type;
    deco,deco_last : Standard_Irreducible_Component_List;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
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

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("MENU to decompose with incremental use of interpolation :");
    put_line("  1. massive interpolate with standard arithmetic;");
    put_line("  2. incremental interpolate with standard arithmetic;");
    put_line("  3.  + determine span with standard arithmetic;");
    put_line("  4.  + use central projections;");
    put_line("  5. massive interpolate with multi-precision arithmetic;");
    put_line("  6. incremental interpolate with multi-precision arithmetic;");
    put_line("  7.  + determine span with multi-precision arithmetic;");
    put_line("  8.  + use central projections;");
    put_line("MENU to decompose with monodromy verified with interpolation :");
    put_line("  9. predict decomposition with monodromy loops;");
    put_line
      ("  A. on given decomposition : use bootstrapping Newton to certify;");
    put_line
      ("  B.                        : use full trace form to certify;");
    put_line
      ("  C.                        : use Newton identities on trace form;");
    put_line
      ("  D.                        : use linear trace only to certify.");
    put("Type 1,2,3,4,5,6,7,8,9,A,B,C or D to make your choice : ");
    Ask_Alternative(ans,"123456789ABCD");
    case ans is
      when '1' => Call_Standard_Interpolate(file,ep,sols,dim,0);
      when '2' => Call_Standard_Interpolate(file,ep,sols,dim,1);
      when '3' => Call_Standard_Interpolate(file,ep,sols,dim,2);
      when '4' => Call_Standard_Interpolate(file,ep,sols,dim,3);
      when '5' => Call_Multprec_Interpolate(file,ep,sols,dim,0);
      when '6' => Call_Multprec_Interpolate(file,ep,sols,dim,1);
      when '7' => Call_Multprec_Interpolate(file,ep,sols,dim,2);
      when '8' => Call_Multprec_Interpolate(file,ep,sols,dim,3);
      when '9' => Call_Monodromy_Breakup(file,ep,sols,dim);
      when 'A' => Call_Newton_Interpolate(file,ep,sols,dim);
      when 'B' => Call_Trace_Interpolate(file,ep,sols,dim);
      when 'C' => Call_Power_Trace_Interpolate(file,ep,sols,dim);
      when 'D' => Call_Linear_Trace_Interpolate(file,ep,sols,dim);
      when others => null;
    end case;
  end Main;

begin
  Main;
end Driver_to_Factor_Components;
