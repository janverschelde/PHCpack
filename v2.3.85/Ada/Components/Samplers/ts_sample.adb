with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;        use Standard_Complex_VecVecs_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Sample_Points,Sample_Points_io;     use Sample_Points,Sample_Points_io;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Sample_Point_Lists_io;              use Sample_Point_Lists_io;
with Sample_Point_Grids;                 use Sample_Point_Grids;
with Sample_Point_Grids_io;              use Sample_Point_Grids_io;
with Sampling_Machine;

procedure ts_sample is

-- DESCRIPTION :
--   Interactive testing of the sampling module.

  procedure Check_Sample ( file : in file_type; spt : Standard_Sample ) is

    n : constant integer32 := Number_of_Variables(spt);
    sol : Standard_Complex_Solutions.Solution(n) := Sample_Point(spt);

  begin
    put(file,"  err : ");
    put(file,sol.err,3);
    put(file,"  rco : ");
    put(file,sol.rco,3);
    put(file,"  res : ");
    put(file,sol.res,3);
    new_line(file);
  end Check_Sample;

  procedure Check_Sample ( file : in file_type; spt : Multprec_Sample ) is

    n : constant integer32 := Number_of_Variables(spt);
    sol : constant Multprec_Complex_Solutions.Solution(n) := Sample_Point(spt);

  begin
    put(file,"  err : ");
    put(file,sol.err,3);
    put(file,"  rco : ");
    put(file,sol.rco,3);
    put(file,"  res : ");
    put(file,sol.res,3);
    new_line(file);
  end Check_Sample;

  procedure Write_Diagnostics
               ( file : in file_type; samples : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Writes diagnostics for each sample on file.

    tmp : Standard_Sample_List := samples;

  begin
    while not Is_Null(tmp) loop
      Check_Sample(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Diagnostics;

  procedure Write_Diagnostics
               ( file : in file_type; samples : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Writes diagnostics for each sample on file.

    tmp : Multprec_Sample_List := samples;

  begin
    while not Is_Null(tmp) loop
      Check_Sample(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Diagnostics;

  procedure Test_List_io ( samples : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the list of samples to a file and the reads it back in.

    infile,outfile : file_type;
    n,k : integer32;
    len : natural32;
    first,last : Standard_Sample_List;

  begin
    put_line("The list of samples : "); put(samples);
    if not Is_Null(samples) then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      len := Length_Of(samples);
      n := Number_of_Variables(Head_Of(samples));
      k := Number_of_Slices(Head_Of(samples));
      put(outfile,len,natural32(n),natural32(k),samples);
      Close(outfile);
      new_line;
      put_line("Reading the name of the file for the samples.");
      Read_Name_and_Open_File(infile);
      get(infile,first,last);
      put_line("The list of samples from file : "); put(first);
    end if;
  end Test_List_io;

  procedure Test_Grid_io ( samples : in Multprec_Sample_Grid ) is

  -- DESCRIPTION :
  --   Writes the grid of samples to a file and the reads it back in.

    infile,outfile : file_type;
    first,last : Multprec_Sample_Grid;

  begin
    put_line("The grid of samples : "); put(samples);
    if not Is_Null(samples) then
      new_line;
      put_line("Reading the name of the output file.");
      Read_Name_and_Create_File(outfile);
      put(outfile,samples);
      Close(outfile);
      new_line;
      put_line("Reading the name of the file for the samples.");
      Read_Name_and_Open_File(infile);
      get(infile,first,last);
      put_line("The grid of samples from file : "); put(first);
    end if;
  end Test_Grid_io;

  procedure Test_Create ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                          sols : in Standard_Complex_Solutions.Solution_List;
                          dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates a list of samples from the slices added to p and the
  --   solutions in the list sols.  Performs a test on input/output.

    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    samples,samples_last : Standard_Sample_List;

  begin
    put_line("The slices added in the embedding : ");
    put_line(hyp);
    for i in 1..Length_of(sols) loop
      put("Sample "); put(i,1); put_line(" :");
      declare
        sample : constant Standard_Sample := Create(Head_Of(tmp).all,hyp);
      begin
        put(sample);
        Check_Sample(Standard_Output,sample);
        Append(samples,samples_last,sample);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Test_List_io(samples);
    Deep_Clear(samples);
  end Test_Create; 

  procedure Test_Sample ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                          sols : in Standard_Complex_Solutions.Solution_list;
                          dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates one new sample.

    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    s1 : constant Standard_Sample := Create(Head_Of(sols).all,hyp);
    ans : character;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner;
    put_line("The current sample : "); put(s1);
    loop
      declare
        s2 : Standard_Sample;
      begin
        Sample(s1,s2);
        put_line("A new sample : "); put(s2);
        Check_Sample(Standard_output,s2);
        Deep_Clear(s2);
      end;
      put("Do you want more samples ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Sampling_Machine.Clear;
  end Test_Sample;

  procedure Test_Samples ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                           sols : in Standard_Complex_Solutions.Solution_list;
                           dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates a whole new list of samples.

    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    newhyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
           := Random_Hyperplanes(natural32(dim),natural32(p'last));
    use Standard_Complex_Solutions;
    s1 : constant Standard_Sample_List := Create(sols,hyp);
    s2,s2_last : Standard_Sample_List;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner;
    Sample(s1,newhyp,s2,s2_last);
    put_line("The new samples : "); put(s2);
    Deep_Clear(s2);
    Sampling_Machine.Clear;
  end Test_Samples;

  procedure Test_Refine ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                          sols : in Standard_Complex_Solutions.Solution_List;
                          dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates a refinement of a given sample.

    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    s1 : Standard_Sample := Create(Head_Of(sols).all,hyp);
    s2 : Multprec_Sample;
    ans : character;

  begin
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,natural32(dim));
    Sampling_Machine.Initialize(p,mp.all,dim,size);
    Sampling_Machine.Interactive_Tune_Refiner(size);
    Refine(s1,s2);
    loop
      put_line("The refined sample : "); put(s2);
      put("Do you wish to refine this sample more ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Sampling_Machine.Interactive_Tune_Refiner(size);
      Refine(s2);
    end loop;
    Deep_Clear(s2);
    Sampling_Machine.Clear;
  end Test_Refine;

  procedure Test_Sample_and_Refine
                   ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in Standard_Complex_Solutions.Solution_List;
                     dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates a new sample, refined with multi-precision.

    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    s1 : constant Standard_Sample := Create(Head_Of(sols).all,hyp);
    ans : character;

  begin
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,natural32(dim));
    Sampling_Machine.Initialize(p,mp.all,dim,size);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    put_line("The current sample : "); put(s1);
    loop
      declare
        s2 : Multprec_Sample;
      begin
        Sample(s1,s2);
        put_line("A new refined sample : "); put(s2);
        Check_Sample(Standard_Output,s2);
        Deep_Clear(s2);
      end;
      put("Do you wish more samples ? (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    Sampling_Machine.Clear;
  end Test_Sample_and_Refine;

  procedure Test_Multi_Sample 
                ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                  sols : in Standard_Complex_Solutions.Solution_list;
                  dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates a list of samples from one sample.

    timer : Timing_Widget;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    s1 : constant Standard_Sample := Create(Head_Of(sols).all,hyp);
    nb : natural32 := 0;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner;
    put_line("The current sample : "); put(s1);
    loop
      put("Give the number of new samples (0 to exit) : "); get(nb);
      exit when (nb = 0);
      declare
        s2,s2_last : Standard_Sample_List;
      begin
        tstart(timer);
        Sample(s1,nb,s2,s2_last);
        tstop(timer);
        put_line("The new samples : "); put(s2);
        Write_Diagnostics(Standard_Output,s2);
        put("The Elapsed User CPU Time : ");
        print_hms(Standard_Output,Elapsed_User_Time(timer)); new_line;
        Deep_Clear(s2);
      end;
    end loop;
    Sampling_Machine.Clear;
  end Test_Multi_Sample;

  procedure Test_Multi_Refine
                 ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                   sols : in Standard_Complex_Solutions.Solution_List;
                   dim : in integer32 ) is

  -- DESCRIPTION :
  --   Refines a given list of samples.

    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    s1 : Standard_Sample_List := Create(sols,hyp);
    s2,s2_last : Multprec_Sample_List;
    ans : character;

  begin
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,natural32(dim));
    Sampling_Machine.Initialize(p,mp.all,dim,size);
    Sampling_Machine.Interactive_Tune_Refiner(size);
    Refine(s1,s2,s2_last);
    loop
      put_line("The refined samples : "); put(s2);
      Write_Diagnostics(Standard_Output,s2);
      put("Do you wish to refine this sample more ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Sampling_Machine.Interactive_Tune_Refiner(size);
      Refine(s2);
    end loop;
    Deep_Clear(s2);
    Sampling_Machine.Clear;
  end Test_Multi_Refine;

  procedure Test_Multi_Sample_and_Refine
                   ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in Standard_Complex_Solutions.Solution_List;
                     dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates many new samples, refined with multi-precision.

    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;
    timer : Timing_Widget;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    s1 : constant Standard_Sample := Create(Head_Of(sols).all,hyp);
    nb : natural32 := 0;

  begin
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,natural32(dim));
    Sampling_Machine.Initialize(p,mp.all,dim,size);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    put_line("The current sample : "); put(s1);
    loop
      put("Give the number of samples (0 to exit) : "); get(nb);
      exit when (nb = 0);
      declare
        s2,s2_last : Multprec_Sample_List;
      begin
        tstart(timer);
        Sample(s1,nb,s2,s2_last);
        tstop(timer);
        put_line("The new refined samples : "); put(s2);
        Write_Diagnostics(Standard_Output,s2);
        put("The Elapsed User CPU Time : ");
        print_hms(Standard_Output,Elapsed_User_Time(timer)); new_line;
        Deep_Clear(s2);
      end;
    end loop;
    Sampling_Machine.Clear;
  end Test_Multi_Sample_and_Refine;

  procedure Test_Create_Grid
                   ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in Standard_Complex_Solutions.Solution_List;
                     dim : in integer32 ) is

  -- DESCRIPTION :
  --   Creates many new samples, refined with multi-precision.

    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    use Standard_Complex_Solutions;
    s1 : constant Standard_Sample := Create(Head_Of(sols).all,hyp);
    nb : natural32 := 0;
    grid,grid_last : Multprec_Sample_Grid;

  begin
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    Get_Multprec_System(p,mp,size,natural32(dim));
    Sampling_Machine.Initialize(p,mp.all,dim,size);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    loop
      put("Give the number of samples (0 to exit) : "); get(nb);
      exit when (nb = 0);
      declare
        s2,s2_last : Multprec_Sample_List;
      begin
        Sample(s1,nb,s2,s2_last);
        put_line("The new refined samples : "); put(s2);
        Append(grid,grid_last,s2);
      end;
    end loop;
    skip_line;
    Test_Grid_io(grid);
    Deep_Clear(grid);
    Sampling_Machine.Clear;
  end Test_Create_Grid;

  procedure Test_Sampling_Membership
                   ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in Standard_Complex_Solutions.Solution_List;
                     dim : in integer32 ) is

  -- DESCRIPTION :
  --   This does some runs on the sampling membership test.

    tol : double_float := 1.0E-8;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    sps : constant Standard_Sample_List := Create(sols,hyp);
    spt : Standard_Sample;
    tmp : Standard_Sample_List;
    ans : character;
    timer,totaltimer : Timing_Widget;
    point : Standard_Complex_Vectors.Vector(p'range);

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner;
    new_line;
    loop
      tmp := sps;
      for i in 1..Length_Of(sps) loop
        put("Sampling from point "); put(i,1); put(" : ");
        spt := Head_Of(tmp);
        declare
          newspt : Standard_Sample;
          isin : natural32;
          s,s_last : Standard_Sample_List;
        begin
          Sample(spt,newspt);
          Membership_Test(sps,Sample_Point(newspt).v,tol,isin,s,s_last);
          put("isin entry : "); put(isin,1); new_line;
          Deep_Clear(s);
          Deep_Clear(newspt);
        end;
        tmp := Tail_Of(tmp);
      end loop;
      put("Do you want to test another time ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    put_line("Testing random points");
    tstart(totaltimer);
    loop
      point := Random_Vector(1,point'last);
      declare
        isin : natural32;
        s,s_last : Standard_Sample_List;
      begin
        tstart(timer);
        Membership_Test(sps,point,tol,isin,s,s_last);
        tstop(timer);
        put("isin entry : "); put(isin,1);
        put("  Elapsed User CPU time : ");
        print_hms(Standard_Output,Elapsed_User_Time(timer)); new_line;
        Deep_Clear(s);
      end;
      put("Do you want to test another random point ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
    tstop(totaltimer);
    put("  Elapsed User CPU time : ");
    print_hms(Standard_Output,Elapsed_User_Time(totaltimer)); new_line;
    Sampling_Machine.Clear;
  end Test_Sampling_Membership;

  procedure Test_Parallel_Sampling
                   ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in Standard_Complex_Solutions.Solution_List;
                     dim : in integer32 ) is

  -- DESCRIPTION :
  --   Allows to generate samples on parallel slices.  There is only one
  --   moving parameter: the "constant" of the last slice.

    solsfile : file_type;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Slices(p,natural32(dim));
    sps : constant Standard_Sample_List := Create(sols,hyp);
    spt,new_spt : Standard_Sample;
    spsols,spsols_last : Standard_Complex_Solutions.Solution_List;
    nb : natural32 := 0;
    step : double_float := 0.01;
    ans : character;

    use Standard_Complex_Solutions;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Interactive_Tune_Refiner;
    spt := Head_Of(sps);
    put("Give number of samples : "); get(nb);
    put("Give the step size : "); get(step);
    for i in 1..nb loop
      put_line("The constant : "); put(hyp(hyp'last)(0)); new_line;
      hyp(hyp'last)(0) := hyp(hyp'last)(0) + Create(step);
     -- put("Give new constant : "); get(hyp(hyp'last)(0));
      put("The new constant : ");  put(hyp(hyp'last)(0)); new_line;
      Sample_on_Slices(spt,hyp,new_spt);
      declare
        sol : constant Solution(p'last+dim) := Sample_Point(new_spt);
      begin
        put_line("New sample :"); put_vector(sol);
        Append(spsols,spsols_last,sol);
      end;
      put("You wish to continue (y/n) "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Sampling_Machine.Change_Slices(hyp);
    end loop;
    Sampling_Machine.Clear;
    put("Collected "); put(Length_Of(spsols),1); put_line(" samples.");
    put("Do you wish to save these samples to separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of file to write the samples.");
      Read_Name_and_Create_File(solsfile);
      put(solsfile,Length_Of(spsols),natural32(p'last+dim),spsols);
      close(solsfile);
    end if;
  end Test_Parallel_Sampling;

  procedure Main is

  -- DESCRIPTION :
  --   Displays the greetings and reads the embedded system.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List; 
    dim : integer32;
    ans : character;

  begin
    new_line;
    put_line("Sampling generic points from components of solutions");
    Standard_Read_Embedding(lp,sols,natural32(dim));
    loop
      new_line;
      put_line("MENU for testing the operations of the sampling module.");
      put_line("  0. exit this program;");
      put_line("  1. creation and input/output;");
      put_line("  2. sampling from one point;");
      put_line("  3. sampling from the list of points;");
      put_line("  4. multi-precision refinement of a sample;");
      put_line("  5. sampling and refinement from one point;");
      put_line("  6. generate many standard samples starting from one point;");
      put_line("  7. many refined samples from a standard samples list;");
      put_line("  8. many samples + refinements from one point;");
      put_line("  9. test creation + io of a grid;");
      put_line("  A. test sampling membership test;");
      put_line("  B. generate samples on parallel slice.");
      put("Type 0,1,2,3,4,5,6,7,8,9,A or B to select : ");
      Ask_Alternative(ans,"0123456789AB");
      exit when ans = '0';
      new_line;
      case ans is
        when '1' => Test_Create(lp.all,sols,dim);
        when '2' => Test_Sample(lp.all,sols,dim);
        when '3' => Test_Samples(lp.all,sols,dim);
        when '4' => Test_Refine(lp.all,sols,dim);
        when '5' => Test_Sample_and_Refine(lp.all,sols,dim);
        when '6' => Test_Multi_Sample(lp.all,sols,dim);
        when '7' => Test_Multi_Refine(lp.all,sols,dim);
        when '8' => Test_Multi_Sample_and_Refine(lp.all,sols,dim);
        when '9' => Test_Create_Grid(lp.all,sols,dim);
        when 'A' => Test_Sampling_Membership(lp.all,sols,dim);
        when 'B' => Test_Parallel_Sampling(lp.all,sols,dim);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_sample;
