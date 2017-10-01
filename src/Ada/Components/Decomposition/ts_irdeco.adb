with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;
with Standard_Natural_Matrices;
with Standard_Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Sampling_Machine;
with Irreducible_Component_Lists;        use Irreducible_Component_Lists;
with Irreducible_Decompositions;         use Irreducible_Decompositions;
with Irreducible_Decompositions_io;      use Irreducible_Decompositions_io;
with Witness_Generate_and_Classify;      use Witness_Generate_and_Classify;
with Drivers_to_Component_Creators;      use Drivers_to_Component_Creators;

procedure ts_irdeco is

  function Column_Sum ( tab : Standard_Natural_Matrices.Matrix;
                        i : integer32 ) return natural32 is

  -- DESCRIPTION :
  --   Returns the sum of all entries in the i-th column.

    res : natural32 := 0;

  begin
    for j in tab'range loop
      res := res + tab(j,i);
    end loop;
    return res;
  end Column_Sum;

  procedure Write_Generate_Summary
              ( file : in file_type; timings : in Array_of_Duration;
                flowtab : in Standard_Natural_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the summary for the sequence of homotopy continuation stages
  --   executed in Witness_Generate.

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b00 : constant string :=
     "  =====================================================================";
    b1 : constant string :=
     "  |       TIMING INFORMATION SUMMARY for Cascade of Homotopies        |";
    b2 : constant string :=
     "  | level | #paths | #nsols | #comps | #infty |     user cpu time     |";

    total : duration := timings(timings'first);

  begin
    for i in timings'first+1..timings'last loop
      total := total + timings(i);
    end loop;
    new_line(file);
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b00);
    for i in reverse flowtab'range(1) loop
      put(file,"  | ");
      put(file,i,4); put(file,"  | ");
      put(file,flowtab(i,1),5); put(file,"  | ");
      put(file,flowtab(i,2),5); put(file,"  | ");
      put(file,flowtab(i,3),5); put(file,"  | ");
      put(file,flowtab(i,4),5); put(file,"  |    ");
      print_hms(file,timings(integer(i))); put_line(file,"     |");
    end loop;
    put_line(file,b0);
    put(file,"  | total | ");
    put(file,Column_Sum(flowtab,1),5);  put(file,"  | ");
    put(file,Column_Sum(flowtab,2),5);  put(file,"  | ");
    put(file,Column_Sum(flowtab,3),5);  put(file,"  | ");
    put(file,Column_Sum(flowtab,4),5); put(file,"  |    ");
    print_hms(file,total); put_line(file,"     |");
    put_line(file,b0);
  end Write_Generate_Summary;

  procedure Write ( file : in file_type; fp : in List; i : in integer32 ) is

  -- DESCRIPTION :
  --   Writes the i-th component of every element in the list.

    tmp : List := fp;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      put(file,lpt(i),4); put(file," |");
      tmp := Tail_Of(tmp);
    end loop;
  end Write;

  procedure Write_Banner ( file : in file_type; 
                           m : in natural32; sep : character ) is
  begin
    for i in 1..m loop
      for j in 1..6 loop
        put(file,sep);
      end loop;
    end loop;
  end Write_Banner;

  procedure Write_Classify_Summary
              ( file : in file_type;
                timings : in Array_of_Duration; fp : in List ) is

    n : integer32;
    m : constant natural32 := Length_Of(fp);
    total : duration := timings(timings'first);

  begin
    for i in timings'first+1..timings'last loop
      total := total + timings(i);
    end loop;
    new_line(file);
    put_line(file,"  -------------------------------------------------------");
    put_line(file,"  |  TIMING INFORMATION SUMMARY for Classifying Points  |");
    put_line(file,"  -------------------------------------------------------");
    if m /= 0 then
      n := Head_Of(fp)'last;
      put(file,"  | dimension "); put(file," |");
      Write(file,fp,n); 
      put_line(file,"  user cpu time  |");
      put(file,"  =============="); Write_Banner(file,m,'=');
      put_line(file,"==================");
      if n > 0 then
        put(file,"  | level ");
        put(file,n-1,1);  put(file,"    |");
        Write(file,fp,n-1);
        put(file," "); print_hms(file,timings(integer(n)-1));
        put_line(file,"  |");
        for i in reverse 0..n-2 loop
          put(file,"  |       ");
          put(file,i,1);  put(file,"    |");
          Write(file,fp,i);
          put(file," "); print_hms(file,timings(integer(i)));
          put_line(file,"  |");
        end loop;
      end if;
      put(file,"  --------------"); Write_Banner(file,m,'-');
      put_line(file,"------------------");
    end if;
    put(file,"  | total time : ");
    Write_Banner(file,m,' ');
    print_hms(file,total); put_line(file,"  |");
    put(file,"  --------------"); Write_Banner(file,m,'-');
    put_line(file,"------------------");
  end Write_Classify_Summary;

  procedure Write_Summary ( file : in file_type;
                            dc : in Standard_Irreducible_Decomposition ) is

    cpl : Standard_Irreducible_Component_List;
    isol : natural32;

  begin
    for i in reverse 1..Top_Dimension(dc) loop
      cpl := Components(dc,i);
      Write_Summary(file,natural32(i),cpl);
    end loop;
    isol := Length_Of(Generic_Points(dc,0));
    put(file,"Found "); put(file,isol,1); 
    if isol = 1
     then put_line(file," isolated solution.");
     else put_line(file," isolated solutions.");
    end if;
  end Write_Summary;

  procedure Write_Summary ( file : in file_type;
                            dc : in Multprec_Irreducible_Decomposition ) is

    cpl : Multprec_Irreducible_Component_List;
    isol : natural32;

  begin
    for i in reverse 1..Top_Dimension(dc) loop
      cpl := Components(dc,i);
      Write_Summary(file,natural32(i),cpl);
    end loop;
    isol := Length_Of(Generic_Points(dc,0));
    put(file,"Found "); put(file,isol,1); 
    if isol = 1
     then put_line(file," isolated solution.");
     else put_line(file," isolated solutions.");
    end if;
  end Write_Summary;

  procedure Test_io ( file : in file_type;
                      p : in Poly_Sys; k : in integer32 ) is

  -- DESCRIPTION :
  --   Tests on the creation and input-output of decompositions.

    ep : constant Array_of_Poly_Sys(0..k) := Slice_and_Embed(p,natural32(k));
    dc : constant Standard_Irreducible_Decomposition := Create(ep);
    infile,outfile : file_type;
    dc_read : Standard_Irreducible_Decomposition;

  begin
    Add_Embed_Symbols(natural32(k));
    put_line(file,"The irreducible decomposition : ");
    put(file,dc);
    new_line;
    put_line("Reading the file name to write the decomposition to.");
    Read_Name_and_Create_File(outfile);
    put(outfile,dc);
    Close(outfile);
    new_line;
    put_line("Reading the name of the file to retrieve the decomposition.");
    Read_Name_and_Open_File(infile);
    get(infile,dc_read);
    Close(infile);
    put_line(file,"The irreducible decomposition read from file : ");
    put(file,dc_read);
  end Test_io;

  procedure Test_Create ( file : in file_type;
                          p : in Poly_Sys; k : in integer32 ) is

  -- DESCRIPTION :
  --   Tests on the creation and input-output of decompositions.

    timer : Timing_Widget;
    dc : Standard_Irreducible_Decomposition;
    degroco : boolean;
    tol : constant double_float := 1.0E-12;
    timings : Array_of_Duration(0..integer(k));
    flowtab : Standard_Natural_Matrices.Matrix(0..k,1..4);
    ans : character;

  begin
    Add_Embed_Symbols(natural32(k));
    new_line;
    put("Full root counting or only based on degrees? (f/d) ");
    Ask_Alternative(ans,"fd");
    degroco := (ans = 'd');
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    Witness_Generate(file,p,k,degroco,tol,timings,flowtab,dc);
    tstop(timer);
    Write_Generate_Summary(file,timings,flowtab);
    put_line(file,"THE IRREDUCIBLE DECOMPOSITION :");
    put(file,dc);
    new_line(file);
    print_times(file,timer,"Witness Generate");
  end Test_Create;

  procedure Display_Breakup_Menu ( method : out natural32 ) is

  -- DESCRIPTION :
  --   Display the menu with the different interpolation strategies
  --   to break up solution sets into irreducible components.

    ans : character;

  begin
    new_line;
    put_line("MENU to break up solution sets into irreducible components : ");
    put_line("  0. massive interpolate, only compares polynomials at end;");
    put_line("  1. incremental interpolate with standard linear projections;");
    put_line("  2. interpolate with exploitation of span of components;");
    put_line("  3. interpolate with central projections;");
    put_line("  4. use monodromy group actions and homotopy filters.");
    put("Type 0,1,2,3 or 4 to select method : ");
    Ask_Alternative(ans,"01234");
    case ans is
      when '0' => method := 0;
      when '1' => method := 1;
      when '2' => method := 2;
      when '3' => method := 3;
      when '4' => method := 4;
      when others => method := 4;
    end case;
  end Display_Breakup_Menu;

  procedure Test_Classify ( file : in file_type;
                            dc : in out Standard_Irreducible_Decomposition ) is

    timer : Timing_Widget;
    method : natural32 := 0;
    full_output : boolean := false;
    stoptol,membtol : double_float;
    k : constant integer32 := Top_Dimension(dc);
    clatims : Array_of_Duration(0..integer(k)) := (0..integer(k) => 0.0);
    fp,fp_last : List;
    deci,size,threshold : natural32 := 0;
    p : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Embedding(dc,k);
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    mdc : Multprec_Irreducible_Decomposition;

  begin
    Display_Breakup_Menu(method);
    if method < 4 then
      new_line;
      put("Give the number of decimal places (16 = standard) : ");
      get(deci); 
      if deci <= 16 then
        size := 0;
        stoptol := 1.0E-8;
        membtol := 1.0E-6;
        Standard_Tuner(file,full_output,stoptol,membtol);
      else
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        Get_Multprec_System(p.all,mp,size,natural32(k));
        mdc := Create(dc);
        Add_Original(mdc,mp);
        stoptol := 10.0**integer(-deci/2);
        membtol := stoptol;
        Multprec_Tuner(file,full_output,size,stoptol,membtol);
      end if;
    else
      size := 0;
      new_line;
      put("Give stabilizing threshold : "); 
      get(threshold); skip_line;
      stoptol := double_float(threshold);
      membtol := 1.0E-6;
      new_line;
      Sampling_Machine.Interactive_Tune_Sampler(file);
      Sampling_Machine.Interactive_Tune_Refiner(file);
    end if;
    put_line(file,"THE INPUT IRREDUCIBLE DECOMPOSITION :");
    put(file,dc);
    new_line;
    put_line("See the output file for results...");
    new_line;
    put_line(file,"Diagnostics of Breakup and Filter :");
    tstart(timer);
    if size = 0 then
       Witness_Classify
         (file,full_output,dc,method,stoptol,membtol,clatims,fp,fp_last);
    else
       Witness_Classify
         (file,full_output,mdc,method,size,stoptol,membtol,clatims,fp,fp_last);
    end if;
    tstop(timer);
    new_line(file);
    put_line(file,"THE OUTPUT IRREDUCIBLE DECOMPOSITION :");
    put(file,dc);
    new_line(file);
    Write_Classify_Summary(file,clatims,fp);
    new_line(file);
    if size = 0
     then Write_Summary(file,dc);
     else Write_Summary(file,mdc);
    end if;
    new_line(file);
    print_times(file,timer,"Witness Classify");
  end Test_Classify;

  procedure Test_Generate_and_Classify
                 ( file : in file_type;
                   p : in Poly_Sys; k : in integer32 ) is

  -- DESCRIPTION :
  --   Tests on the creation and input-output of decompositions.

    timer : Timing_Widget;
    sdc : Standard_Irreducible_Decomposition;
    mdc : Multprec_Irreducible_Decomposition;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    degroco : boolean;
    full_output : boolean := false;
    zerotol : constant double_float := 1.0E-10;
    membtol,stoptol : double_float;
    gentims,clatims
      : Array_of_Duration(0..integer(k)) := (0..integer(k) => 0.0);
    flowtab : Standard_Natural_Matrices.Matrix(0..k,1..4);
    method,deci,size : natural32 := 0;
    ans : character;
    fp,fp_last : List;

  begin
    new_line;
    put("Full root counting or only based on degrees? (f/d) ");
    Ask_Alternative(ans,"fd");
    degroco := (ans = 'd');
    Display_Breakup_Menu(method);
    new_line;
    put("Give the number of decimal places (16 = standard) : "); get(deci);
    if deci <= 16 then
      size := 0;
      stoptol := 1.0E-8;
      membtol := 1.0E-6;
      Standard_Tuner(file,full_output,stoptol,membtol);
    else
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(p,mp,size,0); 
      stoptol := 10.0**integer(-deci/2);
      membtol := stoptol;
      Multprec_Tuner(file,full_output,size,stoptol,membtol);
    end if;
    Add_Embed_Symbols(natural32(k));
    new_line;
    put_line("See the output file for results...");
    new_line;
    tstart(timer);
    if size = 0 then
      Witness_Generate_Classify
        (file,full_output,p,k,method,degroco,zerotol,stoptol,membtol,
         gentims,clatims,flowtab,sdc,fp,fp_last);
    else
      Witness_Generate_Classify
        (file,full_output,p,mp,k,method,size,degroco,zerotol,stoptol,
         membtol,gentims,clatims,flowtab,mdc,fp,fp_last);
    end if;
    tstop(timer);
    if size = 0 then
      put_line(file,"THE IRREDUCIBLE DECOMPOSITION :");
      put(file,sdc);
    end if;
    Write_Generate_Summary(file,gentims,flowtab);
    new_line(file);
    Write_Classify_Summary(file,clatims,fp);
    new_line(file);
    if size = 0
     then Write_Summary(file,sdc);
     else Write_Summary(file,mdc);
    end if;
    new_line(file);
    print_times(file,timer,"Witness Generate and Classify");
  end Test_Generate_and_Classify;

  procedure Main is

    ans : character;
    infile,outfile : file_type;
    lp : Link_to_Poly_Sys;
    k : integer32 := 0;
    dc : Standard_Irreducible_Decomposition;
 
  begin
    new_line;
    put_line("Irreducible decomposition of solution set of polynomial system.");
    new_line;
    put_line("MENU to test irreducible decomposition : ");
    put_line("  1. test input-output facilities;");
    put_line("  2. creation of irreducible decomposition;");
    put_line("  3. classify given candidate lists of generic points;");
    put_line("  4. creation + classification of generic points.");
    put("Type 1,2,3 or 4 to select : "); Ask_Alternative(ans,"1234");
    new_line;
    if ans = '3' then
      put_line("Reading name of the file with candidate generic points.");
      Read_Name_and_Open_File(infile);
      get(infile,dc);
    else
      put_line("Reading the name of the file for the input system.");
      Read_Name_and_Open_File(infile);
      get(infile,lp);
    end if;
    Close(infile);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    if ans /= '3' then
      new_line;
      put("Give the top dimension : "); get(k); skip_line;
    end if;
    case ans is
      when '1' => Test_io(outfile,lp.all,k);
      when '2' => Test_Create(outfile,lp.all,k);
      when '3' => Test_Classify(outfile,dc);
      when '4' => Test_Generate_and_Classify(outfile,lp.all,k);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_irdeco;
