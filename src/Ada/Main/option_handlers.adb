with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Random_Numbers;
with Write_Seed_Number;
with Greeting_Banners;
with Actions_and_Options;
with mainphc,bablphc,bablphc2,bablphc4;
with maingood,mainsymb;
with mainscal,mainred,mainred2,mainred4;
with bablroco,mainroco;
with babldmvc,mainsmvc;
with mainpoco,bablpoco,bablpoco2,bablpoco4;
with mainadep,maintrack;
with mainfac,maindeco;
with mainenum,bablenum,mainfeed;
with mainseries;
with mainvali,bablvali,bablvali2,bablvali4;
with compsolve,compsolve2,compsolve4;
with mainsolve;
with mainhyp,mainhyp2,mainhyp4;
with mainsam,mainwit;
with maindict,mainzip;

package body Option_Handlers is

  fixed_seed : constant natural32 := 68717;

-- BANNERS WITH INFORMATION TO START DIALOGUE WITH USER :

  welcome : constant string := Greeting_Banners.welcome;
  compban : constant string := Greeting_Banners.compban;
  enumban : constant string := Greeting_Banners.enumban;
  facban  : constant string := Greeting_Banners.facban;
  goodban : constant string := Greeting_Banners.goodban;
  symbban : constant string := Greeting_Banners.symbban;
  hypban  : constant string := Greeting_Banners.hypban;
  mvcban  : constant string := Greeting_Banners.mvcban;
  pocoban : constant string := Greeting_Banners.pocoban;
  reduban : constant string := Greeting_Banners.reduban;
  rocoban : constant string := Greeting_Banners.rocoban;
  samban  : constant string := Greeting_Banners.samban;
  scalban : constant string := Greeting_Banners.scalban;
  slvban  : constant string := Greeting_Banners.slvban;
  adepban : constant string := Greeting_Banners.adepban;
  trackban : constant string := Greeting_Banners.trackban;
  seriesban : constant string := Greeting_Banners.seriesban;
  veriban : constant string := Greeting_Banners.veriban;
  witban  : constant string := Greeting_Banners.witban;

-- THE HANDLERS :

  procedure General_Help ( opt : in character ) is
  begin
    case opt is
      when '0' => Greeting_Banners.help4setseed;
      when 'a' => Greeting_Banners.help4eqnbyeqn;
      when 'b' => Greeting_Banners.help4blackbox;
      when 'B' => Greeting_Banners.help4compsolve;
      when 'c' => Greeting_Banners.help4components;
      when 'd' => Greeting_Banners.help4reduction;
      when 'e' => Greeting_Banners.help4enumeration;
      when 'f' => Greeting_Banners.help4factor;
      when 'g' => Greeting_Banners.help4goodformat;
      when 'h' | '-' => Greeting_Banners.help4help;
      when 'j' => Greeting_Banners.help4adepath;
      when 'k' => Greeting_Banners.help4feedback;
      when 'l' => Greeting_Banners.help4hypersurface;
      when 'm' => Greeting_Banners.help4mixvol;
      when 'o' => Greeting_Banners.help4symbols;
      when 'p' => Greeting_Banners.help4continuation;
      when 'q' => Greeting_Banners.help4jumpstart;
      when 'r' => Greeting_Banners.help4rootcounts;
      when 's' => Greeting_Banners.help4scaling;
      when 't' => Greeting_Banners.help4tasking;
      when 'u' => Greeting_Banners.help4series;
      when 'v' => Greeting_Banners.help4verification;
      when 'V' => Greeting_Banners.help4verbose;
      when 'w' => Greeting_Banners.help4witsetinsect;
      when 'x' => Greeting_Banners.help4pythondict;
      when 'y' => Greeting_Banners.help4sampler;
      when 'z' => Greeting_Banners.help4mapleform;
      when others => Greeting_Banners.show_help;
    end case;
  end General_Help;

  procedure General_Help_Handler ( opts : in string ) is
  begin
    if opts'last > 1
     then General_Help(opts(opts'first+1));
     else General_Help(' ');
    end if;
  end General_Help_Handler;

  procedure Help_Version_License
              ( args : in Array_of_Strings; name : in string ) is

    arg : constant string := args(1).all;

  begin
    if arg = "--help" then
      Greeting_Banners.show_help;
    elsif arg = "--version" then
      if name = "" then
        put_line(Greeting_Banners.Version);
      else
        declare
          file : file_type;
        begin
          create(file,out_file,name);
          put_line(file,Greeting_Banners.Version);
          close(file);
        end;
      end if;
    elsif arg = "--license" then
      put_line("PHCpack is free and open source software.");
      put_line("You can redistribute the code and/or modify it under");
      put_line("the GNU General Pulic License as published by");
      put_line("the Free Software Foundation.");
    elsif arg = "--cite" then
      Greeting_Banners.How_to_Cite;
    else
      put_line(arg & " is not recognized.");
    end if;
  end Help_Version_License;

  procedure Full_Mode_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    optstr : constant string := Actions_and_Options.options;
    firstopt : constant character := opts(opts'first);
    pos : constant integer32 := Actions_and_Options.Position(optstr,firstopt);
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if pos < integer32(optstr'first) then
      put("The option '"); put(firstopt);
      put_line("' is not recognised.  Will ignore it...");
    end if;
    put_line(welcome);
    if nt > 0
     then mainphc(nt,infile,outfile,vrblvl);
     else mainphc(0,infile,outfile,vrblvl);
    end if;
  end Full_Mode_Handler;

  procedure Good_Format_Handler
              ( opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4goodformat;
    else
      if infile = "" or outfile = ""  -- still interactive input needed
       then put_line(welcome); put_line(goodban);
      end if;
      maingood(infile,outfile);
    end if;
  end Good_Format_Handler;

  procedure Find_and_Set_Seed
              ( args : in Array_of_Strings; opts : in string ) is

    pos : constant integer32 := Actions_and_Options.Position(opts,'0');
    seed_found : natural32;

  begin
    if pos >= integer32(opts'first) then
      seed_found := Actions_and_Options.Find_Seed(args);
      if seed_found = 0
       then Standard_Random_Numbers.Set_Seed(fixed_seed);
       else Standard_Random_Numbers.Set_Seed(seed_found);
      end if;
    end if;
  end Find_and_Set_Seed;

  procedure EqnByEqn_Solver_Handler
              ( opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4eqnbyeqn;
    else
      put_line(welcome); put_line(slvban);
      mainsolve(infile,outfile);
    end if;
  end EqnByEqn_Solver_Handler;

  procedure BlackBox_Solver_Handler
              ( args : in Array_of_Strings; opts : in string;
                file1,file2,file3 : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    epos : constant integer32 := Actions_and_Options.Position(opts,'e');
    spos : constant integer32 := Actions_and_Options.Position(opts,'s');
    dpos : constant integer32 := Actions_and_Options.Position(opts,'d');
    rpos : constant integer32 := Actions_and_Options.Position(opts,'r');
    mpos : constant integer32 := Actions_and_Options.Position(opts,'m');
    ppos : constant integer32 := Actions_and_Options.Position(opts,'p');
    vpos : constant integer32 := Actions_and_Options.Position(opts,'v');
    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    redprc : constant natural32
           := Actions_and_Options.Scan_Precision(args,'d');
    bbprc : constant natural32
          := Actions_and_Options.Scan_Precision(args,'b');
    valiprc : constant natural32
            := Actions_and_Options.Scan_Precision(args,'v');
    contprc : constant natural32
            := Actions_and_Options.Scan_Precision(args,'p');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4blackbox;
    elsif spos >= integer32(opts'first) then
      mainscal(file1,file2,vrblvl);
    elsif epos >= integer32(opts'first) then
      bablenum(file1,file2,vrblvl);
    elsif dpos >= integer32(opts'first) then
      if redprc = 2 then
        mainred2(file1,file2,vrblvl);
      elsif redprc = 4 then
        mainred4(file1,file2,vrblvl);
      else
        mainred(file1,file2,vrblvl);
      end if;
    elsif rpos >= integer32(opts'first) then
      if nt > 0
       then bablroco(nt,file1,file2,vrblvl);
       else bablroco(0,file1,file2,vrblvl);
      end if;
    elsif mpos >= integer32(opts'first) then
      if nt > 0
       then babldmvc(nt,file1,file2);
       else babldmvc(0,file1,file2);
      end if;
    elsif ppos >= integer32(opts'first) then
      if bbprc = 2 or contprc = 2 then
        bablpoco2(file1,file2,file3,vrblvl);
      elsif bbprc = 4 or contprc = 4 then
        bablpoco4(file1,file2,file3,vrblvl);
      else
        bablpoco(file1,file2,file3,vrblvl);
      end if;
    elsif vpos >= integer32(opts'first) then
      if bbprc = 2 or valiprc = 2 then
        bablvali2(file1,file2,vrblvl);
      elsif bbprc = 4 or valiprc = 4 then
        bablvali4(file1,file2,vrblvl);
      else
        bablvali(file1,file2,vrblvl);
      end if;
    elsif nt > 0 then
      case bbprc is
        when 2 => bablphc2(nt,file1,file2,vrblvl);
        when 4 => bablphc4(nt,file1,file2,vrblvl);
        when others => bablphc(nt,file1,file2,vrblvl);
      end case;
    else
      case bbprc is
        when 2 => bablphc2(0,file1,file2,vrblvl);
        when 4 => bablphc4(0,file1,file2,vrblvl);
        when others => bablphc(0,file1,file2,vrblvl);
      end case;
    end if;
  end BlackBox_Solver_Handler;

  procedure Component_Solver_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    compprc : constant natural32
            := Actions_and_Options.Scan_Precision(args,'B');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);
    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4compsolve;
    else
      if nt > 0 then
        case compprc is
          when 2 => compsolve2(nt,infile,outfile,vrblvl);
          when 4 => compsolve4(nt,infile,outfile,vrblvl);
          when others => compsolve(nt,infile,outfile,vrblvl);
        end case;
      else 
        case compprc is
          when 2 => compsolve2(0,infile,outfile,vrblvl);
          when 4 => compsolve4(0,infile,outfile,vrblvl);
          when others => compsolve(0,infile,outfile,vrblvl);
        end case;
      end if;
    end if;
  end Component_Solver_Handler;

  procedure Scaling_Handler 
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos : constant integer32 := Actions_and_Options.Position(opts,'h');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos >= integer32(opts'first) then
      Greeting_Banners.help4scaling;
    else
      put_line(welcome); put_line(scalban);
      mainscal(infile,outfile,vrblvl);
    end if;
  end Scaling_Handler;

  procedure Reduction_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    redprc : constant natural32
           := Actions_and_Options.Scan_Precision(args,'r');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4reduction;
    else
      put_line(welcome); put_line(reduban);
      if redprc = 2 then
        mainred2(infile,outfile,vrblvl);
      elsif redprc = 4 then
        mainred4(infile,outfile,vrblvl);
      else
        mainred(infile,outfile,vrblvl);
      end if;
    end if;
  end Reduction_Handler;

  procedure Root_Count_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    bpos : constant integer32 := Actions_and_Options.Position(opts,'b');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);
    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4rootcounts;
    elsif bpos >= integer32(opts'first) then
      if nt > 0 
       then bablroco(nt,infile,outfile,vrblvl);
       else bablroco(0,infile,outfile,vrblvl);
      end if;
    elsif nt = 0 then
      put_line(welcome); put_line(rocoban & ", no multitasking.");
      mainroco(0,infile,outfile,vrblvl);
    else    
      declare
        ns : constant string := Convert(integer32(nt));
      begin
        put_line(welcome); put_line(rocoban & ", with " & ns & " tasks.");
        mainroco(nt,infile,outfile,vrblvl);
      end;
    end if;
  end Root_Count_Handler;

  procedure Mixed_Volume_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    bpos : constant integer32 := Actions_and_Options.Position(opts,'b');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);
    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4mixvol;
    elsif bpos >= integer32(opts'first) then
      if nt > 0 
       then babldmvc(nt,infile,outfile);
       else babldmvc(0,infile,outfile);
      end if;
    elsif nt > 0 then
      declare
        ns : constant string := Convert(integer32(nt));
      begin
        put_line(welcome); put_line(mvcban & ", with " & ns & " tasks.");
        mainsmvc(nt,infile,outfile,vrblvl);
      end;
    else
      put_line(welcome); put_line(mvcban & ", no multitasking.");
      mainsmvc(0,infile,outfile,vrblvl);
    end if;
  end Mixed_Volume_Handler;

  procedure Symbols_Handler
              ( opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4symbols;
    else
      put_line(welcome); put_line(symbban);
      mainsymb(infile,outfile);
    end if;
  end Symbols_Handler;

  procedure Continuation_Handler
              ( args : in Array_of_Strings; opts : in string;
                file1,file2,file3 : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    bpos : constant integer32 := Actions_and_Options.Position(opts,'b');
    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    contprc : constant natural32
            := Actions_and_Options.Scan_Precision(args,'p');
    bbprc : constant natural32 := Actions_and_Options.Scan_Precision(args,'b');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4continuation;
    elsif bpos >= integer32(opts'first) then
      if contprc = 2 or bbprc = 2 then
        bablpoco2(file1,file2,file3,vrblvl);
      elsif contprc = 4 or bbprc = 4 then
        bablpoco4(file1,file2,file3,vrblvl);
      else
        bablpoco(file1,file2,file3,vrblvl);
      end if;
    elsif nt = 0 then
      put_line(welcome); put_line(pocoban & ", no multitasking.");
      mainpoco(0,file1,file2,contprc,vrblvl);
    else
      declare
        ns : constant string := Convert(integer32(nt));
      begin
        put_line(welcome);
        put_line(pocoban & ", with " & ns & " tasks.");
        mainpoco(nt,file1,file2,contprc,vrblvl);
      end;
    end if;
  end Continuation_Handler;

  procedure Jumpstart_Handler
              ( args : in Array_of_Strings;
                opts : in string; file1,file2,file3 : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4jumpstart;
    else
      put_line(welcome); put_line(trackban);
      maintrack(file1,file2,file3,vrblvl);
    end if;
  end Jumpstart_Handler;

  procedure Algorithmic_Differentiation_Handler
              ( opts : in string; file1,file2,file3 : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4adepath;
    else
      put_line(welcome); put_line(adepban);
      mainadep(file1,file2,file3);
    end if;
  end Algorithmic_Differentiation_Handler;

  procedure Enumeration_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    bpos : constant integer32 := Actions_and_Options.Position(opts,'b');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4enumeration;
    elsif bpos >= integer32(opts'first) then
      bablenum(infile,outfile,vrblvl);
    else
      put_line(welcome); put_line(enumban);
      mainenum(vrblvl);
    end if;
  end Enumeration_Handler;

  procedure Feedback_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first)
     then Greeting_Banners.help4feedback;
     else mainfeed(infile,outfile,vrblvl);
    end if;
  end Feedback_Handler;

  procedure Decomposition_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4components;
    else
      if nt > 0 then
        declare
          ns : constant string := Convert(integer32(nt));
        begin
          put_line(welcome);
          put_line(compban & " with " & ns & " tasks");
          maindeco(nt,infile,outfile,vrblvl);
        end;
      else
        put_line(welcome); put_line(compban);
        maindeco(0,infile,outfile,vrblvl);
      end if;
    end if;
  end Decomposition_Handler;

  procedure Factorization_Handler
              ( args : in Array_of_Strings; opts : in string;
                infile,outfile : in string ) is

    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4factor;
    else
      put_line(welcome);
      if nt = 0 then
        put_line(facban & ", no multitasking.");
      else
        declare
          ns : constant string := Convert(integer32(nt));
        begin
          put_line(facban & ", with " & ns & " tasks.");
        end;
      end if;
      mainfac(nt,infile,outfile);
    end if;
  end Factorization_Handler;

  procedure Series_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    prc : constant natural32
        := Actions_and_Options.Scan_Precision(args,'u');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4series;
    else
      put_line(welcome);
      case prc is
        when 1 =>
          put_line(seriesban & ", in double precision.");
          mainseries('1',infile,outfile,vrblvl);
        when 2 =>
          put_line(seriesban & ", with double doubles.");
          mainseries('2',infile,outfile,vrblvl);
        when 4 =>
          put_line(seriesban & ", with quad doubles.");
          mainseries('4',infile,outfile,vrblvl);
        when others =>
          put_line(seriesban & ".");
          mainseries('0',infile,outfile,vrblvl);
      end case;
    end if;
  end Series_Handler;

  procedure Verification_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    bpos : constant integer32 := Actions_and_Options.Position(opts,'b');
    veriprc : constant natural32
            := Actions_and_Options.Scan_Precision(args,'v');
    bbprc : constant natural32
          := Actions_and_Options.Scan_Precision(args,'b');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4verification;
    elsif bpos >= integer32(opts'first) then 
      if veriprc = 2 or bbprc = 2 then
        bablvali2(infile,outfile,vrblvl);
      elsif veriprc = 4 or bbprc = 4 then
        bablvali4(infile,outfile,vrblvl);
      else
        bablvali(infile,outfile,vrblvl);
      end if;
    else
      put_line(welcome); put_line(veriban);
      mainvali(infile,outfile,vrblvl);
    end if;
  end Verification_Handler;

  procedure Witness_Set_for_Hypersurface_Handler
              ( args : in Array_of_Strings; opts : in string;
                polyfile,logfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    precision : constant natural32
              := Actions_and_Options.Scan_Precision(args,'l');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4hypersurface;
    else
      if polyfile = "" or logfile = ""
       then put_line(welcome); put_line(hypban);
      end if;
      case precision is
        when 2 => mainhyp2(polyfile,logfile);
        when 4 => mainhyp4(polyfile,logfile);
        when others => mainhyp(polyfile,logfile);
      end case;
    end if;
  end Witness_Set_for_Hypersurface_Handler;

  procedure Witness_Set_Intersection_Handler
              ( opts,witset_one,witset_two,logfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4witsetinsect;
    else
      if witset_one = "" or witset_two = "" or logfile = ""
       then put_line(welcome); put_line(witban);
      end if;
      mainwit(witset_one,witset_two,logfile);
    end if;
  end Witness_Set_Intersection_Handler;

  procedure Witness_Set_Sampler_Handler
              ( opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first) then
      Greeting_Banners.help4sampler;
    else
      put_line(welcome); put_line(samban);
      mainsam(infile,outfile);
    end if;
  end Witness_Set_Sampler_Handler;

  procedure Maple_Format_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first)
     then Greeting_Banners.help4mapleform;
     else mainzip(infile,outfile,vrblvl);
    end if;
  end Maple_Format_Handler;

  procedure Python_Format_Handler
              ( args : in Array_of_Strings;
                opts : in string; infile,outfile : in string ) is

    hpos1 : constant integer32 := Actions_and_Options.Position(opts,'h');
    hpos2 : constant integer32 := Actions_and_Options.Position(opts,'-');
    vrblvl : constant integer32 := Actions_and_Options.Verbose_Level(args);

  begin
    if hpos1 >= integer32(opts'first) or hpos2 >= integer32(opts'first)
     then Greeting_Banners.help4pythondict;
     else maindict(infile,outfile,vrblvl);
    end if;
  end Python_Format_Handler;

-- THE MAIN HANDLERS :

  procedure Handle_Options ( args : in Array_of_Strings; 
                             opts,a1,a2,a3 : in string;
                             verbose : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   The handler after checking the case where the first two
  --   file names in a1 and a2 would be the same nonempty string.

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in option_handlers.Handle_Options ...");
    end if;
    Find_and_Set_Seed(args,opts);
    if opts'length = 0 then
      Handle_no_Options(a1,a2);
    else
      case opts(opts'first) is
        when '-' => Help_Version_License(args,a1);
        when 'g' => Good_Format_Handler(opts,a1,a2);
        when 'h' => General_Help_Handler(opts);
        when 'a' => EqnByEqn_Solver_Handler(opts,a1,a2);
        when 'b' => BlackBox_Solver_Handler(args,opts,a1,a2,a3);
        when 'B' => Component_Solver_Handler(args,opts,a1,a2);
        when 's' => Scaling_Handler(args,opts,a1,a2);
        when 'd' => Reduction_Handler(args,opts,a1,a2);
        when 'r' => Root_Count_Handler(args,opts,a1,a2);
        when 'm' => Mixed_Volume_Handler(args,opts,a1,a2);
        when 'o' => Symbols_Handler(opts,a1,a2);
        when 'p' => Continuation_Handler(args,opts,a1,a2,a3);
        when 'q' => Jumpstart_Handler(args,opts,a1,a2,a3);
        when 'j' => Algorithmic_Differentiation_Handler(opts,a1,a2,a3);
        when 'c' => Decomposition_Handler(args,opts,a1,a2);
        when 'e' => Enumeration_Handler(args,opts,a1,a2);
        when 'k' => Feedback_Handler(args,opts,a1,a2);
        when 'f' => Factorization_Handler(args,opts,a1,a2);
        when 'u' => Series_Handler(args,opts,a1,a2);
        when 'v' => Verification_Handler(args,opts,a1,a2);
        when 'l' => Witness_Set_for_Hypersurface_Handler(args,opts,a1,a2);
        when 'w' => Witness_Set_Intersection_Handler(opts,a1,a2,a3);
        when 'x' => Python_Format_Handler(args,opts,a1,a2);
        when 'y' => Witness_Set_Sampler_Handler(opts,a1,a2);
        when 'z' => Maple_Format_Handler(args,opts,a1,a2);
        when others => Full_Mode_Handler(args,opts,a1,a2);
      end case;
    end if;
  exception
    when others =>
      put_line("Oops, an unhandled exception occurred.");
      Write_Seed_Number(standard_output);
      put_line("Use the seed number to reproduce the error."); raise;
  end Handle_Options;

  procedure Handle ( args : in Array_of_Strings; 
                     opts,a1,a2,a3 : in string;
                     verbose : in integer32 := 0 ) is
  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in option_handlers.Handle ...");
    end if;
    if a1 /= "" and then (a1 = a2) then
      new_line;
      put_line("The first two file names are identical.");
      put_line("Will ignore the second file name.");
      Handle_Options(args,opts,a1,"",a3,verbose-1);
    else
      Handle_Options(args,opts,a1,a2,a3,verbose-1);
    end if;
  end Handle;

  procedure Handle_no_Options ( infile,outfile : in string ) is
  begin
    put_line(welcome);
    mainphc(0,infile,outfile);
  exception
    when others =>
      put_line("Oops, an unhandled exception occurred.");
      Write_Seed_Number(standard_output);
      put_line("Use the seed number to reproduce the error."); raise;
  end Handle_no_Options;

end Option_Handlers;
