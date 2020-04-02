with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Standard_Homotopy;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy;
with QuadDobl_Homotopy_Convolutions_io;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Residual_Convolution_Circuits;      use Residual_Convolution_Circuits;
with Predictor_Corrector_Loops;          use Predictor_Corrector_Loops;
with Standard_Solutions_Queue;
with DoblDobl_Solutions_Queue;
with QuadDobl_Solutions_Queue;
with Multitasking;

procedure ts_mtpcscnv is

-- DESCRIPTION :
--   Development of multitasked tracking with predictor-corrector-shift
--   loops on homotopy systems of convolution circuits.

  procedure Standard_Multitasked_Tracker
              ( file : in file_type; nbtasks : in integer32;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Applies multitasking to track all paths in double precision.

  -- ON ENTRY :
  --   file     to write output information to;
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

    use Standard_Complex_Solutions,Standard_Speelpenning_Convolutions;

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;

    begin
      loop
        myptr := Standard_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        cnt := Standard_Solutions_Queue.Next_Counter;
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks  path "
                         & Multitasking.to_string(cnt));
        delay 0.1;
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    Standard_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
  end Standard_Multitasked_Tracker;

  procedure DoblDobl_Multitasked_Tracker
              ( file : in file_type; nbtasks : in integer32;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Applies multitasking to track all paths in double double precision.

  -- ON ENTRY :
  --   file     to write output information to;
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

    use DoblDobl_Complex_Solutions,DoblDobl_Speelpenning_Convolutions;

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;

    begin
      loop
        myptr := DoblDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        cnt := DoblDobl_Solutions_Queue.Next_Counter;
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks  path "
                         & Multitasking.to_string(cnt));
        delay 0.1;
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    DoblDobl_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
  end DoblDobl_Multitasked_Tracker;

  procedure QuadDobl_Multitasked_Tracker
              ( file : in file_type; nbtasks : in integer32;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                verbose : in boolean := true ) is

  -- DESCRIPTION :
  --   Applies multitasking to track all paths in quad double precision.

  -- ON ENTRY :
  --   file     to write output information to;
  --   nbtasks  the number of tasks;
  --   hom      system of homotopy convolution circuits;
  --   abh      radii as coefficients for mixed residuals;
  --   sols     start solutions;
  --   pars     values for the tolerances and parameters;
  --   verbose  indicates if extra output is requested.
  
  -- ON RETURN :
  --   sols     solutions at the end of the paths.

    use QuadDobl_Complex_Solutions,QuadDobl_Speelpenning_Convolutions;

    procedure Silent_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks without intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
      end loop;
    end Silent_Track;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Track);

    procedure Report_Track ( i,n : in integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n tracks with intermediate output.

      myptr : Solution_List;
      ls : Link_to_Solution;
      cnt : integer32;

    begin
      loop
        myptr := QuadDobl_Solutions_Queue.Next;
        exit when Is_Null(myptr);
        ls := Head_Of(myptr);
        cnt := QuadDobl_Solutions_Queue.Next_Counter;
        put_line("Task " & Multitasking.to_string(i)
                         & " tracks  path "
                         & Multitasking.to_string(cnt));
        delay 0.1;
      end loop;
    end Report_Track;
    procedure report_do_jobs is
      new Multitasking.Reporting_Workers(Report_Track);

  begin
    QuadDobl_Solutions_Queue.Initialize(sols);
    if verbose
     then report_do_jobs(nbtasks);
     else silent_do_jobs(nbtasks);
    end if;
  end QuadDobl_Multitasked_Tracker;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double precision.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg,nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    ans : character;
    file : file_type;
    verbose : boolean;

  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    Standard_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    pars.gamma := Standard_Homotopy.Accessibility_Constant;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Give the number of tasks : "); get(nbt);
    Standard_Multitasked_Tracker(file,nbt,cnvhom,abshom,sols,pars,verbose);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double double precisin.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg,nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    ddgamma : DoblDobl_Complex_Numbers.Complex_Number;
    ans : character;
    file : file_type;
    verbose : boolean;

    use DoblDobl_Complex_Numbers_cv;
  
  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    DoblDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    ddgamma := DoblDobl_Homotopy.Accessibility_Constant;
    pars.gamma := DoblDobl_Complex_to_Standard(ddgamma);
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Give the number of tasks : "); get(nbt);
    DoblDobl_Multitasked_Tracker(file,nbt,cnvhom,abshom,sols,pars,verbose);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in quad double precision.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    idxpar,deg,nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    qdgamma : QuadDobl_Complex_Numbers.Complex_Number;
    ans : character;
    file : file_type;
    verbose : boolean;

    use QuadDobl_Complex_Numbers_cv;

  begin
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    QuadDobl_Homotopy_Convolutions_io.get(deg,cnvhom,sols,idxpar);
    abshom := Residual_Convolution_System(cnvhom);
    qdgamma := QuadDobl_Homotopy.Accessibility_Constant;
    pars.gamma := QuadDobl_Complex_to_Standard(qdgamma);
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put("Give the number of tasks : "); get(nbt);
    QuadDobl_Multitasked_Tracker(file,nbt,cnvhom,abshom,sols,pars,verbose);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the working precision and then launches
  --   the test corresponding to the selected precision.

    precision : constant character := Prompt_for_Precision;

  begin
    case precision is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtpcscnv;
