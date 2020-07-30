with Ada.Calendar;
with Communications_with_User;            use Communications_with_User;
with Write_Seed_Number;
with Time_Stamps;                         use Time_Stamps;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_cv;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_cv;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;    use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;    use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;       use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;       use QuadDobl_Complex_Solutions_io;
with Projective_Transformations;
with Multi_Projective_Transformations;
with Partitions_of_Sets_of_Unknowns;      use Partitions_of_Sets_of_Unknowns;
with Standard_Homotopy;
with Standard_Homotopy_Convolutions_io;
with DoblDobl_Homotopy;
with DoblDobl_Homotopy_Convolutions_io;
with QuadDobl_Homotopy;
with QuadDobl_Homotopy_Convolutions_io;
with Homotopy_Continuation_Parameters_io;
with Standard_Circuit_Makers;
with Standard_Convolution_Splitters;
with Predictor_Corrector_Trackers;        use Predictor_Corrector_Trackers;
with Residual_Convolution_Circuits;       use Residual_Convolution_Circuits;
with Series_Path_Trackers;
with Greeting_Banners;

package body Track_Path_Convolutions is

  procedure Standard_Write_Homotopy
              ( file : in file_type; neq : in integer32;
                sols : in Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                arth : in boolean; verbose : out boolean ) is

    ans : character;

  begin
    if not arth then
      put(file,natural32(neq),natural32(neq+1),
               Standard_Homotopy.Homotopy_System);
    else
      declare
        p : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Homotopy.Target_System;
        q : constant Standard_Complex_Poly_Systems.Poly_Sys
          := Standard_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
             natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
  end Standard_Write_Homotopy;               

  procedure DoblDobl_Write_Homotopy
              ( file : in file_type; neq : in integer32;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                arth : in boolean; verbose : out boolean ) is

    ans : character;

  begin
    if not arth then
      put(file,natural32(neq),natural32(neq+1),
               DoblDobl_Homotopy.Homotopy_System);
    else
      declare
        p : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := DoblDobl_Homotopy.Target_System;
        q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys
          := DoblDobl_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
             natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
  end DoblDobl_Write_Homotopy;

  procedure QuadDobl_Write_Homotopy
              ( file : in file_type; neq : in integer32;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                arth : in boolean; verbose : out boolean ) is

    ans : character;

  begin
    if not arth then
      put(file,natural32(neq),natural32(neq+1),
               QuadDobl_Homotopy.Homotopy_System);
    else
      declare
        p : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := QuadDobl_Homotopy.Target_System;
        q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys
          := QuadDobl_Homotopy.Start_System;
      begin
        put(file,p'last,1); new_line(file); put(file,p);
        new_line(file);
        put_line(file,"THE START SYSTEM :");
        put(file,q'last,1); new_line(file); put(file,q);
      end;
    end if;
    new_line(file);
    put_line(file,"THE START SOLUTIONS :");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
             natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    new_line(file);
    Homotopy_Continuation_Parameters_io.put(file,pars); flush(file);
    new_line;
    put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
    verbose := (ans = 'y');
    new_line;
    put_line("See the output file for results ...");
    new_line;
  end QuadDobl_Write_Homotopy;

  procedure Standard_Write_Solutions 
              ( file : in file_type; arth : in boolean;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

    hcrd : constant boolean := (mhom > 0);

  begin
    new_line(file);
    if arth and hcrd then
      if mhom = 1 then
        put_line(file,"THE 1-HOMOGENEOUS SOLUTIONS :");
      else
        put(file,"THE "); put(file,mhom,1);
        put_line(file,"-HOMOGENEOUS SOLUTIONS :");
      end if;
    else
      put_line(file,"THE SOLUTIONS :");
    end if;
    put(file,Standard_Complex_Solutions.Length_Of(sols),
             natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    if arth and hcrd then
      if mhom = 1
       then Projective_Transformations.Affine_Transformation(sols);
       else Multi_Projective_Transformations.Make_Affine
              (sols,natural32(mhom),idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Standard_Complex_Solutions.Length_Of(sols),
               natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
  end Standard_Write_Solutions;

  procedure DoblDobl_Write_Solutions 
              ( file : in file_type; arth : in boolean;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

    hcrd : constant boolean := (mhom > 0);

  begin
    new_line(file);
    if arth and hcrd then
      if mhom = 1 then
        put_line(file,"THE 1-HOMOGENEOUS SOLUTIONS :");
      else
        put(file,"THE "); put(file,mhom,1);
        put_line(file,"-HOMOGENEOUS SOLUTIONS :");
      end if;
    else
      put_line(file,"THE SOLUTIONS :");
    end if;
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
             natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    if arth and hcrd then
      if mhom = 1
       then Projective_Transformations.Affine_Transformation(sols);
       else Multi_Projective_Transformations.Make_Affine
              (sols,natural32(mhom),idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
               natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
  end DoblDobl_Write_Solutions;

  procedure QuadDobl_Write_Solutions 
              ( file : in file_type; arth : in boolean;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

    hcrd : constant boolean := (mhom > 0);

  begin
    new_line(file);
    if arth and hcrd then
      if mhom = 1 then
        put_line(file,"THE 1-HOMOGENEOUS SOLUTIONS :");
      else
        put(file,"THE "); put(file,mhom,1);
        put_line(file,"-HOMOGENEOUS SOLUTIONS :");
      end if;
    else
      put_line(file,"THE SOLUTIONS :");
    end if;
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
             natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    if arth and hcrd then
      if mhom = 1
       then Projective_Transformations.Affine_Transformation(sols);
       else Multi_Projective_Transformations.Make_Affine
              (sols,natural32(mhom),idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
               natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
  end QuadDobl_Write_Solutions;

-- ON COEFFICIENT CONVOLUTION CIRCUITS :

  procedure Track
              ( file : in file_type;
                hom : in Standard_Coefficient_Convolutions.Link_to_System;
                cfh,abh : in Standard_Coefficient_Circuits.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean; vrb : in integer32 := 0 ) is

    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Track 1 ...");
    end if;
    Standard_Write_Homotopy(file,hom.neq,sols,pars,arth,verbose);
    new_line(file);
    Track_All_Paths(file,hom,cfh,abh,sols,pars,mhom,idz,verbose,vrb-1);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
             natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
  end Track;

-- ON COMPLEX COEFFICIENT CIRCUITS :

  procedure Track
              ( file : in file_type;
                hom : in Standard_Speelpenning_Convolutions.Link_to_System;
                abh : in Standard_Speelpenning_Convolutions.Link_to_System;
                sols : in out Standard_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean; vrb : in integer32 := 0 ) is

    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Track 2 ...");
    end if;
    Standard_Write_Homotopy(file,hom.neq,sols,pars,arth,verbose);
    new_line(file);
    Track_All_Paths(file,hom,abh,sols,pars,mhom,idz,verbose,vrb-1);
    Standard_Write_Solutions(file,arth,mhom,idz,sols);
  end Track;

  procedure Track
              ( file : in file_type;
                hom : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in DoblDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean; vrb : in integer32 := 0 ) is

    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Track 3 ...");
    end if;
    DoblDobl_Write_Homotopy(file,hom.neq,sols,pars,arth,verbose);
    new_line(file);
    Track_All_Paths(file,hom,abh,sols,pars,mhom,idz,verbose,vrb-1);
    DoblDobl_Write_Solutions(file,arth,mhom,idz,sols);
  end Track;

  procedure Track
              ( file : in file_type;
                hom : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : in QuadDobl_Speelpenning_Convolutions.Link_to_System;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                pars : in Homotopy_Continuation_Parameters.Parameters;
                mhom : in integer32;
                idz : in Standard_Natural_Vectors.Link_to_Vector;
                arth : in boolean; vrb : in integer32 := 0 ) is

    verbose : boolean;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Track 4 ...");
    end if;
    QuadDobl_Write_Homotopy(file,hom.neq,sols,pars,arth,verbose);
    new_line(file);
    Track_All_Paths(file,hom,abh,sols,pars,mhom,idz,verbose,vrb-1);
    QuadDobl_Write_Solutions(file,arth,mhom,idz,sols);
  end Track;

  procedure Main
              ( hom : out Standard_Speelpenning_Convolutions.Link_to_System;
                abh : out Standard_Speelpenning_Convolutions.Link_to_System;
                arth : out boolean;
                pars : out Homotopy_Continuation_Parameters.Parameters;
                sols : out Standard_Complex_Solutions.Solution_List;
                mhom : out natural32;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                vrb : in integer32 := 0 ) is

    idxpar,deg : integer32;
    z : Link_to_Partition;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Main 1 ...");
    end if;
    arth := Series_Path_Trackers.Prompt_for_Artificial;
    pars := Homotopy_Continuation_Parameters.Default_Values;
    if not arth
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    Standard_Homotopy_Convolutions_io.get
      (deg,arth,pars.gamma,hom,sols,idxpar,mhom,z,idz);
    abh := Residual_Convolution_System(hom);
    if arth
     then pars.gamma := Standard_Homotopy.Accessibility_Constant;
    end if;
  end Main;

  procedure Main
              ( hom : out DoblDobl_Speelpenning_Convolutions.Link_to_System;
                abh : out DoblDobl_Speelpenning_Convolutions.Link_to_System;
                arth : out boolean;
                pars : out Homotopy_Continuation_Parameters.Parameters;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                mhom : out natural32;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                vrb : in integer32 := 0 ) is

    idxpar,deg : integer32;
    z : Link_to_Partition;
    ddgamma : DoblDobl_Complex_Numbers.Complex_Number;

    use DoblDobl_Complex_Numbers_cv;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Main 2 ...");
    end if;
    arth := Series_Path_Trackers.Prompt_for_Artificial;
    pars := Homotopy_Continuation_Parameters.Default_Values;
    if not arth
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    DoblDobl_Homotopy_Convolutions_io.get
      (deg,arth,pars.gamma,hom,sols,idxpar,mhom,z,idz);
    abh := Residual_Convolution_System(hom);
    if arth then
      ddgamma := DoblDobl_Homotopy.Accessibility_Constant;
      pars.gamma := DoblDobl_Complex_to_Standard(ddgamma);
    end if;
  end Main;

  procedure Main
              ( hom : out QuadDobl_Speelpenning_Convolutions.Link_to_System;
                abh : out QuadDobl_Speelpenning_Convolutions.Link_to_System;
                arth : out boolean;
                pars : out Homotopy_Continuation_Parameters.Parameters;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                mhom : out natural32;
                idz : out Standard_Natural_Vectors.Link_to_Vector;
                vrb : in integer32 := 0 ) is

    idxpar,deg : integer32;
    z : Link_to_Partition;
    qdgamma : QuadDobl_Complex_Numbers.Complex_Number;

    use QuadDobl_Complex_Numbers_cv;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Main 3 ...");
    end if;
    arth := Series_Path_Trackers.Prompt_for_Artificial;
    pars := Homotopy_Continuation_Parameters.Default_Values;
    if not arth
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    new_line;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    QuadDobl_Homotopy_Convolutions_io.get
      (deg,arth,pars.gamma,hom,sols,idxpar,mhom,z,idz);
    abh := Residual_Convolution_System(hom);
    if arth then
      qdgamma := QuadDobl_Homotopy.Accessibility_Constant;
      pars.gamma := QuadDobl_Complex_to_Standard(qdgamma);
    end if;
  end Main;

  procedure Standard_Main ( vrb : in integer32 := 0 ) is

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    file : file_type;
    start_moment,ended_moment : Ada.Calendar.Time;
    ans : character;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.Standard_Main ...");
    end if;
    Main(cnvhom,abshom,artificial,pars,sols,mhom,idz,vrb-1);
    if mhom > 0 then
      ans := 'n'; -- homogenization not yet supported on coefficient conv
    else
     -- new_line;
     -- put("Running with coefficient convolution circuits ? (y/n) ");
     -- Ask_Yes_or_No(ans);
      ans := 'y'; -- make coefficient convolution circuits the default
    end if;
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    start_moment := Ada.Calendar.Clock;
    if ans = 'n' then
      Track(file,cnvhom,abshom,sols,pars,integer32(mhom),idz,artificial,vrb-1);
    else
      declare
        cffhom : Standard_Coefficient_Convolutions.Link_to_System;
        cfs,abh : Standard_Coefficient_Circuits.Link_to_System;
      begin
        cffhom := Standard_Convolution_Splitters.Split(cnvhom);
        cfs := Standard_Circuit_Makers.Make_Coefficient_System(cffhom);
        abh := Standard_Coefficient_Circuits.Copy(cfs);
        Standard_Coefficient_Circuits.AbsVal(abh);
        Track(file,cffhom,cfs,abh,sols,pars,
              integer32(mhom),idz,artificial,vrb-1);
      end;
    end if;
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    put(file,"PHC ran from "); Write_Time_Stamp(file,start_moment);
    put(file," till "); Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    Write_Elapsed_Time(file,start_moment,ended_moment);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end Standard_Main;

  procedure DoblDobl_Main ( vrb : in integer32 := 0 ) is

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    file : file_type;
    start_moment,ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.DoblDobl_Main ...");
    end if;
    Main(cnvhom,abshom,artificial,pars,sols,mhom,idz,vrb-1);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    start_moment := Ada.Calendar.Clock;
    Track(file,cnvhom,abshom,sols,pars,integer32(mhom),idz,artificial,vrb-1);
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    put(file,"PHC ran from "); Write_Time_Stamp(file,start_moment);
    put(file," till "); Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    Write_Elapsed_Time(file,start_moment,ended_moment);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end DoblDobl_Main;

  procedure QuadDobl_Main ( vrb : in integer32 := 0 ) is

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    pars : Homotopy_Continuation_Parameters.Parameters;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    artificial : boolean;
    file : file_type;
    start_moment,ended_moment : Ada.Calendar.Time;

  begin
    if vrb > 0
     then put_line("-> in track_path_convolutions.QuadDobl_Main ...");
    end if;
    Main(cnvhom,abshom,artificial,pars,sols,mhom,idz,vrb-1);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    start_moment := Ada.Calendar.Clock;
    Track(file,cnvhom,abshom,sols,pars,integer32(mhom),idz,artificial,vrb-1);
    ended_moment := Ada.Calendar.Clock;
    new_line(file);
    put(file,"PHC ran from "); Write_Time_Stamp(file,start_moment);
    put(file," till "); Write_Time_Stamp(file,ended_moment);
    put_line(file,".");
    Write_Elapsed_Time(file,start_moment,ended_moment);
    Write_Seed_Number(file);
    put_line(file,Greeting_Banners.Version);
  end QuadDobl_Main;

end Track_Path_Convolutions;
