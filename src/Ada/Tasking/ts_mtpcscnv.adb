with text_io;                             use text_io;
with Ada.Calendar;
with Communications_with_User;            use Communications_with_User;
with Time_Stamps;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;       use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;       use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;       use QuadDobl_Complex_Solutions_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;    use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;    use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;    use QuadDobl_Complex_Poly_Systems_io;
with Projective_Transformations;
with Multi_Projective_Transformations;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Track_Path_Convolutions;
with Multitasked_Path_Convolutions;      use Multitasked_Path_Convolutions;

procedure ts_mtpcscnv is

-- DESCRIPTION :
--   Tests multitasked tracking with predictor-corrector-shift
--   loops on homotopy systems of convolution circuits.

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double precision.

    sols : Standard_Complex_Solutions.Solution_List;
    cnvhom,abshom : Standard_Speelpenning_Convolutions.Link_to_System;
    nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters;
    ans : character;
    file : file_type;
    artificial,hcrd,verbose : boolean;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    startmoment,stopmoment : Ada.Calendar.Time;
  
  begin
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    hcrd := (mhom > 0);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    if not artificial then
      put(file,natural32(cnvhom.neq),natural32(cnvhom.neq+1),
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
    put("Give the number of tasks : "); get(nbt);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    startmoment := Ada.Calendar.Clock;
    Standard_Multitasked_Tracker
      (nbt,cnvhom,abshom,sols,pars,integer32(mhom),idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    new_line(file);
    if artificial and hcrd then
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
    if artificial and hcrd then
      if mhom = 1
       then Projective_Transformations.Affine_Transformation(sols);
       else Multi_Projective_Transformations.Make_Affine(sols,mhom,idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,Standard_Complex_Solutions.Length_Of(sols),
               natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in double double precisin.

    sols : DoblDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : DoblDobl_Speelpenning_Convolutions.Link_to_System;
    nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters;
    ans : character;
    file : file_type;
    artificial,hcrd,verbose : boolean;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    startmoment,stopmoment : Ada.Calendar.Time;

  begin
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    hcrd := (mhom > 0);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    if not artificial then
      put(file,natural32(cnvhom.neq),natural32(cnvhom.neq+1),
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
    put("Give the number of tasks : "); get(nbt);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    startmoment := Ada.Calendar.Clock;
    DoblDobl_Multitasked_Tracker
      (nbt,cnvhom,abshom,sols,pars,integer32(mhom),idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    new_line(file);
    if artificial and hcrd then
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
    if artificial and hcrd then
      if mhom = 1
       then Projective_Transformations.Affine_Transformation(sols);
       else Multi_Projective_Transformations.Make_Affine(sols,mhom,idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
               natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy in quad double precision.

    sols : QuadDobl_Complex_Solutions.Solution_List;
    cnvhom,abshom : QuadDobl_Speelpenning_Convolutions.Link_to_System;
    nbt : integer32 := 0;
    pars : Homotopy_Continuation_Parameters.Parameters;
    ans : character;
    file : file_type;
    artificial,hcrd,verbose : boolean;
    mhom : natural32 := 0;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    startmoment,stopmoment : Ada.Calendar.Time;

  begin
    Track_Path_Convolutions.Main(cnvhom,abshom,artificial,pars,sols,mhom,idz);
    hcrd := (mhom > 0);
    new_line;
    put_line("Reading the name of the output file ...");
    Read_Name_and_Create_File(file);
    if not artificial then
      put(file,natural32(cnvhom.neq),natural32(cnvhom.neq+1),
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
    put("Give the number of tasks : "); get(nbt);
    new_line;
    put_line("See the output file for results ...");
    new_line;
    startmoment := Ada.Calendar.Clock;
    QuadDobl_Multitasked_Tracker
      (nbt,cnvhom,abshom,sols,pars,integer32(mhom),idz,verbose);
    stopmoment := Ada.Calendar.Clock;
    new_line(file);
    if artificial and hcrd then
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
    if artificial and hcrd then
      if mhom = 1
       then Projective_Transformations.Affine_Transformation(sols);
       else Multi_Projective_Transformations.Make_Affine(sols,mhom,idz.all);
      end if;
      new_line(file);
      put_line(file,"THE SOLUTIONS :");
      put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
               natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
    end if;
    new_line(file);
    put(file,"Elapsed wall clock time with ");
    put(file,nbt,1); put_line(file," tasks :");
    Time_Stamps.Write_Elapsed_Time(file,startmoment,stopmoment);
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
