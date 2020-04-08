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
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
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
    Standard_Multitasked_Tracker(nbt,cnvhom,abshom,sols,pars,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,Standard_Complex_Solutions.Length_Of(sols),
        natural32(Standard_Complex_Solutions.Head_Of(sols).n),sols);
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
    DoblDobl_Multitasked_Tracker(nbt,cnvhom,abshom,sols,pars,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,DoblDobl_Complex_Solutions.Length_Of(sols),
        natural32(DoblDobl_Complex_Solutions.Head_Of(sols).n),sols);
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
    QuadDobl_Multitasked_Tracker(nbt,cnvhom,abshom,sols,pars,verbose);
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(file,QuadDobl_Complex_Solutions.Length_Of(sols),
        natural32(QuadDobl_Complex_Solutions.Head_Of(sols).n),sols);
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
