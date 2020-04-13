with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Homotopy_Continuation_Parameters;
with Homotopy_Continuation_Parameters_io;
with Series_Path_Trackers;               use Series_Path_Trackers;
with Standard_Speelpenning_Convolutions;
with Standard_Homotopy_Convolutions_io;

procedure ts_scalecnv is

-- DESCRIPTION :
--   Development of homotopy convolution circuits with scaling of
--   the solutions, after a projective coordinate transformation.

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the settings of the homotopy.

    art : constant boolean := Prompt_for_Artificial;
    deg,idxpar : integer32;
    pars : Homotopy_Continuation_Parameters.Parameters
         := Homotopy_Continuation_Parameters.Default_Values;
    cnvhom : Standard_Speelpenning_Convolutions.Link_to_System;
    mhom : natural32 := 0;
    z : Link_to_Partition;
    idz : Standard_Natural_Vectors.Link_to_Vector;
    sols : Standard_Complex_Solutions.Solution_List;

  begin
    new_line;
    if not art
     then pars.gamma := Standard_Complex_Numbers.Create(1.0);
    end if;
    Homotopy_Continuation_Parameters_io.Tune(pars);
    deg := integer32(pars.numdeg + pars.dendeg + 2);
    Standard_Homotopy_Convolutions_io.get
      (deg,art,pars.gamma,cnvhom,sols,idxpar,mhom,z,idz);
  end Main;

begin
  Main;
end ts_scalecnv;
