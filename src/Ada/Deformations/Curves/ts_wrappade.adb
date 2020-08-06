with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Homotopy;
with Standard_Complex_Solutions;
with Standard_System_and_Solutions_io;
with Standard_Parameter_Systems;
with Wrapped_Pade_Trackers;

procedure ts_wrappade is

-- DESCRIPTION :
--   Development of the wrapping of the Pade trackers.

  procedure Standard_Track
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Wrapped_Pade_Trackers.Run(file,p'last,p,ls.v,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Track;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Reads a natural parameter in double precision.
  --   Tunes the parameters and then calls the path tracker.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    outfile : file_type;
    ans : character;
    nbequ,nbvar,idxpar,nbpar : integer32;

    use Standard_Parameter_Systems;

  begin
    Read_Parameter_Homotopy(lp,sols,nbequ,nbvar,nbpar);
    declare
      par : Standard_Integer_Vectors.Vector(1..nbpar);
    begin
      par := Define_Parameters(nbequ,nbvar,nbpar);
      idxpar := par(1);
    end;
    Standard_Homotopy.Create(lp.all,idxpar);
    new_line;
    put("Output to file ? "); Ask_Yes_or_No(ans);
    if ans /= 'y' then
      new_line;
      Wrapped_Pade_Trackers.Set_Parameters;
      Standard_Track(standard_output,lp.all,sols);
    else
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
      Standard_System_and_Solutions_io.put(outfile,lp.all,sols);
      new_line(outfile);
      new_line;
      Wrapped_Pade_Trackers.Set_Parameters(outfile);
      Standard_Track(outfile,lp.all,sols);
      close(outfile);
    end if;
  end Standard_Main;

begin
  Standard_Main;
end ts_wrappade;
