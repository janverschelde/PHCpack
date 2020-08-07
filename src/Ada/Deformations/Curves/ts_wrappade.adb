with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Homotopy;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with Standard_Parameter_Systems;
with Wrapped_Solution_Vectors;
with Wrapped_Pade_Trackers;

procedure ts_wrappade is

-- DESCRIPTION :
--   Development of the wrapping of the Pade trackers.

  procedure Standard_Track
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                verbose : in boolean; vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Tracks all paths writing output to file,
  --   prompting the user in each step whether to continue or not.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    ans : character;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      new_line(file);
      Wrapped_Pade_Trackers.Run(file,p'last+1,p,ls.v,ls,verbose,vrblvl);
      new_line;
      put("Continue ? "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
    new_line(file);
    put_line(file,"THE SOLUTIONS :");
    put(File,Length_Of(sols),natural32(Head_Of(sols).n),sols);
  end Standard_Track;

  procedure Standard_Track
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

  -- DESCRIPTION :
  --   Tracks all paths without output,
  --   prompting the user in each step whether to continue or not.

    use Standard_Complex_Solutions;

    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    ans : character;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Wrapped_Pade_Trackers.Run(p'last+1,p,ls.v,ls,vrblvl);
      new_line;
      put("Continue ? "); Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Track;

  procedure Standard_Main is

  -- DESCRIPTION :
  --   Reads a natural parameter in double precision.
  --   Tunes the parameters and then calls the path tracker.

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols,xtsols : Standard_Complex_Solutions.Solution_List;
    outfile : file_type;
    ans : character;
    nbequ,nbvar,idxpar,nbpar : integer32;
    verbose : boolean;

    use Standard_Complex_Solutions;
    use Standard_Parameter_Systems;

  begin
    Read_Parameter_Homotopy(lp,xtsols,nbequ,nbvar,nbpar);
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
      new_line;
      put("Track one after the other ? "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Standard_Track(lp.all,xtsols,99);
        sols := Wrapped_Solution_Vectors.Create(xtsols);
      else
        sols := Wrapped_Solution_Vectors.Create(xtsols);
        Wrapped_Pade_Trackers.Run(lp'last+1,lp.all,sols,99);
      end if;
      new_line;
      put_line("THE SOLUTIONS :");
      put(standard_output,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    else
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(outfile);
      Standard_System_and_Solutions_io.put(outfile,lp.all,sols);
      new_line(outfile);
      new_line;
      Wrapped_Pade_Trackers.Set_Parameters(outfile);
      new_line;
      put("Verbose ? (y/n) "); Ask_Yes_or_No(ans);
      verbose := (ans = 'y');
      put("Track one after the other ? "); Ask_Yes_or_No(ans);
      if ans = 'y' then
        Standard_Track(outfile,lp.all,xtsols,verbose,99);
      else
        sols := Wrapped_Solution_Vectors.Create(xtsols);
        Wrapped_Pade_Trackers.Run(outfile,lp'last+1,lp.all,sols,verbose,99);
        new_line(outfile);
        put_line(outfile,"THE SOLUTIONS :");
        put(outfile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      end if;
      close(outfile);
    end if;
  end Standard_Main;

begin
  Standard_Main;
end ts_wrappade;
