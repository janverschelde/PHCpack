with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Main_Eqn_by_Eqn_Solvers;            use Main_Eqn_by_Eqn_Solvers;

procedure bablsolve
            ( p : in Poly_Sys; outname : in string;
              outfile : in file_type; to_file : in boolean;
              verbose : in integer32 := 0 ) is

  file : file_type;
  name : Link_to_String;

begin
  if verbose > 0
   then put_line("-> in bablsolve for a polynomial system ...");
  end if;
  new_line;
  put_line("Calling the equation-by-equation solver ...");
  if to_file then -- output file is already given in outfile
    Shuffle_Polynomials_and_Solve(outfile,outname,p);
  else
    if outname = "" then
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file,name);
      Shuffle_Polynomials_and_Solve(file,name.all,p);
    else
      Create_Output_File(file,outname,name);
      if name /= null
       then Shuffle_Polynomials_and_Solve(file,name.all,p);
       else Shuffle_Polynomials_and_Solve(file,outname,p);
      end if;
    end if;
  end if;
end bablsolve;
