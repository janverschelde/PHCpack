with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_Laur_Poly_Convertors;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Black_Box_Root_Refiners;
with Prompt_for_Systems;
with Prompt_for_Solutions;

procedure bablvali ( infilename,outfilename : in string;
                     verbose : in integer32 := 0 ) is

  procedure Refine ( file : in out file_type; lp : in Link_to_Laur_Sys;
                     sysonfile : in boolean ) is

  -- DESCRIPTION :
  --   Calls the root refiner on the system lp,
  --   reading the solutions from file if sysonfile.

    outfile : file_type;
    sols : Solution_List;
    nbvar : constant natural32
          := Standard_Complex_Laurentials.Number_of_Unknowns(lp(lp'first));

  begin
    Create_Output_File(outfile,outfilename);
    if lp'last = integer32(nbvar)
     then put(outfile,natural32(lp'last),lp.all);
     else put(outfile,natural32(lp'last),nbvar,lp.all);
    end if;
    Prompt_for_Solutions.Read_Solutions(file,sysonfile,sols);
    if Standard_Laur_Poly_Convertors.Is_Genuine_Laurent(lp.all) then
      Black_Box_Root_Refiners.Refine_Roots(outfile,lp.all,sols);
    else
      declare
        use Standard_Laur_Poly_Convertors;
        p : constant Poly_Sys(lp'range)
          := Positive_Laurent_Polynomial_System(lp.all);
      begin
        Black_Box_Root_Refiners.Refine_Roots(outfile,p,sols);
      end;
    end if;
  end Refine;

  procedure Main is

  -- DESCRIPTION :
  --   Reads the system, the solutions,
  --   and then calls the black box root refiner.

    infile : file_type;
    sysonfile : boolean;
    lp : Link_to_Laur_Sys := null;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in bablvali.Main ...");
    end if;
    Prompt_for_Systems.Read_System(infile,infilename,lp,sysonfile);
    if lp /= null
     then Refine(infile,lp,sysonfile);
    end if;
  end Main;

begin
  Main;
end bablvali;
