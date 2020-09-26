with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;
with Black_Box_Helpers;
with Black_Box_Univariate_Solvers;       use Black_Box_Univariate_Solvers;
with Black_Box_Factorization;            use Black_Box_Factorization;

package body Black_Box_Single_Solvers is

  procedure Solve ( infilename,outfilename : in string;
                    p : in Standard_Complex_Polynomials.Poly;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions;

    n : constant natural32
      := Standard_Complex_Polynomials.Number_of_Unknowns(p);
    sols : Solution_List;
    outfile : file_type;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line("-> in black_box_single_solvers.Solve 1 ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if n = 1 then
      if verbose > 0 then
        put_line("-> calling the method of Weierstrass (Durand-Kerner) ...");
      end if;
      if output_to_file then
        Black_Box_Durand_Kerner(outfile,p,sols);
        Black_Box_Helpers.Append_Solutions_to_Input_file
          (infilename,sols,append_sols);
      else
        Black_Box_Durand_Kerner(p,sols);
        new_line;
        put_line("THE SOLUTIONS :");
        put(Length_Of(sols),1,sols);
      end if;
    elsif n > 1 then
      if verbose > 0 then
        put_line("-> calling the absolute factorization methods ...");
      end if;
      Standard_Black_Box_Factorization(infilename,outfile,p);
    else
      if output_to_file then
        put(outfile,"Number of unknowns = "); put(outfile,n,1);
        put_line(outfile,"...");
      else 
        put("Number of unknowns = "); put(n,1); put_line("...");
      end if;
    end if;
  end Solve;

  procedure Solve ( infilename,outfilename : in string;
                    p : in DoblDobl_Complex_Polynomials.Poly;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions;

    n : constant natural32
      := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p);
    sols : Solution_List;
    outfile : file_type;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line("-> in black_box_single_solvers.Solve 2 ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if n = 1 then
      if verbose > 0 then
        put_line("-> calling the method of Weierstrass (Durand-Kerner) ...");
      end if;
      if output_to_file then
        Black_Box_Durand_Kerner(outfile,p,sols);
        Black_Box_Helpers.Append_Solutions_to_Input_file
          (infilename,sols,append_sols);
      else
        Black_Box_Durand_Kerner(p,sols);
        new_line;
        put_line("THE SOLUTIONS :");
        put(Length_Of(sols),1,sols);
      end if;
    elsif n > 1 then
      if verbose > 0 then
        put_line("-> calling the absolute factorization methods ...");
      end if;
      DoblDobl_Black_Box_Factorization(infilename,outfile,p);
    else
      if output_to_file then
        put(outfile,"Number of unknowns = "); put(outfile,n,1);
        put_line(outfile,"...");
      else 
        put("Number of unknowns = "); put(n,1); put_line("...");
      end if;
    end if;
  end Solve;

  procedure Solve ( infilename,outfilename : in string;
                    p : in QuadDobl_Complex_Polynomials.Poly;
                    append_sols : in boolean; verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions;

    n : constant natural32
      := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p);
    sols : Solution_List;
    outfile : file_type;
    output_to_file : boolean;

  begin
    if verbose > 0 then
      put_line("-> in black_box_single_solvers.Solve 2 ...");
    end if;
    Black_Box_Helpers.Ask_Output_File(outfile,outfilename,output_to_file);
    if n = 1 then
      if verbose > 0 then
        put_line("-> calling the method of Weierstrass (Durand-Kerner) ...");
      end if;
      if output_to_file then
        Black_Box_Durand_Kerner(outfile,p,sols);
        Black_Box_Helpers.Append_Solutions_to_Input_file
          (infilename,sols,append_sols);
      else
        Black_Box_Durand_Kerner(p,sols);
        new_line;
        put_line("THE SOLUTIONS :");
        put(Length_Of(sols),1,sols);
      end if;
    elsif n > 1 then
      if verbose > 0 then
        put_line("-> calling the absolute factorization methods ...");
      end if;
      QuadDobl_Black_Box_Factorization(infilename,outfile,p);
    else
      if output_to_file then
        put(outfile,"Number of unknowns = "); put(outfile,n,1);
        put_line(outfile,"...");
      else 
        put("Number of unknowns = "); put(n,1); put_line("...");
      end if;
    end if;
  end Solve;

end Black_Box_Single_Solvers;
