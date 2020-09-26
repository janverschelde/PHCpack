with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with QuadDobl_Complex_Solutions_io;      use QuadDobl_Complex_Solutions_io;

package body Black_Box_Helpers is

  function Is_Constant_In
              ( p : Standard_Complex_Polynomials.Poly ) return boolean is

    use Standard_Complex_Numbers;

    nvr : constant natural32
        := Standard_Complex_Polynomials.Number_of_Unknowns(p);
    z32 : constant natural32 := natural32(0); -- make sure it is 32-bit
    zdg : Standard_Complex_Polynomials.Degrees
        := new Standard_Natural_Vectors.Vector'(1..integer32(nvr) => z32);
    cff : constant Complex_Number := Standard_Complex_Polynomials.Coeff(p,zdg);

  begin
    Standard_Complex_Polynomials.Clear(zdg);
    if REAL_PART(cff) = 0.0 and IMAG_PART(cff) = 0.0
     then return false;
     else return true;
    end if;
  end Is_Constant_In;

  function Is_Constant_In
              ( p : DoblDobl_Complex_Polynomials.Poly ) return boolean is

    use DoblDobl_Complex_Numbers;

    nvr : constant natural32
        := DoblDobl_Complex_Polynomials.Number_of_Unknowns(p);
    z32 : constant natural32 := natural32(0); -- make sure it is 32-bit
    zdg : DoblDobl_Complex_Polynomials.Degrees
        := new Standard_Natural_Vectors.Vector'(1..integer32(nvr) => z32);
    cff : constant Complex_Number := DoblDobl_Complex_Polynomials.Coeff(p,zdg);
    zero : constant double_double := create(0.0);

  begin
    DoblDobl_Complex_Polynomials.Clear(zdg);
    if REAL_PART(cff) = zero and IMAG_PART(cff) = zero
     then return false;
     else return true;
    end if;
  end Is_Constant_In;

  function Is_Constant_In
              ( p : QuadDobl_Complex_Polynomials.Poly ) return boolean is

    use QuadDobl_Complex_Numbers;

    nvr : constant natural32
        := QuadDobl_Complex_Polynomials.Number_of_Unknowns(p);
    z32 : constant natural32 := natural32(0); -- make sure it is 32-bit
    zdg : QuadDobl_Complex_Polynomials.Degrees
        := new Standard_Natural_Vectors.Vector'(1..integer32(nvr) => z32);
    cff : constant Complex_Number := QuadDobl_Complex_Polynomials.Coeff(p,zdg);
    zero : constant quad_double := create(0.0);

  begin
    QuadDobl_Complex_Polynomials.Clear(zdg);
    if REAL_PART(cff) = zero and IMAG_PART(cff) = zero
     then return false;
     else return true;
    end if;
  end Is_Constant_In;

  function Are_Constants_In
              ( p : Standard_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Constant_In(p(i))
       then return false;
      end if;
    end loop;
    return true;
  --exception
   -- when others => put_line("Are_Constants_In() raised exception.");
   --                raise;
  end Are_Constants_In;

  function Are_Constants_In
              ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Constant_In(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Are_Constants_In;

  function Are_Constants_In
              ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys ) return boolean is
  begin
    for i in p'range loop
      if not Is_Constant_In(p(i))
       then return false;
      end if;
    end loop;
    return true;
  end Are_Constants_In;

  procedure Timing_Summary
              ( file : in file_type;
                roco,hoco,poco,total : in duration ) is

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b1 : constant string :=
     "  |                    TIMING INFORMATION SUMMARY                     |";
    b2 : constant string :=
     "  |   root counts  |  start system  |  continuation  |   total time   |";

  begin
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b0);
    put(file,"  | ");
    print_hms(file,roco); put(file," | ");
    print_hms(file,hoco); put(file," | ");
    print_hms(file,poco); put(file," | ");
    print_hms(file,total); put_line(file," |");
    put_line(file,b0);
  end Timing_Summary;

  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in Standard_Complex_Solutions.Solution_list;
                append_sols : in boolean ) is

    use Standard_Complex_Solutions;

    infile : file_type;

  begin
    if not Is_Null(sols) and append_sols then
      Open_Append_File(infile,infilename);
      new_line(infile);
      put_line(infile,"THE SOLUTIONS :");
      put(infile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Close(infile);
    end if;
  end Append_Solutions_to_Input_File;

  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in DoblDobl_Complex_Solutions.Solution_list;
                append_sols : in boolean ) is

    use DoblDobl_Complex_Solutions;

    infile : file_type;

  begin
    if not Is_Null(sols) and append_sols then
      Open_Append_File(infile,infilename);
      new_line(infile);
      put_line(infile,"THE SOLUTIONS :");
      put(infile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Close(infile);
    end if;
  end Append_Solutions_to_Input_File;

  procedure Append_Solutions_to_Input_File
              ( infilename : in string;
                sols : in QuadDobl_Complex_Solutions.Solution_list;
                append_sols : in boolean ) is

    use QuadDobl_Complex_Solutions;

    infile : file_type;

  begin
    if not Is_Null(sols) and append_sols then
      Open_Append_File(infile,infilename);
      new_line(infile);
      put_line(infile,"THE SOLUTIONS :");
      put(infile,Length_Of(sols),natural32(Head_Of(sols).n),sols);
      Close(infile);
    end if;
  end Append_Solutions_to_Input_File;

  procedure Ask_Output_File
              ( outfile : out file_type; outfilename : in string;
                output_to_file : out boolean ) is

    outnewname : Link_to_String;

  begin
    Ask_Output_File(outfile,outfilename,output_to_file,outnewname);
  end Ask_Output_file;

  procedure Ask_Output_File
              ( outfile : out file_type; outfilename : in string;
                output_to_file : out boolean;
                outnewname : out Link_to_String ) is

    ans : character := 'y';

  begin
    if outfilename = "" then
      new_line;
      put("Do you want the output to file ? (y/n) ");
      Ask_Yes_or_No(ans);
    end if;
    if ans = 'y'
     then Create_Output_File(outfile,outfilename,outnewname);
    end if;
    output_to_file := (ans = 'y');
  end Ask_Output_file;

end Black_Box_Helpers;
