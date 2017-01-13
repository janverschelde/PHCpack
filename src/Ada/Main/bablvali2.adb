with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with File_Scanning;                      use File_Scanning;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Solutions;         use DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;      use DoblDobl_Complex_Solutions_io;
with Black_Box_Root_Refiners;

procedure bablvali2 ( infilename,outfilename : in string ) is

  lp : Link_to_Poly_Sys;
  sysonfile : boolean;

  procedure Read_System ( file : in out file_type; filename : in string ) is
  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      sysonfile := true;
    else
      sysonfile := false;
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   sysonfile := false;
                   lp := null; return;
  end Read_System;

  procedure Main is

    infile,outfile : file_type;
    nbvar : natural32;
    ans : character;
    sols : Solution_List;
    n : integer32 := 0;
    found : boolean;

  begin
    Read_System(infile,infilename);
    if lp = null then
      new_line;
      put("Is the system on file ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        put_line("Reading the name of the input file.");
        Read_Name_and_Open_File(infile);
        get(infile,lp);
        sysonfile := true;
      else
        put("Give the dimension : "); get(n);
        lp := new Poly_Sys(1..n);
        put("Give "); put(n,1); put(" "); put(n,1); 
        put_line("-variate polynomials :");
       -- get(natural32(n),lp.all);
        get(standard_input,lp.all);
        skip_line;
        sysonfile := false;
      end if;
    end if;
    nbvar := DoblDobl_Complex_Polynomials.Number_of_Unknowns(lp(lp'first));
    Create_Output_File(outfile,outfilename);
    if lp'last = integer32(nbvar)
     then put(outfile,natural32(lp'last),lp.all);
     else put(outfile,natural32(lp'last),nbvar,lp.all);
    end if;
    if sysonfile then
      Scan_and_Skip(infile,"THE SOLUTIONS",found);
      if found
       then get(infile,sols);
      end if;
      Close(infile);
    else
      found := false;
    end if;
    if not found
     then new_line; Read(sols);
    end if;
    Black_Box_Root_Refiners.Refine_Roots(outfile,lp.all,sols);
  end Main;

begin
  Main;
end bablvali2;
