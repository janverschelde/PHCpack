with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with String_Splitters;                  use String_Splitters;
with Characters_and_Numbers;            use Characters_and_Numbers;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;

procedure ts_rdisol is

-- DESCRIPTION :
--   Recursive application of the intrinsic diagonal homotopies
--   to solve a system of polynomial equation.

  procedure Initialize ( file : in file_type; name : in string;
                         nq,nv : in integer32; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Creates witness sets for all polynomials in p
  --   and writes each witness set to a separate file.

  -- ON ENTRY :
  --   file     main output file type;
  --   name     name of the file for the input system;
  --   nq       number of equations in the system p;
  --   nv       number of variables in the system p;
  --   p        the polynomials defining the system.

    d : constant integer32 := nv - 1;

  begin
    for i in 1..nq loop
      put_line("Writing " & name & "_" & Convert(i) & "e" & convert(d));
    end loop;
  end Initialize;

  procedure Main is

    infile,outfile : file_type;
    name : Link_to_String;
    lp : Link_to_Poly_Sys;
    nq,nv : integer32 := 0;

  begin
    new_line;
    put_line("recursive solving with intrinsic diagonal homotopies");
    new_line;
    put_line("Reading the name of the input file.");
    Read_Name_and_Open_File(infile,name);
    get(infile,lp);
    declare
      outname : constant string := name.all & ".out";
    begin
      Create_Output_File(outfile,outname);
      nq := lp'last;
      nv := integer32(Number_of_Unknowns(lp(lp'last)));
      put(outfile,natural32(nq),natural32(nv),lp.all);
      new_line;
      put_line("See the file " & outname & " for results ...");
      new_line;
    end;
    Initialize(outfile,name.all,nq,nv,lp.all);
  end Main;

begin
  Main;
end ts_rdisol;
