with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with String_Splitters;                  use String_Splitters;
with Communications_with_User;          use Communications_with_User;
with File_Scanning;                     use File_Scanning;
with Greeting_Banners;
with DoblDobl_Complex_Poly_Systems;     use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;     use DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Laur_Poly_Convertors;
with DoblDobl_System_Readers;
with Drivers_to_Cascade_Filtering;      use Drivers_to_Cascade_Filtering;

procedure compsolve2
            ( nt : in natural32; infilename,outfilename : in string ) is

  procedure Write_Greeting ( nbtasks : in natural32 ) is

  -- DESCRIPTION :
  --   Writes the greeting and the number of tasks in nbtasks,
  --   or writes no tasking, followed by the precision.

  begin
    put_line(Greeting_Banners.welcome & ".");
    put("Numerical irreducible decomposition solver");
    if nbtasks = 0
     then put(", no tasking");
     else put(", with "); put(nbtasks,1); put(" tasks");
    end if;
    put_line(", in double precision.");
  end Write_Greeting;

  procedure Main is

  -- DESCRIPTION :
  --   Processing of the input and the output file name arguments.

    infile,outfile : file_type;
    outname : Link_to_String;
    append_sols : boolean := false;
    greeted : boolean := false;
    p : Link_to_Poly_Sys;
    q : Link_to_Laur_Sys;
    tofile : character;

  begin
    DoblDobl_System_Readers.Read_System(infile,infilename,q);
    if q = null then
      greeted := true;
      Write_Greeting(nt);
      new_line; get(q);
      new_line;
    else
      Scan_and_Skip(infile,"SOLUTIONS",append_sols);
      append_sols := not append_sols;
      close(infile);
    end if;
    if outfilename /= "" then
      tofile := 'y';
    else
      if not greeted then
        Write_Greeting(nt); new_line;
        greeted := true;
      end if;
      put("Do you want the output to file ? (y/n) ");
      Ask_Yes_or_No(tofile); -- will be 'y' if yes
    end if;
    if tofile = 'y'
     then Create_Output_File(outfile,outfilename,outname);
    end if;
    if DoblDobl_Laur_Poly_Convertors.Is_Genuine_Laurent(q.all) then
      if tofile = 'y'
       then DoblDobl_Embed_and_Cascade(outfile,outname.all,nt,q.all);
       else DoblDobl_Embed_and_Cascade(nt,q.all,true);
      end if;
    else
      declare
        use DoblDobl_Laur_Poly_Convertors;
        t : constant Poly_Sys(q'range)
          := Positive_Laurent_Polynomial_System(q.all);
      begin
        p := new Poly_Sys'(t);
        if tofile = 'y'
         then DoblDobl_Embed_and_Cascade(outfile,outname.all,nt,p.all);
         else DoblDobl_Embed_and_Cascade(nt,p.all,true);
        end if;
      end;
    end if;
  end Main;

begin
  Main;
end compsolve2;
