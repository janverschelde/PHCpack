with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Drivers_for_Reduction;              use Drivers_for_Reduction;

procedure mainred ( infilename,outfilename : in string ) is

  lp : Link_to_Poly_Sys;

  procedure Read_System ( filename : in string ) is

    file : file_type;

  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,lp);
      Close(file);
    end if;
  exception
    when others => 
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null; return;
  end Read_System;

  procedure Main is

    outfile : file_type;
    d : natural32;
    ans : character;

  begin
    Read_System(infilename);
    if lp = null
     then new_line; get(lp);
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,natural32(lp'last),lp.all); new_line(outfile);
    Driver_for_Reduction(outfile,lp.all,d,false);
    Close(outfile);
    new_line;
    put("Do you want the reduced system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        redfile : file_type;
      begin
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(redfile);
        put(redfile,natural32(lp'last),lp.all);
        Close(redfile);
      end;
    end if;
  end Main;

begin
  Main;
end mainred;
