with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Communications_with_User;           use Communications_with_User;
with DoblDobl_Complex_Poly_Systems;      use DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_System_Readers;
with Drivers_for_Reduction;              use Drivers_for_Reduction;

procedure mainred2 ( infilename,outfilename : in string;
                     verbose : in integer32 := 0 ) is

  procedure Main is

    lp : Link_to_Poly_Sys;
    infile,outfile : file_type;
    d : natural32;
    ans : character;

  begin
    if verbose > 0 then
      put("At verbose level "); put(verbose,1);
      put_line(", in mainred2.Main ...");
    end if;
    DoblDobl_System_Readers.Read_System(infile,infilename,lp);
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
end mainred2;
