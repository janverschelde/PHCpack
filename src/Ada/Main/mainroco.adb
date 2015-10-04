with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Drivers_for_Root_Counts;            use Drivers_for_Root_Counts;
with Write_Seed_Number;
with Greeting_Banners;
with Bye_Bye_Message;

procedure mainroco ( infilename,outfilename : in string ) is
 
  n : integer32 := 0;
  outft : file_type;
  lp : Link_to_Poly_Sys;

  procedure Read_System ( filename : in string ) is

    file : file_type;

  begin
    if filename /= "" then
      Open(file,in_file,filename);
      get(file,n);
      lp := new Poly_Sys(1..n);
      get(file,natural32(n),lp.all);
      Close(file);
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null; return;
  end Read_System;

  procedure Main is

  begin
    Read_System(infilename);
    if lp = null
     then new_line; get(lp);
    end if;
    Create_Output_File(outft,outfilename);
    put(outft,lp.all);
    declare
      q : Poly_Sys(lp'range);
      qsols : Solution_List;
      rc : natural32;
    begin
      Driver_for_Root_Counts(outft,lp.all,q,false,qsols,rc);
    end; 
    new_line(outft);
    put_line(outft,Bye_Bye_Message);
    Write_Seed_Number(outft);
    put_line(outft,Greeting_Banners.Version);
    Close(outft);
  end Main;

begin
  Main;
end mainroco;
