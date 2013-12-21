with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Black_Box_Root_Counters;            use Black_Box_Root_Counters;

procedure bablroco ( nt : in natural32; infilename,outfilename : in string ) is

  procedure Read_System ( file : in out file_type; filename : in string;
                          lp : out Link_to_Poly_Sys ) is
  begin
    if filename /= "" then
      Open_Input_File(file,filename);
      get(file,lp);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   lp := null; return;
  end Read_System;

  procedure Main is

    lp,lq : Link_to_Poly_Sys;
    infile,outfile : file_type;
    rc : natural32;
    roco,poco : duration;
    qsols,qsols0 : Solution_List;

  begin
    Read_System(infile,infilename,lp);
    if lp = null then
      new_line;
      get(lp);
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,lp.all);
    lq := new Poly_Sys(lp'range);
    Black_Box_Root_Counting
      (outfile,integer32(nt),lp.all,false,rc,lq.all,qsols,qsols0,roco,poco);
  end Main;

begin
  Main;
end bablroco;
