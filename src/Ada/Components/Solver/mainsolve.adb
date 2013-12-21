with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Drivers_to_Eqn_by_Eqn_Solvers;      use Drivers_to_Eqn_by_Eqn_Solvers;

procedure mainsolve ( infilename,outfilename : in string ) is

  procedure Read_System ( name : in string; p : out Link_to_Poly_Sys ) is

    file : file_type;

  begin
    Open(file,in_file,name);
    get(file,p);
  exception
    when others => new_line;
                   put_line("Exception occurred with file " & name & ".");
                   get(p);
  end Read_System;

  procedure Interactive_Create_Output_File
               ( file : in out file_type; name : in string;
                 new_name : out Link_to_String ) is

    procedure Ask_File is
    begin
      new_line;
      put_line("Reading the name of the output file ...");
      Read_Name_and_Create_File(file,new_name);
    end Ask_File;

  begin
    if name = ""
     then Ask_File;
     else Create_Output_File(file,name,new_name); 
    end if;
  end Interactive_Create_Output_File;

  procedure Main is

    file : file_type;
    lp : Link_to_Poly_Sys;
    name : Link_to_String;

  begin
    if infilename /= ""
     then Read_System(infilename,lp);
     else new_line; get(lp);
    end if;
    Interactive_Create_Output_File(file,outfilename,name);
    if name /= null
     then Shuffle_Polynomials_and_Solve(file,name.all,lp.all);
     else Shuffle_Polynomials_and_Solve(file,outfilename,lp.all);
    end if;
  end Main;

begin
  Main;
end mainsolve;
