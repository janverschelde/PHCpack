with Communications_with_User;           use Communications_with_User;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;   use DoblDobl_Complex_Laur_Systems_io;

package body DoblDobl_System_Readers is

  procedure Read_System
               ( file : in out file_type; filename : in string;
                 p : out DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    if filename /= "" then
      Open_Input_File(file,filename);
      get(file,p);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   p := null; return;
  end Read_System;

  procedure Read_System
               ( file : in out file_type; filename : in string;
                 p : out DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys ) is
  begin
    if filename /= "" then
      Open_Input_File(file,filename);
      get(file,p);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   p := null; return;
  end Read_System;

end DoblDobl_System_Readers;
