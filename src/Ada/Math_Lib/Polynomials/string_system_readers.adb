with Communications_with_User;           use Communications_with_User;

package body String_System_Readers is

  procedure Read_System
               ( file : in out file_type; name : in string;
                 n,m : out natural32; p : out Link_to_Array_of_Strings ) is
  begin
    if name /= "" then
      Open_Input_File(file,name);
      get(file,natural(n),natural(m),p);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   p := null; return;
  end Read_System;

end String_System_Readers;
