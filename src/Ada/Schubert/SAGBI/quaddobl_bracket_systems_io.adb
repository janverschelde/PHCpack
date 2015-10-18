with QuadDobl_Bracket_Polynomials_io; use QuadDobl_Bracket_Polynomials_io;

package body QuadDobl_Bracket_Systems_io is

  procedure put ( s : in Bracket_System ) is
  begin
    put(Standard_Output,s);
  end put;

  procedure put ( file : in file_type; s : in Bracket_System ) is
  begin
    for i in s'range loop
      put(file,s(i));
    end loop;
  end put;

  procedure put_line ( s : in Bracket_System ) is
  begin
    put_line(Standard_Output,s);
  end put_line;

  procedure put_line ( file : in file_type; s : in Bracket_System ) is
  begin
    for i in s'range loop
      put_line(file,s(i));
    end loop;
  end put_line;

end QuadDobl_Bracket_Systems_io;
