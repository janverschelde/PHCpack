with text_io;                            use text_io;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Sets_of_Unknowns_io;                use Sets_of_Unknowns_io;

package body Degree_Sets_Tables_io is

  procedure put ( ase : in Array_of_Sets ) is
  begin
    put(standard_output,ase);
  end put;

  procedure put ( file : in file_type; ase : in Array_of_Sets ) is
  begin
    for i in ase'range loop
      put(file,ase(i)); new_line(file);
    end loop;
  end put;

  procedure put ( dst : in Degree_Sets_Table ) is
  begin
    put(standard_output,dst);
  end put;

  procedure put ( file : in file_type; dst : in Degree_Sets_Table ) is
  begin
    put(file,dst.s);
    put(file,dst.a);
  end put;

end Degree_Sets_Tables_io;
