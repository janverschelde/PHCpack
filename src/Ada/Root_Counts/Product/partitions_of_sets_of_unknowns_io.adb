with Sets_of_Unknowns_io;                use Sets_of_Unknowns_io;

package body Partitions_of_Sets_of_Unknowns_io is

-- DESCRIPTION :
--   This package provides i/o operations for partitions of
--   set of unknowns.

  procedure get ( p : in out Partition ) is
  begin
    get(Standard_Input,p);
  end get;

  procedure get ( file : in file_type; p : in out Partition ) is
  begin
    for i in p'range loop
      get(file,p(i));
    end loop;
  end get;

  procedure put ( p : in Partition ) is
  begin
    put(Standard_Output,p);
  end put;

  procedure put ( file : in file_type; p : in Partition ) is
  begin
    for i in p'range loop
      put(file,p(i));
    end loop;
  end put;
 
end Partitions_of_Sets_of_Unknowns_io;
