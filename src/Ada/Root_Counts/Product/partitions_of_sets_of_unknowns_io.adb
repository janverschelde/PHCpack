with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Numbers_io;
with Symbol_Table,Symbol_Table_io;
with Sets_of_Unknowns;
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

-- INTERACTIVE INPUT :

  function iget ( m : natural32 ) return Standard_Natural_Vectors.Vector is

    nbr : constant integer32 := integer32(Symbol_Table.Number);
    res : Standard_Natural_Vectors.Vector(1..nbr) := (1..nbr => 0);
    idx : natural32 := 0;

  begin
    for i in 1..nbr loop
      declare
        sb : constant Symbol_Table.Symbol := Symbol_Table.get(natural32(i));
      begin
        loop
          put("-> which set will "); Symbol_Table_io.put(sb);
          put(" be in ? ");
          Numbers_io.Read_Natural(idx);
          exit when ((idx > 0) and (idx <= m));
          put("-> index must be larger than 0 and less than "); 
          put(m+1,1); put_line(".  Please try again.");
        end loop;
        res(i) := idx;
      end;
    end loop;
    return res;
  end iget;

  function Make_Partition
             ( n,m : natural32; p : Standard_Natural_Vectors.Vector )
             return Partition is

    res : Partition(1..m);

  begin
    for i in res'range loop
      res(i) := Sets_of_Unknowns.Create(n); -- initialize each set
    end loop;
    for i in p'range loop
      Sets_of_Unknowns.Add(res(p(i)),natural32(i));
    end loop;
    return res;
  end Make_Partition;

  function iget ( m : natural32 ) return Partition is

    p : constant Standard_Natural_Vectors.Vector := iget(m);
    res : constant Partition(1..m)
        := Make_Partition(natural32(p'last),m,p);

  begin
    return res;
  end iget;
 
end Partitions_of_Sets_of_Unknowns_io;
