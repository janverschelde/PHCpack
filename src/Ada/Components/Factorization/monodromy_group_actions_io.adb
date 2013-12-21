with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;

package body Monodromy_Group_Actions_io is

  procedure put ( ic : in Irreducible_Components; i : in integer32 ) is
  begin
    put(Standard_Output,ic,i);
  end put;

  procedure put ( file : in file_type;
                  ic : in Irreducible_Components; i : in integer32 ) is
  begin
    if not Empty(ic,i) then
      declare
        s : constant Vector := Component(ic,i);
      begin
        put(file,"{");
        for i in s'range loop
          if i > s'first
           then put(file,",");
          end if;
          put(file,s(i),1);
        end loop;
        put(file,"}");
      end;
    end if;
  end put;

  procedure put ( ic : in Irreducible_Components ) is
  begin
    put(Standard_Output,ic);
  end put;

  procedure put ( file : in file_type; ic : in Irreducible_Components ) is
  begin
    for i in 1..Sum_of_Degrees(ic) loop
      put(file,ic,i);
    end loop;
    new_line(file);
  end put;

  procedure put_labels ( ic : in Irreducible_Components ) is
  begin
    put(Standard_Output,ic);
  end put_labels;

  procedure put_labels ( file : in file_type;
                         ic : in Irreducible_Components ) is
  begin
    put(file,Cardinality(ic),1);
    new_line(file);
    for i in 1..Sum_of_Degrees(ic) loop
      if not Empty(ic,i) then
        declare
          s : constant Vector := Component(ic,i);
        begin
          put(file,integer32(s'length),1);
          put(file," : ");
          put(file,s); new_line(file);
        end;
      end if;
    end loop;
  end put_labels;

end Monodromy_Group_Actions_io;
