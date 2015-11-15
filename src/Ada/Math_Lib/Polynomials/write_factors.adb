with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Symbol_Table_io;

package body Write_Factors is

  procedure Write_Factor ( file : in file_type; d : in natural32;
                           sb : in Symbol; pow : in power ) is
  begin
    Symbol_Table_io.put(file,sb);
    if d > 1 then
      if pow = '^'
       then put(file,'^');
       else put(file,"**");
      end if;
      put(file,d,1);
    end if;
  end Write_Factor;

  procedure Write_Factor ( file : in file_type; d,i : in natural32;
                           standard : in boolean; pow : in power ) is

    sb : Symbol;

  begin
    if standard then
      put(file,'x');
      put(file,i,1);
    else 
      sb := Symbol_Table.get(i);
      Symbol_Table_io.put(file,sb);
    end if;
    if d > 1 then
      if pow = '^'
       then put(file,'^');
       else put(file,"**");
      end if;
      put(file,d,1);
    end if;
  end Write_Factor;

end Write_Factors;
