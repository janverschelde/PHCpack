with Symbol_Table,Symbol_Table_io;       use Symbol_Table;
with Set_Structure;

package body Set_Structure_io is

  procedure get is
  begin
    Set_Structure_io.get(Standard_Input);
  end get;

  procedure get ( file : in file_type ) is

    ch : character;
    sb : symbol;

  begin
    for i in 1..Set_Structure.Dimension loop
     -- Read sets for ith equation :
      for j in 1..Set_Structure.Number_of_Sets(i) loop
       -- Read the jth set of the ith equation :
        get(file,ch); while (ch = ' ') loop get(file,ch); end loop;
	if ch = '{' then
	  loop
            get(file,ch);
            while (ch = ' ') loop get(file,ch); end loop;
              exit when (ch = '}');
              sb := (sb'range => ' ');
              Symbol_Table_io.get(file,ch,sb,'}');
              Set_Structure.Add(i,j,Symbol_Table.Get(sb));
              exit when (ch = '}');
            end loop;
        end if;
      end loop;
    end loop;
  end get;

  procedure put ( i,j : in natural32 ) is
  begin
    Set_Structure_io.put(Standard_Output,i,j);
  end put;

  procedure put ( file : in file_type; i,j : in natural32 ) is

    n : constant natural32 := Set_Structure.Dimension;
    sb : Symbol;

  begin
    put(file,"{ ");
    for k in 1..n loop
      if Set_Structure.Is_In(i,j,k)
       then sb := Symbol_Table.Get(k);
            Symbol_Table_io.put(file,sb);
            put(file,' ');
      end if;
    end loop;
    put(file,'}');
  end put;

  procedure put ( i : in natural32 ) is
  begin
    Set_Structure_io.put(Standard_Output,i);
  end put;

  procedure put ( file : in file_type; i : in natural32 ) is
  begin
    for j in 1..Set_Structure.Number_of_Sets(i) loop
      Set_Structure_io.put(file,i,j);
    end loop;
  end put;

  procedure put is
  begin
    Set_Structure_io.put(Standard_Output);
  end put;

  procedure put ( file : in file_type ) is

    n : constant natural32 := Set_Structure.Dimension;

  begin
    for i in 1..n loop
      put(file,"     ");
      Set_Structure_io.put(file,i);
      new_line(file);
    end loop;
  end put;

end Set_Structure_io;
