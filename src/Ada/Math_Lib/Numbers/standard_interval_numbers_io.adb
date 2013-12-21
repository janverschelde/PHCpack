with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;

package body Standard_Interval_Numbers_io is

-- AUXILIARY :

  procedure Skip_Spaces ( file : in file_type; ch : in out character ) is

  -- DESCRIPTION :
  --   As long as ch is a space, a new character ch will be read from file.

  begin
    while ch = ' ' loop
      get(file,ch);
    end loop;
  end Skip_Spaces;

-- INPUT OPERATIONS :

  procedure get ( i : in out Interval ) is
  begin
    get(standard_input,i);
  end get;

  procedure get ( file : file_type; i : in out Interval) is

    ch : character;
    a,b : double_float := 0.0;

  begin
    get(file,ch);
    Skip_Spaces(file,ch);      -- read opening bracket
    if ch = '[' then
      get(file,a);
      get(file,ch);
      Skip_Spaces(file,ch);
      if ch = ',' then
        get(file,b);
      end if;
      Skip_Spaces(file,ch);    -- also read closing bracket
    end if;
    i := Create(a,b);
  end get;

-- OUTPUT OPERATIONS :

  procedure put ( i : in Interval ) is
  begin
    put(standard_output,i);
  end put;

  procedure put ( i : in Interval; d : in natural32 ) is
  begin
    put(standard_output,i,d);
  end put;

  procedure put ( file : in file_type; i : in Interval ) is
  begin
    put(file,'[');
    put(file,left(i)); put(file,','); put(file,right(i));
    put(file,']');
  end put;

  procedure put ( file : in file_type; i : in Interval; d : in natural32 ) is
  begin
    put(file,'[');
    put(file,left(i),d); put(file,','); put(file,right(i),d);
    put(file,']');
  end put;

end Standard_Interval_Numbers_io;
