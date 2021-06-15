with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_integer_VecVecs;          use Standard_integer_VecVecs;
with Standard_integer_VecVecs_io;       use Standard_integer_VecVecs_io;

package body Integer_Faces_of_Polytope_io is

  procedure put ( f : in Face ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Face ) is
  begin
    put(file," spanned by "); put(file,integer32(f.all'length),1);
    put_line(file," points :"); put(file,f.all);
  end put;

  procedure put ( f : in Faces ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Faces ) is

    cnt : integer32 := 0;
    tmp : Faces := f;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      put(file,"Face "); put(file,cnt,1); put(file,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( f : in Array_of_Faces ) is
  begin
    put(Standard_Output,f);
  end put;

  procedure put ( file : in file_type; f : in Array_of_Faces ) is
  begin
    for i in f'range loop
      if not Is_Null(f(i)) then
        put(file,"faces at component "); put(file,i,1);
        put_line(file," :"); put(file,f(i));
      end if;
    end loop;
  end put;

end Integer_Faces_of_Polytope_io;
