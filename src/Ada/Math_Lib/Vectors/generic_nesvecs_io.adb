with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;

package body Generic_NesVecs_io is

  use Vectors_io;

  procedure get ( v : in out NesVec ) is
  begin
    get(Standard_Input,v);
  end get;

  procedure get ( file : in file_type; v : in out NesVec ) is
  begin
    case v.n is
      when 1 => get(file,v.v);
      when others => get(file,v.w);
    end case;
  end get;

  procedure get ( v : in out Link_to_NesVec ) is
  begin
    get(Standard_Input,v);
  end get;

  procedure get ( file : in file_type; v : in out Link_to_NesVec ) is

    n : natural32 := 0;
    a,b : integer32 := 0;

  begin
    get(file,n); get(file,a); get(file,b);
    declare
      v_rep : NesVec(n,a,b);
    begin
      get(file,v_rep);
      v := new NesVec'(v_rep);
    end;
  end get;

  procedure get ( v : in out Array_of_NesVecs ) is
  begin
    get(Standard_Input,v);
  end get;

  procedure get ( file : in file_type; v : in out Array_of_NesVecs ) is
  begin
    for i in v'range loop
      get(file,v(i));
    end loop;
  end get;

  procedure put ( v : in NesVec ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in NesVec ) is
  begin
    put(file,v.n,1);
    put(file," "); put(file,v.a,1);
    put(file," "); put(file,v.b,1);
    new_line(file);
    case v.n is
      when 1 => put_line(file,v.v);
      when others => put(file,v.w);
    end case;
  end put;

  procedure put ( v : in Link_to_NesVec ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_NesVec ) is
  begin
    if v /= null
     then put(file,v.all);
    end if;
  end put;

  procedure put ( v : in Array_of_NesVecs ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Array_of_NesVecs ) is
  begin
    for i in v'range loop
      put(file,v(i));
    end loop;
  end put;

end Generic_NesVecs_io;
