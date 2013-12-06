with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package body Generic_Vectors_io is

  use Ring_io;

  procedure get ( v : in out Vector ) is
  begin
    get(Standard_Input,v);
  end get;

  procedure get ( file : in file_type; v : in out Vector ) is
  begin
    for i in v'range loop
      get(file,v(i));
    end loop;
  end get;

  procedure get ( n : in natural32; v : in out Link_to_Vector ) is
  begin
    get(Standard_Input,n,v);
  end get;

  procedure get ( file : in file_type; n : in natural32;
                  v : in out Link_to_Vector ) is
  begin
    v := new Vector(1..integer32(n));
    for i in 1..integer32(n) loop
      get(file,v(i));
    end loop;
  end get;

  procedure put ( v : in Vector ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Vector ) is
  begin
    for i in v'range loop
      put(file,' '); put(file,v(i));
    end loop;
  end put;

  procedure put ( v : in Link_to_Vector ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_Vector ) is
  begin
    if v /= null
     then put(file,v.all);
    end if;
  end put;

  procedure put_line ( v : in Vector ) is
  begin
    put_line(Standard_Output,v);
  end put_line;

  procedure put_line ( file : in file_type; v : in Vector ) is
  begin
    for i in v'range loop
      put(file,v(i)); new_line(file);
    end loop;
  end put_line;

  procedure put_line ( v : in Link_to_Vector ) is
  begin
    put_line(Standard_Output,v);
  end put_line;

  procedure put_line ( file : in file_type; v : in Link_to_Vector ) is
  begin
    if v /= null
     then put_line(file,v.all);
    end if;
  end put_line;

  procedure put ( v : in Vector; dp : in natural32 ) is
  begin
    put(Standard_Output,v,dp);
  end put;

  procedure put ( file : in file_type; v : in Vector; dp : in natural32 ) is
  begin
    for i in v'range loop
      put(file,' '); put(file,v(i),dp);
    end loop;
  end put;

  procedure put ( v : in Link_to_Vector; dp : in natural32 ) is
  begin
    put(Standard_Output,v,dp);
  end put;

  procedure put ( file : in file_type;
                  v : in Link_to_Vector; dp : in natural32 ) is
  begin
    if v /= null
     then put(file,v.all,dp);
    end if;
  end put;

  procedure put_line ( v : in Vector; dp : in natural32 ) is
  begin
    put_line(Standard_Output,v,dp);
  end put_line;

  procedure put_line ( file : in file_type; v : in Vector;
                       dp : in natural32 ) is
  begin
    for i in v'range loop
      put(file,v(i),dp); new_line(file);
    end loop;
  end put_line;

  procedure put_line ( v : in Link_to_Vector; dp : in natural32 ) is
  begin
    put_line(Standard_Output,v,dp);
  end put_line;

  procedure put_line ( file : in file_type;
                       v : in Link_to_Vector; dp : in natural32 ) is
  begin
    if v /= null
     then put_line(file,v.all,dp);
    end if;
  end put_line;

end Generic_Vectors_io;
