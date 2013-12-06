with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package body Generic_VecVecs_io is

  use Vectors,Vectors_io;

  procedure get ( n : in natural32; v : in out VecVec ) is
  begin
    get(Standard_Input,n,v);
  end get;

  procedure get ( file : in file_type; n : in natural32; v : in out VecVec ) is
  begin
    for i in v'range loop
      get(file,n,v(i));
    end loop;
  end get;

  procedure get ( n1,n2 : in natural32; v : in out Link_to_VecVec ) is
  begin
    get(Standard_Input,n1,n2,v);
  end get;

  procedure get ( file : in file_type; n1,n2 : in natural32;
                  v : in out Link_to_VecVec ) is
  begin
    v := new VecVec(1..integer32(n1));
    get(file,n2,v.all);
  end get;

  procedure put ( v : in VecVec ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in VecVec ) is
  begin
    for i in v'range loop
      put(file,v(i)); new_line(file);
    end loop;
  end put;

  procedure put ( v : in Link_to_VecVec ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_VecVec ) is
  begin
    if v /= null
     then put(file,v.all);
    end if;
  end put;

  procedure put_line ( v : in VecVec ) is
  begin
    put_line(Standard_Output,v);
  end put_line;

  procedure put_line ( file : in file_type; v : in VecVec ) is
  begin
    for i in v'range loop
      put_line(file,v(i)); new_line(file);
    end loop;
  end put_line;

  procedure put_line ( v : in Link_to_VecVec ) is
  begin
    put_line(Standard_Output,v);
  end put_line;

  procedure put_line ( file : in file_type; v : in Link_to_VecVec ) is
  begin
    if v /= null
     then put_line(file,v.all);
    end if;
  end put_line;

  procedure put ( v : in VecVec; dp : in natural32 ) is
  begin
    put(Standard_Output,v,dp);
  end put;

  procedure put ( file : in file_type; v : in VecVec; dp : in natural32 ) is
  begin
    for i in v'range loop
      put(file,v(i),dp); new_line(file);
    end loop;
  end put;

  procedure put ( v : in Link_to_VecVec; dp : in natural32 ) is
  begin
    put(Standard_Output,v,dp);
  end put;

  procedure put ( file : in file_type;
                  v : in Link_to_VecVec; dp : in natural32 ) is
  begin
    if v /= null
     then put(file,v.all,dp);
    end if;
  end put;

  procedure put_line ( v : in VecVec; dp : in natural32 ) is
  begin
    put_line(Standard_Output,v,dp);
  end put_line;

  procedure put_line ( file : in file_type; v : in VecVec;
                       dp : in natural32 ) is
  begin
    for i in v'range loop
      put_line(file,v(i),dp); new_line(file);
    end loop;
  end put_line;

  procedure put_line ( v : in Link_to_VecVec; dp : in natural32 ) is
  begin
    put_line(Standard_Output,v,dp);
  end put_line;

  procedure put_line ( file : in file_type;
                       v : in Link_to_VecVec; dp : in natural32 ) is
  begin
    if v /= null
     then put_line(file,v.all,dp);
    end if;
  end put_line;

end Generic_VecVecs_io;
