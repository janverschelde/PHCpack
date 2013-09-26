with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Generic_VecMats_io is

  use Matrices,Matrices_io;

  procedure get ( v : in out VecMat ) is
  begin
    get(Standard_Input,v);
  end get;

  procedure get ( file : in file_type; v : in out VecMat ) is

    n1,n2 : natural32 := 0;

  begin
    for i in v'range loop
      get(file,n1);
      get(file,n2);
      v(i) := new Matrix(1..integer32(n1),1..integer32(n2));
      get(file,v(i).all);
    end loop;
  end get;

  procedure get ( n,n1,n2 : in natural32; v : in out Link_to_VecMat ) is
  begin
    get(Standard_Input,n,n1,n2,v);
  end get;

  procedure get ( file : in file_type; n,n1,n2 : in natural32;
                  v : in out Link_to_VecMat ) is
  begin
    v := new VecMat(1..integer32(n));
    for i in v'range loop
      v(i) := new Matrix(1..integer32(n1),1..integer32(n2));
      get(file,v(i).all);
    end loop;
  end get;

  procedure put ( v : in VecMat ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in VecMat ) is
  begin
    for i in v'range loop
      put(file,v(i).all); new_line(file);
    end loop;
  end put;

  procedure put ( v : in Link_to_VecMat ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_VecMat ) is
  begin
    if v /= null
     then put(file,v.all);
    end if;
  end put;

  procedure put ( v : in VecMat; dp : in natural32 ) is
  begin
    put(Standard_Output,v,dp);
  end put;

  procedure put ( file : in file_type; v : in VecMat; dp : in natural32 ) is
  begin
    for i in v'range loop
      put(file,v(i).all,dp); new_line(file);
    end loop;
  end put;

  procedure put ( v : in Link_to_VecMat; dp : in natural32 ) is
  begin
    put(Standard_Output,v,dp);
  end put;

  procedure put ( file : in file_type;
                  v : in Link_to_VecMat; dp : in natural32 ) is
  begin
    if v /= null
     then put(file,v.all,dp);
    end if;
  end put;

end Generic_VecMats_io;
