package body Generic_Matrices_io is

  use Ring_io;

  procedure get ( m : in out Matrix ) is
  begin
    get(Standard_Input,m);
  end get;

  procedure get ( m : in out Matrix; rw1,rw2 : in integer32 ) is
  begin
    get(Standard_Input,m,rw1,rw2);
  end get;

  procedure get ( file : in file_type; m : in out Matrix ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        get(file,m(i,j));
      end loop;
    end loop;
  end get;

  procedure get ( file : in file_type;
                  m : in out Matrix; rw1,rw2 : in integer32 ) is
  begin
    for i in rw1..rw2 loop
      for j in m'range(2) loop
        get(file,m(i,j));
      end loop;
    end loop;
  end get;

  procedure put ( m : in Matrix ) is
  begin
    put(Standard_Output,m);
  end put;

  procedure put ( m : in Matrix; rw1,rw2 : in integer32 ) is
  begin
    put(Standard_Output,m,rw1,rw2);
  end put;

  procedure put ( file : in file_type; m : in Matrix ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j));
      end loop;
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  m : in Matrix; rw1,rw2 : in integer32 ) is
  begin
    for i in rw1..rw2 loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j));
      end loop;
      new_line(file);
    end loop;
  end put;

  procedure put ( m : in Matrix; dp : in natural32 ) is
  begin
    put(Standard_Output,m,dp);
  end put;

  procedure put ( m : in Matrix; rw1,rw2 : in integer32; dp : in natural32 ) is
  begin
    put(Standard_Output,m,rw1,rw2,dp);
  end put;

  procedure put ( file : in file_type; m : in Matrix; dp : in natural32 ) is
  begin
    for i in m'range(1) loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j),dp);
      end loop;
      new_line(file);
    end loop;
  end put;

  procedure put ( file : in file_type;
                  m : in Matrix; rw1,rw2 : in integer32; dp : in natural32 ) is
  begin
    for i in rw1..rw2 loop
      for j in m'range(2) loop
        put(file,' '); put(file,m(i,j),dp);
      end loop;
      new_line(file);
    end loop;
  end put;

end Generic_Matrices_io;
