with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors;

package body Brackets_io is

  procedure get ( b : in out Bracket ) is
  begin
    get(Standard_Input,b);
  end get;

  procedure get ( b : in out Bracket; sign : out integer32 ) is
  begin
    get(Standard_Input,b,sign);
  end get;

  procedure get ( file : in file_type; b : in out Bracket ) is

    sign : integer32;

  begin
    get(file,b,sign);
  end get;

  procedure get ( file : in file_type; b : in out Bracket;
                  sign : out integer32 ) is

    v : Standard_Natural_Vectors.Vector(b'range);

  begin
    for i in b'range loop
      v(i) := 0;
      get(file,v(i));
    end loop;
    Create(v,b,sign);
  end get;

  procedure put ( b : in Bracket ) is
  begin
    put(Standard_Output,b);
  end put;

  procedure put ( file : in file_type; b : in Bracket ) is
  begin
    put(file,"[");
    for i in b'first..b'last-1 loop
      put(file,b(i),1);
      put(file," ");
    end loop;
    put(file,b(b'last),1);
    put(file,"]");
  end put;

end Brackets_io;
