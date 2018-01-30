with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Boolean_Numbers_io is

  package integer32_io is new text_io.integer_io(integer32);

  procedure get ( b : in out boolean ) is

    i : integer32;

  begin
    integer32_io.get(i);
    b := (i = 1);
  end get;

  procedure get ( file : in file_type; b : in out boolean ) is

    i : integer32;

  begin
    integer32_io.get(file,i);
    b := (i = 1);
  end get;

  procedure put ( b : in boolean ) is

    i : integer32;

  begin
    if b
     then i := 1;
     else i := 0;
    end if;
    integer32_io.put(i,1);
  end put;

  procedure put ( file : in file_type; b : in boolean ) is

    i : integer32;

  begin
    if b
     then i := 1;
     else i := 0;
    end if;
    integer32_io.put(file,i,1);
  end put;

  procedure put ( b : in boolean; dp : in natural32 ) is

    i : integer32;

  begin
    if b
     then i := 1;
     else i := 0;
    end if;
    integer32_io.put(i,natural(dp));
  end put;

  procedure put ( file : in file_type; b : in boolean; dp : in natural32 ) is

    i : integer32;

  begin
    if b
     then i := 1;
     else i := 0;
    end if;
    integer32_io.put(file,i,natural(dp));
  end put;

end Boolean_Numbers_io;
