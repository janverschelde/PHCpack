package body Standard_Integer_Numbers_io is

  package integer64_io is new text_io.integer_io(integer64);
  package integer32_io is new text_io.integer_io(integer32);

  procedure get ( i : in out integer32 ) is
  begin
    integer32_io.get(i);
  end get;

  procedure get ( file : in file_type; i : in out integer32 ) is
  begin
    integer32_io.get(file,i);
  end get;

  procedure get ( i : in out integer64 ) is
  begin
    integer64_io.get(i);
  end get;

  procedure get ( file : in file_type; i : in out integer64 ) is
  begin
    integer64_io.get(file,i);
  end get;

  procedure put ( i : in integer32 ) is
  begin
    integer32_io.put(i,1);
  end put;

  procedure put ( file : in file_type; i : in integer32 ) is
  begin
    integer32_io.put(file,i,1);
  end put;

  procedure put ( i : in integer64 ) is
  begin
    integer64_io.put(i,1);
  end put;

  procedure put ( file : in file_type; i : in integer64 ) is
  begin
    integer64_io.put(file,i,1);
  end put;

  procedure put ( i : in integer32; dp : in natural32 ) is
  begin
    integer32_io.put(i,natural(dp));
  end put;

  procedure put ( file : in file_type; i : in integer32; dp : in natural32 ) is
  begin
    integer32_io.put(file,i,natural(dp));
  end put;

  procedure put ( file : in file_type; i : in integer64; dp : in natural32 ) is
  begin
    integer64_io.put(file,i,natural(dp));
  end put;

  procedure put ( i : in integer64; dp : in natural32 ) is
  begin
    integer64_io.put(i,natural(dp));
  end put;

end Standard_Integer_Numbers_io;
