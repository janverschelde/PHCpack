package body Standard_Natural_Numbers_io is

  package natural32_io is new text_io.integer_io(natural32);
  package natural64_io is new text_io.integer_io(natural64);

  procedure get ( n : in out natural32 ) is
  begin
    natural32_io.get(n);
  end get;

  procedure get ( file : in file_type; n : in out natural32 ) is
  begin
    natural32_io.get(file,n);
  end get;

  procedure get ( n : in out natural64 ) is
  begin
    natural64_io.get(n);
  end get;

  procedure get ( file : in file_type; n : in out natural64 ) is
  begin
    natural64_io.get(file,n);
  end get;

  procedure put ( n : in natural32 ) is
  begin
    natural32_io.put(n,1);
  end put;

  procedure put ( file : in file_type; n : in natural32 ) is
  begin
    natural32_io.put(file,n,1);
  end put;

  procedure put ( n : in natural64 ) is
  begin
    natural64_io.put(n,1);
  end put;

  procedure put ( file : in file_type; n : in natural64 ) is
  begin
    natural64_io.put(file,n,1);
  end put;

  procedure put ( n,dp : in natural32 ) is
  begin
    natural32_io.put(n,natural(dp));
  end put;

  procedure put ( file : in file_type; n,dp : in natural32 ) is
  begin
    natural32_io.put(file,n,natural(dp));
  end put;

  procedure put ( n : in natural64; dp : in natural32 ) is
  begin
    natural64_io.put(n,natural(dp));
  end put;

  procedure put ( file : in file_type; n : in natural64; dp : in natural32 ) is
  begin
    natural64_io.put(file,n,natural(dp));
  end put;

end Standard_Natural_Numbers_io;
