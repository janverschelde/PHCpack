with Standard_Integer_Numbers_io;
with Standard_Complex_Vectors_io;

package body Standard_Dense_Series2_io is

  procedure get ( s : in out Series ) is
  begin
    get(standard_input,s);
  end get;

  procedure get ( file : in file_type; s : in out Series ) is
  begin
    Standard_Complex_Vectors_io.get(file,s.cff(0..s.deg));
  end get;

  procedure put ( s : in Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type; s : in Series ) is
  begin
    Standard_Complex_Vectors_io.put_line(file,s.cff(0..s.deg));
  end put;

  procedure put ( s : in Series; dp : in natural32 ) is
  begin
    put(standard_output,s,dp);
  end put;

  procedure put ( file : in file_type;
                  s : in Series; dp : in natural32 ) is
  begin
    Standard_Complex_Vectors_io.put_line(file,s.cff(0..s.deg),dp);
  end put;

end Standard_Dense_Series2_io;
