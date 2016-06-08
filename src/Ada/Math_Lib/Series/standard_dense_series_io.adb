with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;

package body Standard_Dense_Series_io is

  procedure put ( s : in Series ) is
  begin
    put(standard_output,s);
  end put;

  procedure put ( file : in file_type; s : in Series ) is
  begin
    Standard_Complex_Vectors_io.put_line
      (file,Standard_Complex_Vectors.Vector(s));
  end put;

end Standard_Dense_Series_io;
