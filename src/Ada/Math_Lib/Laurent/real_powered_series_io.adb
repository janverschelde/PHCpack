with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Strings_and_Numbers;
with Standard_Write_Numbers;

package body Real_Powered_Series_IO is

  function append_term ( k : integer32; s : string;
                         c : Standard_Complex_Vectors.Vector;
                         p : Standard_Floating_Vectors.Vector;
                         t : character := 't' ) return string is

  -- DESCRIPTION :
  --   Appends the k-th term to the string s, if k <= c'last.

  begin
    if k > c'last then
      return s;
    else
      declare
        cff : constant string := Strings_and_Numbers.Convert(c(k));
        pwr : constant string := Strings_and_Numbers.Convert(p(k));
        trm : constant string := s & " + " & cff & t & "**" & pwr;
      begin
        return append_term(k+1,trm,c,p,t);
      end;
    end if;
  end append_term;

  function to_string ( c : Standard_Complex_Vectors.Vector;
                       p : Standard_Floating_Vectors.Vector;
                       t : character := 't' ) return string is

    cst : constant string := Strings_and_Numbers.Convert(c(0));

  begin
    return append_term(1,cst,c,p,t);
  end to_string;

  procedure put ( c : in Standard_Complex_Vectors.Vector;
                  p : in Standard_Floating_Vectors.Vector;
                  t : in character := 't' ) is
  begin
    put(standard_output,c,p,t);
  end put;

  procedure put ( file : in file_type;
                  c : in Standard_Complex_Vectors.Vector;
                  p : in Standard_Floating_Vectors.Vector;
                  t : in character := 't' ) is

    cnt : natural32 := 0;

  begin
    Standard_Write_Numbers.Write_Number(file,c(0),cnt);
    for i in p'range loop
      put(file," + ");
      Standard_Write_Numbers.Write_Coefficient(file,c(i),cnt);
      put(file,t);
      put(file,"**");
      put(file,p(i),1,14,3);
    end loop;
  end put;

  procedure put_line ( c : in Standard_Complex_Vectors.Vector;
                       p : in Standard_Floating_Vectors.Vector;
                       t : in character := 't' ) is
  begin
    put_line(standard_output,c,p,t);
  end put_line;

  procedure put_line ( file : in file_type;
                       c : in Standard_Complex_Vectors.Vector;
                       p : in Standard_Floating_Vectors.Vector;
                       t : in character := 't' ) is

    cnt : natural32 := 0;
  
  begin
    Standard_Write_Numbers.Write_Number(file,c(0),cnt); new_line(file);
    for i in p'range loop
      put(file," + ");
      Standard_Write_Numbers.Write_Coefficient(file,c(i),cnt);
      put(file,t);
      put(file,"**");
      put(file,p(i),1,14,3); new_line(file);
    end loop;
  end put_line;

end Real_Powered_Series_IO;
