with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with QuadDobl_Complex_Polynomials;       use QuadDobl_Complex_Polynomials;

package body QuadDobl_Complex_Term_Lists_io is

  procedure put ( p : in Term_List ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Term_List ) is

    tmp : Term_List := p;
    t : Term;

  begin
    while not Is_Null(tmp) loop
      t := Head_Of(tmp);
      put(t.cf);
      put(t.dg.all);
      new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( p : in Array_of_Term_Lists ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Array_of_Term_Lists ) is
  begin
    put(file,p(p'first));
    for i in p'first+1..p'last loop
      new_line(file);
      put(file,p(i));
    end loop;
  end put;

end QuadDobl_Complex_Term_Lists_io;
