with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with DoblDobl_Complex_Polynomials;       use DoblDobl_Complex_Polynomials;

package body DoblDobl_Complex_Term_Lists_io is

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

end DoblDobl_Complex_Term_Lists_io;
