with Multprec_Complex_Numbers_io;        use Multprec_Complex_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Multprec_Complex_Polynomials;       use Multprec_Complex_Polynomials;

package body Multprec_Complex_Term_Lists_io is

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

end Multprec_Complex_Term_Lists_io;
