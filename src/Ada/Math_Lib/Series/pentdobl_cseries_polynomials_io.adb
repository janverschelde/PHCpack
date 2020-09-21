with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with PentDobl_Complex_Polynomials_io;
with PentDobl_Complex_Series;           use PentDobl_Complex_Series;

package body PentDobl_CSeries_Polynomials_io is

  procedure put ( t : in Term ) is
  begin
    put(standard_output,t);
  end put;

  procedure put ( file : in file_type; t : in Term ) is

    s : constant Link_to_Series := t.cf;
    dt : constant natural32 := Standard_Natural_Vectors.Sum(t.dg.all);

  begin
    put(file," + ( ");
    for i in s.cff'range loop
      if i > 0
       then put(file,"   + ");
      end if;
      PentDobl_Complex_Polynomials_io.Write_Number(file,s.cff(i));
      if i > 0 then
        put(file,"*t");
        if i > 1
         then put(file,"^"); put(file,i,1);
        end if;
      end if;
      if i = s.cff'last
       then put(file," )");
       else new_line(file);
      end if;
    end loop;
    if dt > 0 then
      for i in t.dg'range loop
        if t.dg(i) > 0 then
          put(file,"*x"); put(file,i,1);
          if t.dg(i) > 1
           then put(file,"^"); put(file,i,1);
          end if;
        end if;
      end loop;
    end if;
  end put;

  procedure put ( p : in Poly ) is
  begin
    put(standard_output,p);
  end put;

  procedure put ( file : in file_type; p : in Poly ) is

    nbr : constant natural32 := Number_of_Terms(p);
    cnt : natural32 := 0;

    procedure Write_Term ( t : in Term; continue : out boolean ) is
    begin
      cnt := cnt + 1;
      put(file,t);
      if cnt < nbr
       then new_line(file);
      end if;
      continue := true;
    end Write_Term;
    procedure Write_Terms is new Visiting_Iterator(Write_Term);

  begin
    Write_Terms(p);
    put_line(file,";");
  end put;

end PentDobl_CSeries_Polynomials_io;
