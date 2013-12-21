with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Power_Transformations;     use Standard_Power_Transformations;

package body Standard_Puiseux_Certificates_io is

  procedure put ( g : in Germ ) is
  begin
    put(standard_output,g);
  end put;

  procedure put ( file : in file_type; g : in Germ ) is

    t : constant Term := g.t;
    k : constant integer32 := Pivot(t.dg.all);
    m : constant Matrix(1..2,1..2) := Eliminate(t.dg.all,k);

  begin
    if k = 2 then
      put(file,"[ x = "); put(file,t.cf); new_line(file);
      put(file,"    + ("); put(file,g.c);
      put(file,"*i)*t^"); put(file,g.w,1);
      put(file,", y = t^"); put(file,t.dg(2),1); put_line(file," ]");
    else
      put(file,"[ x = t^"); put(file,t.dg(1),1);
      if m(2,1) = 0 then
        put(file,",");
      else
        put(file,"*{ ("); put(file,t.cf); put_line(file,"*i)");
        put(file,"    +("); put(file,g.c); put(file,"*i)*t^");
        put(file,g.w,1); put(file," }^"); put(file,m(2,1),1);
        if m(2,2) = 0
         then put(file,",");
         else put_line(file,",");
        end if;
      end if;
      put(file," y = t^"); put(file,t.dg(2),1);
      if m(2,2) /= 0 then
        put(file,"*{ ("); put(file,t.cf); put_line(file,"*i)");
        put(file,"    +("); put(file,g.c);  put(file,"*i)*t^");
        put(file,g.w,1); put(file," }^"); put(file,m(2,2),1);
      end if;
      put_line(file," ]");
    end if;
  end put;

  procedure put ( g : in List_of_Germs ) is
  begin
    put(standard_output,g);
  end put;

  procedure put ( file : in file_type; g : in List_of_Germs ) is

    p : List_of_Germs := g;

  begin
    while not Is_Null(p) loop
      put(file,Head_Of(p));
      p := Tail_Of(p);
    end loop;
  end put;

end Standard_Puiseux_Certificates_io;
