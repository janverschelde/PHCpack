with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;

package body Normal_Cone_Intersections_io is

  procedure put ( ima : in Intersection_Matrix ) is
  begin
    put(Standard_Output,ima);
  end put;

  procedure put ( file : in file_type; ima : in Intersection_Matrix ) is

    cnt : integer32;

  begin
    put(file,"offset vector : "); put(file,ima.sv); new_line(file);
    put_line(file,"intersection matrix : ");
    for i in ima.im'range(1) loop
      cnt := ima.sv'first+1;
      for j in ima.im'range(2) loop
        if j = ima.sv(cnt) then
          put(file," &");
          if cnt < ima.sv'last
           then cnt := cnt + 1;
          end if;
        end if;
        put(file," "); put(file,ima.im(i,j),1);
      end loop;
      new_line(file);
    end loop;
  end put;

end Normal_Cone_Intersections_io;
