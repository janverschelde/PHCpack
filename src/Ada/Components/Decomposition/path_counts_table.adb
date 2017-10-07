with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;

package body Path_Counts_Table is

  procedure Update_Path_Counts
              ( cnts : in out Standard_Natural_VecVecs.VecVec;
                dim,nsols,nsols0,nsols1 : in natural32 ) is

    counts : Standard_Natural_Vectors.Vector(1..3);

  begin
    counts(1) := nsols;
    counts(2) := nsols0;
    counts(3) := nsols1;
    cnts(integer32(dim)) := new Standard_Natural_Vectors.Vector'(counts);
  end Update_Path_Counts;

  procedure Write_Path_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec ) is
  begin
    new_line(file);
    put(file,"dim | ");
    put(file," #paths | ");
    put(file,"slack=0 | ");
    put(file,"slack!=0");
    new_line(file);
    put_line(file,"----+---------+---------+---------");
    for i in reverse cnts'range loop
      put(file,i,3);
      put(file," | "); put(file,cnts(i)(1),7);
      put(file," | "); put(file,cnts(i)(2),7);
      put(file," | "); put(file,cnts(i)(3),7);
      new_line(file);
    end loop;
  end Write_Path_Counts;

  procedure Write_Path_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec;
                times : in Array_of_Duration; totaltime : in duration ) is

    sum_path,sum_slack0,sum_nonsol : natural32 := 0;

  begin
    new_line(file);
    put(file,"dim | ");
    put(file," #paths | ");
    put(file,"slack=0 | ");
    put(file,"slack!=0 | ");
    put(file," CPU user time");
    new_line(file);
    put_line(file,"----+---------+---------+----------+---------------");
    for i in reverse cnts'range loop
      put(file,i,3);
      put(file," | "); put(file,cnts(i)(1),7);
      put(file," | "); put(file,cnts(i)(2),7);
      put(file," | "); put(file,cnts(i)(3),8);
      put(file," | "); print_hms(file,times(integer(i)));
      new_line(file);
      sum_path := sum_path + cnts(i)(1);
      sum_slack0 := sum_slack0 + cnts(i)(2);
      sum_nonsol := sum_nonsol + cnts(i)(3);
    end loop;
    put_line(file,"----+---------+---------+----------+---------------");
    put(file,"sum | "); put(file,sum_path,7);
    put(file," | "); put(file,sum_slack0,7);
    put(file," | "); put(file,sum_nonsol,8);
    put(file," | "); print_hms(file,totaltime);
    new_line(file);
  end Write_Path_Counts;

end Path_Counts_Table;
