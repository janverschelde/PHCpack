with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;

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

    use Standard_Natural_Vectors;

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
      if cnts(i) = null then
        put(file," | "); put(file,natural32(0),7);
        put(file," | "); put(file,natural32(0),7);
        put(file," | "); put(file,natural32(0),7);
        new_line(file);
      else
        put(file," | "); put(file,cnts(i)(1),7);
        put(file," | "); put(file,cnts(i)(2),7);
        put(file," | "); put(file,cnts(i)(3),7);
        new_line(file);
      end if;
    end loop;
  end Write_Path_Counts;

  procedure Write_Path_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec;
                times : in Array_of_Duration; totaltime : in duration ) is

    sum_path,sum_slack0,sum_nonsol : natural32 := 0;

    use Standard_Natural_Vectors;

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
      if cnts(i) = null then
        put(file," | "); put(file,natural32(0),7);
        put(file," | "); put(file,natural32(0),7);
        put(file," | "); put(file,natural32(0),8);
        put(file," | "); print_hms(file,times(integer(i)));
        new_line(file);
      else
        put(file," | "); put(file,cnts(i)(1),7);
        put(file," | "); put(file,cnts(i)(2),7);
        put(file," | "); put(file,cnts(i)(3),8);
        put(file," | "); print_hms(file,times(integer(i)));
        new_line(file);
        sum_path := sum_path + cnts(i)(1);
        sum_slack0 := sum_slack0 + cnts(i)(2);
        sum_nonsol := sum_nonsol + cnts(i)(3);
      end if;
    end loop;
    put_line(file,"----+---------+---------+----------+---------------");
    put(file,"sum | "); put(file,sum_path,7);
    put(file," | "); put(file,sum_slack0,7);
    put(file," | "); put(file,sum_nonsol,8);
    put(file," | "); print_hms(file,totaltime);
    new_line(file);
  end Write_Path_Counts;

  procedure Write_Filter_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec ) is
  begin
    new_line(file);
    new_line(file);
    put(file,"dim | ");
    put(file," solutions after filter");
    new_line(file);
    put_line(file,"----+-----------------------");
    for i in reverse cnts'range loop
      put(file,i,3);
      put(file," | ");
      put(file,cnts(i)(0),1);
      for j in 1..cnts(i)'last loop
        put(file," -> ");
        put(file,cnts(i)(j),1);
      end loop;
      new_line(file);
    end loop;
  end Write_Filter_Counts;

  procedure Write_Filter_Counts
              ( file : in file_type;
                cnts : in Standard_Natural_VecVecs.VecVec;
                times : in Array_of_Duration; totaltime : in duration ) is
  begin
    new_line(file);
    new_line(file);
    put(file,"dim | ");
    put(file," CPU user time |");
    put(file," solutions after filter");
    new_line(file);
    put_line(file,"----+----------------+-----------------------");
    for i in reverse cnts'range loop
      put(file,i,3);
      put(file," | "); print_hms(file,times(integer(i)));
      put(file," | ");
      put(file,cnts(i)(0),1);
      for j in 1..cnts(i)'last loop
        put(file," -> ");
        put(file,cnts(i)(j),1);
      end loop;
      new_line(file);
    end loop;
    put_line(file,"----+----------------+");
    put(file,"sum | "); 
    print_hms(file,totaltime);
    new_line(file);
  end Write_Filter_Counts;

  procedure Write_Factor_Counts
              ( file : in file_type;
                deco : in Standard_Natural_VecVecs.Array_of_VecVecs;
                times : in Array_of_Duration; totaltime : in duration ) is

    use Standard_Natural_VecVecs;
    use Standard_Natural_Vectors;

  begin
    new_line(file);
    new_line(file);
    put(file,"dim | ");
    put(file," CPU user time |");
    put(file," degrees of factors");
    new_line(file);
    put_line(file,"----+----------------+-------------------");
    for i in reverse 1..deco'last loop
      put(file,i,3);
      put(file," | "); print_hms(file,times(integer(i)));
      put(file," |");
      if deco(i) /= null then
        for j in deco(i)'range loop
          if(deco(i)(j) /= null)
           then put(file," "); put(file,deco(i)(j)'last,1);
          end if;
        end loop;
      end if;
      new_line(file);
    end loop;
    put_line(file,"----+----------------+");
    put(file,"sum | "); 
    print_hms(file,totaltime);
    new_line(file);
  end Write_Factor_Counts;

  procedure Write_Decomposition
              ( file : in file_type;
                deco : in Standard_Natural_VecVecs.Array_of_VecVecs ) is

    use Standard_Natural_VecVecs;
    use Standard_Natural_Vectors;

    cnt : natural32;

  begin
    new_line(file);
    for i in reverse 1..deco'last loop
      if deco(i) /= null then
        put(file,"factors of dimension "); put(file,i,1); put_line(file," :");
        cnt := 0;
        for j in deco(i)'range loop
          if(deco(i)(j) /= null) then
            cnt := cnt + 1;
            put(file,"  points in factor "); put(file,cnt,1); put(file," :");
            put(file,deco(i)(j).all); new_line(file);
          end if;
        end loop;
      end if;
    end loop;
  end Write_Decomposition;

end Path_Counts_Table;
