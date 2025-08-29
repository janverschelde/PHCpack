with text_io;                           use text_io;
with String_Splitters;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;

package body DEMiCs_Input_Main is

  procedure read_data_from_file
              ( data : out demics_input_data.class_dataSet.dataSet;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use demics_input_data;
    use demics_input_data.class_dataSet;

    filename : String_Splitters.Link_to_String;

  begin
    if vrblvl > 0
     then put_line("-> in demics_input_main.Read_Data_From_File ...");
    end if;
    data := new_dataSet;
    put_line("Reading the name of the input file ...");
    declare
      name : constant string := String_Splitters.Read_String;
    begin
      filename := new string'(name);
    end;
    if vrblvl > 0
     then put_line(" The name of the input file : " & filename.all);
    end if;
    getInputFile(data,filename,fail);
    if fail then
      put_line("Reading from file failed!");
    else
      if vrblvl > 0 then
        new_line;
        put_line("The dimension, the number of distinct support sets,");
        put_line("the number of points in each support set, and");
        put_line("the number of occurrences of each support set :");
        new_line;
        info_preamble(data);
        new_line;
        put_line("The points in the support sets : ");
        new_line;
        info_supports(data);
        put_line("The name of the output file : " & data.outFile.all);
        new_line;
        put("termMax : "); put(data.termMax,1); new_line;
        put("typeMax : "); put(data.typeMax,1); new_line;
        put("termStart :");
        for i in 1..data.supN+1 loop
          put(" "); put(data.termStart(i));
        end loop;
        new_line;
      end if;
    end if;
  end read_data_from_file;

  procedure interactive_input_data
              ( data : out demics_input_data.class_dataSet.dataSet;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    use demics_input_data;
    use demics_input_data.class_dataSet;

    size,offset : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in demIcs_input_main.Interactive_Input_Data ...");
    end if;
    data := new_dataSet;
    fail := false;
    put("Give the dimension : "); get(data.dim);
    put("Give the number of distinct supports : "); get(data.supN);
    data.supType := new Standard_Integer_Vectors.Vector(1..data.supN);
    for k in 1..data.supN loop
      put("Give the number of occurrences of support ");
      put(k,1); put(" : ");
      get(data.supType(k));
    end loop;
    data.termSet := new Standard_Integer_Vectors.Vector(1..data.supN);
    for k in 1..data.supN loop
      put("Give the number of points in support ");
      put(k,1); put(" : ");
      get(data.termSet(k));
    end loop;
    if vrblvl > 0 then
      new_line;
      put_line("The dimension, the number of distinct support sets,");
      put_line("the number of points in each support set, and");
      put_line("the number of occurrences of each support set :");
      new_line;
      info_preamble(data);
    end if;
    data.termSumNum := 0;
    for k in 1..data.supN loop
      data.termSumNum := data.termSumNum + data.termSet(k);
    end loop;
    size := data.termSumNum*data.dim;
    data.support := new Standard_Floating_Vectors.Vector(1..size);
    offset := 0;
    new_line;
    for i in 1..data.supN loop
      put("Reading the points of support ");
      put(i,1); put_line(" ...");
      for j in 1..data.termSet(i) loop
        put("Give coordinates of point ");
        put(j,1); put(" : ");
        for k in 1..data.dim loop
          declare
            x : double_float := 0.0;
          begin
            get(x);
            support_in(data,offset+j,k,x);
          end;
        end loop;
      end loop;
      offset := offset + data.termSet(i);
    end loop;
    if vrblvl > 0 then
      new_line;
      put_line("The points in the support sets : ");
      new_line;
      info_supports(data);
    end if;
  exception
    when others => fail := true;
  end interactive_input_data;

end DEMiCs_Input_Main;
