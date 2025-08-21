with text_io;                           use text_io;
with Communications_with_User;
with String_Splitters;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with demics_input_data;                 use demics_input_data;

package body Test_DEMiCs_Input is

  procedure read_data_from_file is

    use demics_input_data.class_dataSet;

    filename : String_Splitters.Link_to_String;
    data : dataSet := new_dataSet;
    fail : boolean;

  begin
    put_line("Reading the name of the input file ...");
    declare
      name : constant string := String_Splitters.Read_String;
    begin
      filename := new string'(name);
    end;
    put_line(" The name of the input file : " & filename.all);
    getInputFile(data,filename,fail);
    if fail then
      put_line("Reading from file failed!");
    else
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
  end read_data_from_file;

  procedure interactive_input_data is

    use demics_input_data.class_dataSet;

    data : dataSet := new_dataSet;
    size,offset : integer32;

  begin
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
    new_line;
    put_line("The dimension, the number of distinct support sets,");
    put_line("the number of points in each support set, and");
    put_line("the number of occurrences of each support set :");
    new_line;
    info_preamble(data);
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
    new_line;
    put_line("The points in the support sets : ");
    new_line;
    info_supports(data);
  end interactive_input_data;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Testing input data for DEMiCs ...");
    new_line;
    put("Reading data from file ? (y/n) ");
    Communications_with_User.Ask_Yes_or_No(ans);
    new_line;
    if ans = 'y'
     then read_data_from_file;
     else interactive_input_data;
    end if;
  end Main;

end Test_DEMiCs_Input;
