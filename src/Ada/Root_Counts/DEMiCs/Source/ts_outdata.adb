with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with String_Splitters;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;      use Standard_Floating_VecVecs_io;
with Lists_of_Strings;
with DEMiCs_Output_Data;

procedure ts_outdata is

-- DESCRIPTION :
--   Test on the operations in DEMiCs_Output_Data.

  procedure Test_on_Lifting_Values ( r : in integer32 ) is

  -- DESCRIPTION :
  --   Test on the data stored by the package DEMiCs_Output_Data.
  --   The number r on input is the number of distinct supports.
  --   Prompts the user for the cardinalities of each support
  --   and then asks for lifting values for some points.

    crd : Standard_Integer_Vectors.Vector(1..r) := (1..r => 0);
    lif : Standard_Floating_VecVecs.Link_to_VecVec;
    ans : character;
    idxsup,idxpnt : integer32 := 0;
    lifval,retval : double_float := 0.0;

  begin
    for i in 1..r loop
      put("Give the number of points in support "); put(i,1);
      put(" : "); get(crd(i));
    end loop;
    DEMiCs_Output_Data.Initialize_Lifting(crd);
    loop
      lif := DEMiCs_Output_Data.Lifting_Values;
      put_line("The lifting values : "); put(lif.all);
      put("Continue ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("Give an index for a support : "); get(idxsup);
      put("Give an index for a point : "); get(idxpnt);
      put("Give a lifting value : "); get(lifval);
      DEMiCs_Output_Data.Assign_Lifting(idxsup,idxpnt,lifval);
      retval := DEMiCs_Output_Data.Retrieve_Lifting(idxsup,idxpnt);
      put("Retrieved lifting value : "); put(retval); new_line;
    end loop;
  end Test_on_Lifting_Values;

  procedure Test_Lifting is

  -- DESCRIPTION :
  --   Prompts the user for the number of supports
  --   and then tests the storing of lifting values.

    r : integer32 := 0;

  begin
    put("Give the number of distinct supports : "); get(r);
    Test_on_Lifting_Values(r);
  end Test_Lifting;

  procedure Write_Strings is

  -- DESCRIPTION :
  --   Writes the strings stored in DEMiCs_Output_Data.

    data : constant Lists_of_Strings.List
         := DEMiCs_Output_Data.Retrieve_Cell_Indices;
    tmp : Lists_of_Strings.List := data;
    sls : String_Splitters.Link_to_String;
    cnt : integer32 := 0;

  begin
    while not Lists_of_Strings.Is_Null(tmp) loop
      sls := Lists_of_Strings.Head_Of(tmp);
      cnt := cnt + 1;
      put(cnt,1); put(" : "); put_line(sls.all);
      tmp := Lists_of_Strings.Tail_Of(tmp);
    end loop;
  end Write_Strings;

  procedure Prompt_for_Index is

  -- DESCRIPTION :
  --   Prompts for an index and then shows the string
  --   at the location defined by the index.

    index : integer32 := 0;
    ls : String_Splitters.Link_to_String;

    use String_Splitters;

  begin
    loop
      put("Give an index to a string (0 to exit) : ");
      get(index);
      exit when (index = 0);
      ls := DEMiCs_Output_Data.Get_Cell_Indices(index);
      if ls /= null then
        put("The string at position "); put(index);
        put_line(" :"); put_line(ls.all);
      end if;
    end loop;
  end Prompt_for_Index;

  procedure Test_on_Strings is

  -- DESCRIPTION :
  --   Prompts the user for strings and stores the strings
  --   then in DEMiCs_Output_Data.

    ans : character;

  begin
    loop
      put_line("Reading a string ...");
      declare
        s : constant String := String_Splitters.Read_String;
      begin
        DEMiCs_Output_Data.Add_Cell_Indices(s);
      end;
      Write_Strings;
      Prompt_for_Index;
      put("Continue ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_on_Strings;

  procedure Main is

  -- DESCRIPTION :
  --   Tests the operations in the package DEMiCs_Output_Data.

    ans : character;

  begin
    new_line;
    put_line("MENU to test DEMiCs_Output_Data :");
    put_line("  1. test on adding lifting values;");
    put_line("  2. test on storing strings.");
    put("Type 1 or 2 for a test : ");
    Ask_Alternative(ans,"12");
    new_line;
    case ans is
      when '1' => Test_Lifting;
      when '2' => Test_on_Strings;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_outdata;
