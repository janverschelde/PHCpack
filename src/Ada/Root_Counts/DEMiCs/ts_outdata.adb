with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_VecVecs_io;      use Standard_Floating_VecVecs_io;
with DEMiCs_Output_Data;

procedure ts_outdata is

-- DESCRIPTION :
--   Test on the operations in DEMiCs_Output_Data.

  procedure Test is

  -- DESCRIPTION :
  --   Test on the data stored by the package DEMiCs_Output_Data.

    lif : Standard_Floating_VecVecs.Link_to_VecVec;
    ans : character;
    idxsup,idxpnt : integer32 := 0;
    lifval,retval : double_float := 0.0;

  begin
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
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of supports
  --   and the cardinalities of each support.

    r : integer32 := 0;

  begin
    put("Give the number of distinct supports : "); get(r);
    declare
      crd : Standard_Integer_Vectors.Vector(1..r) := (1..r => 0);
    begin
      for i in 1..r loop
        put("Give the number of points in support "); put(i,1);
        put(" : "); get(crd(i));
      end loop;
      DEMiCs_Output_Data.Initialize_Lifting(crd);
    end;
    Test;
  end Main;

begin
  Main;
end ts_outdata;
