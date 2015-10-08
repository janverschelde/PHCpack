with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Numbers_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors;
with Continuation_Parameters_io;
with Pack_Continuation_Parameters;

procedure ts_pactun is

-- DESCRIPTION :
--   This procedure provides an interactive tuning of the continuation
--   parameters, using the operations in Pack_Continuation_Parameters.

  procedure Main is

  -- DESCRIPTION :
  --   Displays the values of the continuation parameters
  --   and the user can adjust the values via the packed representation
  --   as a vector of 34 floats.

    v : Standard_Floating_Vectors.Vector(1..34);
    k : integer32 := 0;
    val : double_float;
    ans : character;

  begin
    loop
      new_line;
      put_line("The current settings of the continuation parameters : ");
      Continuation_Parameters_io.put;
      put_line("The values in vector format :");
      v := Pack_Continuation_Parameters.Get;
      Pack_Continuation_Parameters.Write(v);
      put("Type in an index to change (0 to exit) : ");
      Numbers_io.Read_Integer(k);
      exit when (k = 0);
      put("Give value for v("); put(k,1); put(") : ");
      Numbers_io.Read_Double_Float(v(k));
      Pack_Continuation_Parameters.Set(v);
      loop
        put("-> the value for index "); put(k,1); put(" : ");
        put(Pack_Continuation_Parameters.Get_Value(natural32(k)),3);
        new_line;
        put("-> change this value ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
        put("-> give a value for index "); put(k,1); put(" : ");
        Numbers_io.Read_Double_Float(val);
        Pack_Continuation_Parameters.Set_Value(natural32(k),val);
      end loop;
    end loop;
  end Main;

begin
  Main;
end ts_pactun;
