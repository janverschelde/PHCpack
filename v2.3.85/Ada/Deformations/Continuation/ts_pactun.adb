with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Floating_Vectors;
with Continuation_Parameters_io;
with Pack_Continuation_Parameters;

procedure ts_pactun is

-- DESCRIPTION :
--   This procedure provides an interactive tuning of the continuation
--   parameters, using the operations in Pack_Continuation_Parameters.

  procedure Main is

    v : Standard_Floating_Vectors.Vector(1..34);
    k : integer32 := 0;

  begin
    loop
      new_line;
      put_line("The current settings of the continuation parameters : ");
      Continuation_Parameters_io.put;
      put_line("The values in vector format :");
      v := Pack_Continuation_Parameters.Get;
      Pack_Continuation_Parameters.Write(v);
      put("Type in an index to change (0 to exit) : ");
      get(k);
      exit when (k = 0);
      put("Give value for v("); put(k,1); put("): "); get(v(k));
      Pack_Continuation_Parameters.Set(v);
    end loop;
  end Main;

begin
  Main;
end ts_pactun;
