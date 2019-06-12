with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Unix_Command_Line;
with String_Splitters;                   use String_Splitters;
with Actions_and_Options;
with Option_Handlers;

procedure ts_opthand is

-- DESCRIPTION :
--   Development of the handlers of the options of an executable.

  function Unix_Command_Line_Arguments
             ( argc : integer32 ) return Array_of_Strings is

    res : Array_of_Strings(1..integer(argc));

  begin
    for i in 1..argc loop
      declare
        arg : constant string := Unix_Command_Line.Argument(natural(i));
      begin
        res(integer(i)) := new string'(arg);
      end;
    end loop;
    return res;
  end Unix_Command_Line_Arguments;

  procedure Test_Options ( args : in Array_of_Strings ) is

  -- DESCRIPTION :
  --   Prints the values of the options.

    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    vrb : constant integer32 := Actions_and_Options.Verbose_Level(args);
    seed : constant natural32 := Actions_and_Options.Find_Seed(args);
    opts : constant string := Actions_and_Options.Scan_Options(args);
    sortopts : constant string := Actions_and_Options.Sort_Options(opts);
    arg1 : constant string := Actions_and_Options.Get_Argument(args,1);
    arg2 : constant string := Actions_and_Options.Get_Argument(args,2);
    arg3 : constant string := Actions_and_Options.Get_Argument(args,3);
  
  begin
    put("The number of tasks : "); put(nt,1); new_line;
    put("The verbose level : "); put(vrb,1); new_line;
    put("The seed : "); put(seed,1); new_line;
    put("The options : "); put_line(opts);
    put("The sorted options : "); put_line(sortopts);
    put("Argument 1 : ");  put_line(arg1);
    put("Argument 2 : ");  put_line(arg2);
    put("Argument 3 : ");  put_line(arg3);
    Option_Handlers.Handle(args,sortopts,arg1,arg2,arg3);
  end Test_Options;

  procedure Main is

  -- DESCRPTION :
  --   Tests the processing of the actions and options.

    argc : constant integer32
         := integer32(Unix_Command_Line.Number_of_Arguments);

  begin
    new_line;
    put_line("Testing the processing of actions and options ...");
    put("The number of arguments : "); put(argc,1); new_line;
    if argc = 0 then
      Option_Handlers.Handle_no_Options("","");
    else -- argc > 0
      declare
        args : constant Array_of_Strings := Unix_Command_Line_Arguments(argc);
      begin
        for i in args'range loop
          put("argument "); put(integer32(i),1); put(" : ");
          put_line(args(i).all);
        end loop;
        Test_Options(args);
      end;
    end if;
  end Main;

begin
  Main;
end ts_opthand;
