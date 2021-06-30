with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Unix_Command_Line;
with String_Splitters;                   use String_Splitters;
with Actions_and_Options;
with Option_Handlers;
with use_c2phc;    -- to force compilation on Windows
with pieri_solver; -- to force compilation

procedure phctop is

  function Unix_Command_Line_Arguments
             ( argc : integer32 ) return Array_of_Strings is

  -- DESCRIPTION :
  --   Returns the command line arguments.

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

  procedure Handle_Options ( args : in Array_of_Strings ) is

  -- DESCRIPTION :
  --   Handles the command line arguments, sorting the options
  --   so the actions come before the pure options.
  --   If the verbose level is positive, then extra information is written.

    nt : constant natural32 := Actions_and_Options.Number_of_Tasks(args);
    vrb : constant integer32 := Actions_and_Options.Verbose_Level(args);
    seed : constant natural32 := Actions_and_Options.Find_Seed(args);
    opts : constant string := Actions_and_Options.Scan_Options(args);
    sortopts : constant string := Actions_and_Options.Sort_Options(opts);
    arg1 : constant string := Actions_and_Options.Get_Argument(args,1);
    arg2 : constant string := Actions_and_Options.Get_Argument(args,2);
    arg3 : constant string := Actions_and_Options.Get_Argument(args,3);
  
  begin
    if vrb > 0 then
      put("At verbose level "); put(vrb,1);
      put_line(", in phctop.Handle_Options ...");
      for k in args'range loop
        put("argument "); put(natural32(k),1); put(" : ");
        put_line(args(k).all);
      end loop;
      put("The number of tasks : "); put(nt,1); new_line;
      put("The verbose level : "); put(vrb,1); new_line;
      put("The seed : "); put(seed,1); new_line;
      put("The options : "); put_line(opts);
      put("The sorted options : "); put_line(sortopts);
      put("Argument 1 : ");  put_line(arg1);
      put("Argument 2 : ");  put_line(arg2);
      put("Argument 3 : ");  put_line(arg3);
    end if;
    Option_Handlers.Handle(args,sortopts,arg1,arg2,arg3,vrb-1);
  end Handle_Options;

  procedure Main is

  -- DESCRPTION :
  --   Gets the command line arguments and passes those to the handlers.

    argc : constant integer32
         := integer32(Unix_Command_Line.Number_of_Arguments);

  begin
    if argc = 0 then
      Option_Handlers.Handle_no_Options("","");
    else -- argc > 0
      declare
        args : constant Array_of_Strings := Unix_Command_Line_Arguments(argc);
      begin
        Handle_Options(args);
      end;
    end if;
  end Main;

begin
  Main;
end phctop;
