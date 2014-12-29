with text_io;                            use text_io;
with String_Splitters;                   use String_Splitters;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Symbol_Table,Symbol_Table_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;   use Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;   use DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;   use QuadDobl_System_and_Solutions_io;
with Multprec_Floating_Numbers;
with Multprec_Complex_Solutions;
with Multprec_Complex_Solutions_io;
with Multprec_Solution_Strings;
with Multprec_System_and_Solutions_io;   use Multprec_System_and_Solutions_io;
with Standard_to_Multprec_Convertors;
with Standard_Path_Tracker;
with DoblDobl_Path_Tracker;
with QuadDobl_Path_Tracker;
with Multprec_Path_Tracker;
with Verification_of_Solutions;
with Varbprec_Path_Tracker;

-- for testing purposes
-- with Process_io;
with Standard_Homotopy;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;
 use Standard_Complex_Vectors_io;
with DoblDobl_Homotopy;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;
 use DoblDobl_Complex_Vectors_io;
with QuadDobl_Homotopy;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;
 use QuadDobl_Complex_Vectors_io;
with Multprec_Homotopy;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;
 use Multprec_Complex_Vectors_io;

procedure ts_nxtsol is

-- DESCRIPTION :
--   Allows tracking of one path via a generator "next"
--   which gives the next solution point on the path.

  procedure Standard_Natural_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy and start solutions to
  --   initialize the standard path tracker.

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;

    hom : Link_to_Poly_Sys;
    idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the homotopy ...");
    get(hom);
    new_line;
    put("The symbols in the homotopy : "); Symbol_Table_io.Write; new_line;
    put("-> give the index of the deformation parameter t : ");
    get(idx); skip_line;
    if idx /= hom'last+1
     then put_line("Works only when deformation parameter is last symbol!");
    end if;
    Symbol_Table.Remove(natural32(idx)); -- remove t for reading of solutions
    new_line;
    Standard_Complex_Solutions_io.Read(sols);
    Standard_Path_Tracker.Init(hom,idx,Head_Of(sols));
   -- only good when standard_path_tracker has reporting corrector
   -- declare
   --   use Process_io;
   --   oc : constant Process_io.Output_Code := c;
   -- begin
   --   Process_io.Set_Output_Code(oc);
   -- end;
    declare
      y : constant Standard_Complex_Vectors.Vector
        := Standard_Homotopy.Eval(Head_Of(sols).v,Head_Of(sols).t);
    begin
      put_line("The start solution evaluated : ");
      put_line(y);
    end;
  end Standard_Natural_Parameter_Initialize;

  procedure DoblDobl_Natural_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy and start solutions to
  --   initialize the double double path tracker.

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;

    hom : Link_to_Poly_Sys;
    idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the homotopy ...");
    get(hom);
    new_line;
    put("The symbols in the homotopy : "); Symbol_Table_io.Write; new_line;
    put("-> give the index of the deformation parameter t : ");
    get(idx); skip_line;
    if idx /= hom'last+1
     then put_line("Works only when deformation parameter is last symbol!");
    end if;
    Symbol_Table.Remove(natural32(idx)); -- remove t for reading of solutions
    new_line;
    DoblDobl_Complex_Solutions_io.Read(sols);
    DoblDobl_Path_Tracker.Init(hom,idx,Head_Of(sols));
   -- only good when standard_path_tracker has reporting corrector
   -- declare
   --   use Process_io;
   --   oc : constant Process_io.Output_Code := c;
   -- begin
   --   Process_io.Set_Output_Code(oc);
   -- end;
    declare
      y : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Homotopy.Eval(Head_Of(sols).v,Head_Of(sols).t);
    begin
      put_line("The start solution evaluated : ");
      put_line(y);
    end;
  end DoblDobl_Natural_Parameter_Initialize;

  procedure QuadDobl_Natural_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy and start solutions to
  --   initialize the quad double path tracker.

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;

    hom : Link_to_Poly_Sys;
    idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the homotopy ...");
    get(hom);
    new_line;
    put("The symbols in the homotopy : "); Symbol_Table_io.Write; new_line;
    put("-> give the index of the deformation parameter t : ");
    get(idx); skip_line;
    if idx /= hom'last+1
     then put_line("Works only when deformation parameter is last symbol!");
    end if;
    Symbol_Table.Remove(natural32(idx)); -- remove t for reading of solutions
    new_line;
    QuadDobl_Complex_Solutions_io.Read(sols);
    QuadDobl_Path_Tracker.Init(hom,idx,Head_Of(sols));
   -- only good when standard_path_tracker has reporting corrector
   -- declare
   --   use Process_io;
   --   oc : constant Process_io.Output_Code := c;
   -- begin
   --   Process_io.Set_Output_Code(oc);
   -- end;
    declare
      y : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Homotopy.Eval(Head_Of(sols).v,Head_Of(sols).t);
    begin
      put_line("The start solution evaluated : ");
      put_line(y);
    end;
  end QuadDobl_Natural_Parameter_Initialize;

  procedure Multprec_Natural_Parameter_Initialize
              ( deci : in natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy and start solutions to
  --   initialize the multiprecision path tracker.
  --   The number of decimal places in the working precision equals deci.

    use Multprec_Complex_Poly_Systems,Multprec_Complex_Solutions;

    hom : Link_to_Poly_Sys;
    idx : integer32 := 0;
    sols : Solution_List;

  begin
    new_line;
    put_line("Reading the homotopy ...");
    get(hom);
    new_line;
    put("The symbols in the homotopy : "); Symbol_Table_io.Write; new_line;
    put("-> give the index of the deformation parameter t : ");
    get(idx); skip_line;
    if idx /= hom'last+1
     then put_line("Works only when deformation parameter is last symbol!");
    end if;
    Symbol_Table.Remove(natural32(idx)); -- remove t for reading of solutions
    new_line;
    Multprec_Complex_Solutions_io.Read(sols);
    Multprec_Path_Tracker.Init(hom,idx,0,deci,Head_Of(sols));
   -- only good when standard_path_tracker has reporting corrector
   -- declare
   --   use Process_io;
   --   oc : constant Process_io.Output_Code := c;
   -- begin
   --   Process_io.Set_Output_Code(oc);
   -- end;
    declare
      y : constant Multprec_Complex_Vectors.Vector
        := Multprec_Homotopy.Eval(Head_Of(sols).v,Head_Of(sols).t);
    begin
      put_line("The start solution evaluated : ");
      put_line(y);
    end;
  end Multprec_Natural_Parameter_Initialize;

  procedure Varbprec_Natural_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a homotopy and start solutions to
  --   initialize the variable precision path tracker.

    hom : Link_to_Array_of_Strings;
    sol : Link_to_String;
    sols : Multprec_Complex_Solutions.Solution_List;
    ls : Multprec_Complex_Solutions.Link_to_Solution;
    nq,nv,len : natural32;
    idx : integer32 := 0;

  begin
    new_line;
    put_line("Reading the homotopy and start solutions from file ...");
    Multprec_System_and_Solutions_io.get(nq,nv,hom,sols);
    new_line;
    put("-> read "); put(nq,1); put(" polynomials in "); put(nv,1);
    put_line(" variables.");
    len := Multprec_Complex_Solutions.Length_Of(sols);
    put("-> read "); put(len,1); put_line(" solutions from file.");
    if len > 0 then
      ls := Multprec_Complex_Solutions.Head_Of(sols);
      sol := new string'(Multprec_Solution_Strings.Write(ls.all));
      put_line("-> the first solution :"); put_line(sol.all);

      put("The symbols in the homotopy : "); Symbol_Table_io.Write; new_line;
      put("-> give the index of the deformation parameter t : ");
      get(idx); skip_line;
      if idx /= integer32(hom'last)+1
       then put_line("Works only when deformation parameter is last symbol!");
      end if;
      Varbprec_Path_Tracker.Init(hom,idx,sol);
    end if;
  end Varbprec_Natural_Parameter_Initialize;

  procedure Standard_Artificial_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the standard path tracker.

    use Standard_Complex_Poly_Systems,Standard_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' 
     then Standard_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols));
     else Standard_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols));
    end if;
  end Standard_Artificial_Parameter_Initialize;

  procedure DoblDobl_Artificial_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the double double path tracker.

    use DoblDobl_Complex_Poly_Systems,DoblDobl_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols));
     else DoblDobl_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols));
    end if;
  end DoblDobl_Artificial_Parameter_Initialize;

  procedure QuadDobl_Artificial_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the double double path tracker.

    use QuadDobl_Complex_Poly_Systems,QuadDobl_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols));
     else QuadDobl_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols));
    end if;
  end QuadDobl_Artificial_Parameter_Initialize;

  procedure Multprec_Artificial_Parameter_Initialize
              ( deci : in natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system,
  --   and start solutions and initializes the multiprecision path tracker.

    use Multprec_Complex_Poly_Systems,Multprec_Complex_Solutions;

    tgt_sys,sta_sys : Link_to_Poly_Sys;
    sols : Solution_List;
    ans : character;
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(deci);

  begin
    new_line;
    put_line("Reading the target system ...");
    get(tgt_sys);
    Standard_to_Multprec_Convertors.Set_Size(tgt_sys.all,size);
    new_line;
    put_line("Reading the start system and its solutions...");
    get(sta_sys,sols);
    Standard_to_Multprec_Convertors.Set_Size(sta_sys.all,size);
    Multprec_Complex_Solutions.Set_Size(sols,size);
    new_line;
    put("Fixed gamma constant ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' 
     then Multprec_Path_Tracker.Init(tgt_sys,sta_sys,true,Head_Of(sols),deci);
     else Multprec_Path_Tracker.Init(tgt_sys,sta_sys,false,Head_Of(sols),deci);
    end if;
  end Multprec_Artificial_Parameter_Initialize;

  procedure Varbprec_Artificial_Parameter_Initialize is

  -- DESCRIPTION :
  --   Prompts the user for a target system, a start system, start solutions,
  --   and initializes the variable precision path tracker.

    pfile : file_type;
    p,q : Link_to_Array_of_Strings;
    sol : Link_to_String;
    sols : Multprec_Complex_Solutions.Solution_List;
    ls : Multprec_Complex_Solutions.Link_to_Solution;
    nq,nv,len : natural32;
    ans : character;

  begin
    new_line;
    put_line("Reading the target system ...");
    Read_Name_and_Open_File(pfile);
    get(pfile,natural(nq),natural(nv),p);
    new_line;
    put("-> read "); put(nq,1); put(" polynomials in "); put(nv,1);
    put_line(" variables.");
    new_line;
    put_line("Reading the start system and start solutions from file ...");
    Multprec_System_and_Solutions_io.get(nq,nv,q,sols);
    new_line;
    put("-> read "); put(nq,1); put(" polynomials in "); put(nv,1);
    put_line(" variables.");
    len := Multprec_Complex_Solutions.Length_Of(sols);
    put("-> read "); put(len,1); put_line(" solutions from file.");
    if len > 0 then
      ls := Multprec_Complex_Solutions.Head_Of(sols);
      sol := new string'(Multprec_Solution_Strings.Write(ls.all));
      put_line("-> the first solution :"); put_line(sol.all);
      new_line;
      put("Fixed gamma constant ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then Varbprec_Path_Tracker.Init(p,q,true,sol);
       else Varbprec_Path_Tracker.Init(p,q,false,sol);
      end if;
    end if;
  end Varbprec_Artificial_Parameter_Initialize;

  procedure Standard_Initialize_Path_Tracker is

  -- DESCRIPTION :
  --   Prompts the user for the choice between natural parameter,
  --   or artificial parameter homotopy, with start and target system.

    ans : character;

  begin
    new_line;
    put_line("MENU to choose between artificial or natural parameter :");
    put_line("  1. artificial parameter t with start and target system;");
    put_line("  2. natural parameter t is a variable in the homotopy system.");
    put("Type 1 or 2 to select type of homotopy : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Standard_Artificial_Parameter_Initialize;
     else Standard_Natural_Parameter_Initialize;
    end if;
  end Standard_Initialize_Path_Tracker;

  procedure DoblDobl_Initialize_Path_Tracker is

  -- DESCRIPTION :
  --   Prompts the user for the choice between natural parameter,
  --   or artificial parameter homotopy, with start and target system.

    ans : character;

  begin
    new_line;
    put_line("MENU to choose between artificial or natural parameter :");
    put_line("  1. artificial parameter t with start and target system;");
    put_line("  2. natural parameter t is a variable in the homotopy system.");
    put("Type 1 or 2 to select type of homotopy : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then DoblDobl_Artificial_Parameter_Initialize;
     else DoblDobl_Natural_Parameter_Initialize;
    end if;
  end DoblDobl_Initialize_Path_Tracker;

  procedure QuadDobl_Initialize_Path_Tracker is

  -- DESCRIPTION :
  --   Prompts the user for the choice between natural parameter,
  --   or artificial parameter homotopy, with start and target system.

    ans : character;

  begin
    new_line;
    put_line("MENU to choose between artificial or natural parameter :");
    put_line("  1. artificial parameter t with start and target system;");
    put_line("  2. natural parameter t is a variable in the homotopy system.");
    put("Type 1 or 2 to select type of homotopy : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then QuadDobl_Artificial_Parameter_Initialize;
     else QuadDobl_Natural_Parameter_Initialize;
    end if;
  end QuadDobl_Initialize_Path_Tracker;

  procedure Multprec_Initialize_Path_Tracker ( deci : in natural32 ) is

  -- DESCRIPTION :
  --   Prompts the user for the choice between natural parameter,
  --   or artificial parameter homotopy, with start and target system.

    ans : character;

  begin
    new_line;
    put_line("MENU to choose between artificial or natural parameter :");
    put_line("  1. artificial parameter t with start and target system;");
    put_line("  2. natural parameter t is a variable in the homotopy system.");
    put("Type 1 or 2 to select type of homotopy : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Multprec_Artificial_Parameter_Initialize(deci);
     else Multprec_Natural_Parameter_Initialize(deci);
    end if;
  end Multprec_Initialize_Path_Tracker;

  procedure Varbprec_Initialize_Path_Tracker is

  -- DESCRIPTION :
  --   Prompts the user for the choice between natural parameter,
  --   or artificial parameter homotopy, with start and target system.

    ans : character;

  begin
    new_line;
    put_line("MENU to choose between artificial or natural parameter :");
    put_line("  1. artificial parameter t with start and target system;");
    put_line("  2. natural parameter t is a variable in the homotopy system.");
    put("Type 1 or 2 to select type of homotopy : ");
    Ask_Alternative(ans,"12");
    if ans = '1'
     then Varbprec_Artificial_Parameter_Initialize;
     else Varbprec_Natural_Parameter_Initialize;
    end if;
  end Varbprec_Initialize_Path_Tracker;

  procedure Standard_Run_Path_Tracker is

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with standard double precision arithmetic.

    s : Standard_Complex_Solutions.Link_to_Solution
      := Standard_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    Standard_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := Standard_Path_Tracker.get_next;
      put_line("The current solution : ");
      Standard_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end Standard_Run_Path_Tracker;

  procedure DoblDobl_Run_Path_Tracker is

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with double double precision arithmetic.
    
    s : DoblDobl_Complex_Solutions.Link_to_Solution
      := DoblDobl_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    DoblDobl_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := DoblDobl_Path_Tracker.get_next;
      put_line("The current solution : ");
      DoblDobl_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end DoblDobl_Run_Path_Tracker;

  procedure QuadDobl_Run_Path_Tracker is

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with quad double precision arithmetic.
    
    s : QuadDobl_Complex_Solutions.Link_to_Solution
      := QuadDobl_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    QuadDobl_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := QuadDobl_Path_Tracker.get_next;
      put_line("The current solution : ");
      QuadDobl_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end QuadDobl_Run_Path_Tracker;

  procedure Multprec_Run_Path_Tracker is

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with multiprecision arithmetic.
    
    s : Multprec_Complex_Solutions.Link_to_Solution
      := Multprec_Path_Tracker.get_current;
    ans : character;

  begin
    new_line;
    put_line("The current solution : ");
    Multprec_Complex_Solutions_io.put(s.all); new_line;
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := Multprec_Path_Tracker.get_next;
      put_line("The current solution : ");
      Multprec_Complex_Solutions_io.put(s.all); new_line;
    end loop;
  end Multprec_Run_Path_Tracker;

  procedure Varbprec_Run_Path_Tracker is

  -- DESCRIPTION :
  --   Shows the current solution and prompts the user each time to make
  --   one extra predictor-corrector step to compute the next solution
  --   on the path with variable precision arithmetic.

    use Verification_of_Solutions;

    s : Link_to_String := VarbPrec_Path_Tracker.get_current;
    ans : character;
    wanted,maxitr,maxprc : natural32;
    verbose : boolean;

  begin
    Menu_to_Set_Parameters(wanted,maxitr,maxprc,verbose);
    new_line;
    put_line("The current solution : "); put_line(s.all);
    loop
      put("Do predictor-corrector step ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      s := Varbprec_Path_Tracker.get_next(wanted,maxitr,maxprc,verbose);
      put_line("The current solution : "); put_line(s.all);
    end loop;
  end Varbprec_Run_Path_Tracker;

  procedure Main is

    ans : character;
    deci : natural32 := 0;

  begin
    new_line;
    put_line("Tracking a path with a generator ...");
    put_line("  1. in complex standard floating-point arithmetic;");
    put_line("  2. in complex double double arithmetic;");
    put_line("  3. in complex quad double arithmetic;");
    put_line("  4. in complex multiprecision arithmetic;");
    put_line("  5. in variable multiprecision arithmetic.");
    put("Type 1, 2, 3, 4, or 5 to select precision : ");
    Ask_Alternative(ans,"12345");
    case ans is
      when '1' =>
        Standard_Initialize_Path_Tracker;
        Standard_Run_Path_Tracker;
        Standard_Path_Tracker.Clear;
      when '2' =>
        DoblDobl_Initialize_Path_Tracker;
        DoblDobl_Run_Path_Tracker;
        DoblDobl_Path_Tracker.Clear;
      when '3' =>
        QuadDobl_Initialize_Path_Tracker;
        QuadDobl_Run_Path_Tracker;
        QuadDobl_Path_Tracker.Clear;
      when '4' =>
        new_line;
        put("Give the number of decimal places in the working precision : ");
        get(deci);
        Multprec_Initialize_Path_Tracker(deci);
        Multprec_Run_Path_Tracker;
        Multprec_Path_Tracker.Clear;
      when '5' =>
        Varbprec_Initialize_Path_Tracker;
        Varbprec_Run_Path_Tracker;
        Varbprec_Path_Tracker.Clear;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_nxtsol;
