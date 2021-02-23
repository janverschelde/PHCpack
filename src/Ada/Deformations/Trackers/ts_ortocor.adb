with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Complex_Solutions;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;
with Standard_System_and_Solutions_io;
with DoblDobl_Complex_Solutions;
with DoblDobl_Complex_Solutions_io;     use DoblDobl_Complex_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions_io;     use QuadDobl_Complex_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Continuation_Parameters;           use Continuation_Parameters;
with Standard_Continuation_Data;
with Standard_Orthogonal_Correctors;
with DoblDobl_Continuation_Data;
with DoblDobl_Orthogonal_Correctors;
with QuadDobl_Continuation_Data;
with QuadDobl_Orthogonal_Correctors;
with Process_io;
with Main_Poly_Continuation;            use Main_Poly_Continuation;

procedure ts_ortocor is

-- DESCRIPTION :
--   Interactive tests on the orthogonal correctors.

  procedure Run_Standard_Correctors
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in out Standard_Continuation_Data.Solu_Info_Array;
                output,conditioned : in boolean ) is

  -- DESCRIPTION :
  --   Runs the standard correctors on the solutions in s.

    use Standard_Complex_Numbers;
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;
    use Standard_Complex_Norms_Equals;
    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;
    use Standard_Orthogonal_Correctors;

    n : constant integer32 := p'last;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    cp : constant Corr_Pars := Create_End_Game;

    function Eval ( x : Vector; t : Complex_Number ) return Vector is
    begin
      return Eval(f,x);
    end Eval;

    function Diff ( x : Vector; t : Complex_Number ) return Matrix is
    begin
      return Eval(jf,x);
    end Diff;

    procedure QRLS_Silent_Correct is new
      Silent_QRLS_Corrector(Max_Norm,Eval,Diff);
    procedure QRLS_Reporting_Correct is new
      Reporting_QRLS_Corrector(Max_Norm,Eval,Diff);
    procedure SVD_Silent_Correct is new
      Silent_SVD_Corrector(Max_Norm,Eval,Diff);
    procedure SVD_Reporting_Correct is new
      Reporting_SVD_Corrector(Max_Norm,Eval,Diff);

  begin
    if output then
      Process_io.Set_Output_Code(Process_io.c);
      if conditioned then
        for i in s'range loop
          put("*** correcting solution "); put(i,1); put_line(" ***");
          SVD_Reporting_Correct(standard_output,n,s(i),cp);
        end loop;
      else
        for i in s'range loop
          put("*** correcting solution "); put(i,1); put_line(" ***");
          QRLS_Reporting_Correct(standard_output,n,s(i),cp);
        end loop;
      end if;
    else
      if conditioned then
        for i in s'range loop
          SVD_Silent_Correct(n,s(i),cp);
        end loop;
      else
        for i in s'range loop
          QRLS_Silent_Correct(n,s(i),cp);
        end loop;
      end if;
    end if;
    Clear(f); Clear(jm); Clear(jf);
  end Run_Standard_Correctors;

  procedure Run_DoblDobl_Correctors
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Continuation_Data.Solu_Info_Array;
                output,conditioned : in boolean ) is

  -- DESCRIPTION :
  --   Runs the dobldobl correctors on the solutions in s.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Vector_Norms;
    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;
    use DoblDobl_Orthogonal_Correctors;

    n : constant integer32 := p'last;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    cp : constant Corr_Pars := Create_End_Game;

    function Eval ( x : Vector; t : Complex_Number ) return Vector is
    begin
      return Eval(f,x);
    end Eval;

    function Diff ( x : Vector; t : Complex_Number ) return Matrix is
    begin
      return Eval(jf,x);
    end Diff;

    procedure QRLS_Silent_Correct is new
      Silent_QRLS_Corrector(Max_Norm,Eval,Diff);
    procedure QRLS_Reporting_Correct is new
      Reporting_QRLS_Corrector(Max_Norm,Eval,Diff);
    procedure SVD_Silent_Correct is new
      Silent_SVD_Corrector(Max_Norm,Eval,Diff);
    procedure SVD_Reporting_Correct is new
      Reporting_SVD_Corrector(Max_Norm,Eval,Diff);

  begin
    if output then
      Process_io.Set_Output_Code(Process_io.c);
      if conditioned then
        for i in s'range loop
          put("*** correcting solution "); put(i,1); put_line(" ***");
          SVD_Reporting_Correct(standard_output,n,s(i),cp);
        end loop;
      else
        for i in s'range loop
          put("*** correcting solution "); put(i,1); put_line(" ***");
          QRLS_Reporting_Correct(standard_output,n,s(i),cp);
        end loop;
      end if;
    else
      if conditioned then
        for i in s'range loop
          SVD_Silent_Correct(n,s(i),cp);
        end loop;
      else
        for i in s'range loop
          QRLS_Silent_Correct(n,s(i),cp);
        end loop;
      end if;
    end if;
    Clear(f); Clear(jm); Clear(jf);
  end Run_DoblDobl_Correctors;

  procedure Run_QuadDobl_Correctors
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Continuation_Data.Solu_Info_Array;
                output,conditioned : in boolean ) is

  -- DESCRIPTION :
  --   Runs the QuadDobl correctors on the solutions in s.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Vector_Norms;
    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;
    use QuadDobl_Orthogonal_Correctors;

    n : constant integer32 := p'last;
    f : Eval_Poly_Sys(p'range) := Create(p);
    jm : Jaco_Mat(p'range,p'range) := Create(p);
    jf : Eval_Jaco_Mat(p'range,p'range) := Create(jm);
    cp : constant Corr_Pars := Create_End_Game;

    function Eval ( x : Vector; t : Complex_Number ) return Vector is
    begin
      return Eval(f,x);
    end Eval;

    function Diff ( x : Vector; t : Complex_Number ) return Matrix is
    begin
      return Eval(jf,x);
    end Diff;

    procedure QRLS_Silent_Correct is new
      Silent_QRLS_Corrector(Max_Norm,Eval,Diff);
    procedure QRLS_Reporting_Correct is new
      Reporting_QRLS_Corrector(Max_Norm,Eval,Diff);
    procedure SVD_Silent_Correct is new
      Silent_SVD_Corrector(Max_Norm,Eval,Diff);
    procedure SVD_Reporting_Correct is new
      Reporting_SVD_Corrector(Max_Norm,Eval,Diff);

  begin
    if output then
      Process_io.Set_Output_Code(Process_io.c);
      if conditioned then
        for i in s'range loop
          put("*** correcting solution "); put(i,1); put_line(" ***");
          SVD_Reporting_Correct(standard_output,n,s(i),cp);
        end loop;
      else
        for i in s'range loop
          put("*** correcting solution "); put(i,1); put_line(" ***");
          QRLS_Reporting_Correct(standard_output,n,s(i),cp);
        end loop;
      end if;
    else
      if conditioned then
        for i in s'range loop
          SVD_Silent_Correct(n,s(i),cp);
        end loop;
      else
        for i in s'range loop
          QRLS_Silent_Correct(n,s(i),cp);
        end loop;
      end if;
    end if;
    Clear(f); Clear(jm); Clear(jf);
  end Run_QuadDobl_Correctors;

  function Menu_of_Correctors return character is

  -- DESCRIPTION :
  --   Displays the menu of the four available types of correctors.
  --   Returns the corresponding character code.

    ans : character;

  begin
    put_line("MENU for types of corrector : ");
    put_line("  1. QRLS and silent: no output, no condition number;");
    put_line("  2. SVD and silent: no output and condition number;");
    put_line("  3. QRLS and reporting: output, no condition number;");
    put_line("  4. SVD and reporting: output and condition number.");
    put("Type 1, 2, 3, or 4 to choose type : ");
    Ask_Alternative(ans,"1234");
    return ans;
  end Menu_of_Correctors;

  procedure Set_Corrector_Flags
              ( a : in character; output,condition : out boolean ) is

  -- DESCRIPTION :
  --   Depending on the type of corrector in the character a,
  --   the flags output and condition are set.

  begin
    case a is
      when '1' => output := false; condition := false;
      when '2' => output := false; condition := true;
      when '3' => output := true; condition := false;
      when '4' => output := true; condition := true;
      when others => null;
    end case;
  end Set_Corrector_Flags;

  procedure Call_Standard_Correctors
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                s : in out Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prepares the data formats to call the correctors for the system p,
  --   on a list of solutions in s and calls correctors with single doubles.

    m : constant natural32 := Standard_Complex_Solutions.Length_Of(s);
    n : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Read "); put(m,1); put_line(" solutions.");
    if m > 0 then
      put("Do you want to see the initial solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        n := natural32(Standard_Complex_Solutions.Head_Of(s).n);
        put_line("The solutions : "); put(standard_output,m,n,s);
      end if;
      ans := Menu_of_Correctors;
      declare
        use Standard_Continuation_Data;
        a : Solu_Info_Array(1..integer32(m)) := Shallow_Create(s);
        output,condition : boolean;
      begin
        Set_Corrector_Flags(ans,output,condition);   
        Driver_for_Continuation_Parameters;
        Run_Standard_Correctors(p,a,output,condition);
        s := Shallow_Create(a);
      end;
      put("Do you want to see the corrected solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        n := natural32(Standard_Complex_Solutions.Head_Of(s).n);
        put_line("The solutions : "); put(standard_output,m,n,s);
      end if;
    end if;
  end Call_Standard_Correctors;

  procedure Call_DoblDobl_Correctors
              ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prepares the data formats to call the correctors for the system p,
  --   on a list of solutions in s and calls correctors with double doubles.

    m : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(s);
    n : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Read "); put(m,1); put_line(" solutions.");
    if m > 0 then
      put("Do you want to see the initial solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        n := natural32(DoblDobl_Complex_Solutions.Head_Of(s).n);
        put_line("The solutions : "); put(standard_output,m,n,s);
      end if;
      ans := Menu_of_Correctors;
      declare
        use DoblDobl_Continuation_Data;
        a : Solu_Info_Array(1..integer32(m)) := Shallow_Create(s);
        output,condition : boolean;
      begin
        Set_Corrector_Flags(ans,output,condition);   
        Driver_for_Continuation_Parameters;
        Run_DoblDobl_Correctors(p,a,output,condition);
        s := Shallow_Create(a);
      end;
      put("Do you want to see the corrected solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        n := natural32(DoblDobl_Complex_Solutions.Head_Of(s).n);
        put_line("The solutions : "); put(standard_output,m,n,s);
      end if;
    end if;
  end Call_DoblDobl_Correctors;

  procedure Call_QuadDobl_Correctors
              ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                s : in out QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Prepares the data formats to call the correctors for the system p,
  --   on a list of solutions in s and calls correctors with quad doubles.

    m : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(s);
    n : natural32 := 0;
    ans : character;

  begin
    new_line;
    put("Read "); put(m,1); put_line(" solutions.");
    if m > 0 then
      put("Do you want to see the initial solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        n := natural32(QuadDobl_Complex_Solutions.Head_Of(s).n);
        put_line("The solutions : "); put(standard_output,m,n,s);
      end if;
      ans := Menu_of_Correctors;
      declare
        use QuadDobl_Continuation_Data;
        a : Solu_Info_Array(1..integer32(m)) := Shallow_Create(s);
        output,condition : boolean;
      begin
        Set_Corrector_Flags(ans,output,condition);   
        Driver_for_Continuation_Parameters;
        Run_QuadDobl_Correctors(p,a,output,condition);
        s := Shallow_Create(a);
      end;
      put("Do you want to see the corrected solutions ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y' then
        n := natural32(QuadDobl_Complex_Solutions.Head_Of(s).n);
        put_line("The solutions : "); put(standard_output,m,n,s);
      end if;
    end if;
  end Call_QuadDobl_Correctors;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the precision level
  --   and then calls the proper procedures.

    ans : character;
    p_sd : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    s_sd : Standard_Complex_Solutions.Solution_List;
    p_dd : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    s_dd : DoblDobl_Complex_Solutions.Solution_List;
    p_qd : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    s_qd : QuadDobl_Complex_Solutions.Solution_List;

  begin
    new_line;
    put_line("MENU to test the correctors :");
    put_line("  1. use standard double arithmetic to correct;");
    put_line("  2. use double double arithmetic to correct;");
    put_line("  3. use quad double arithmetic to correct.");
    put("Type 1, 2, or 3 to choose type of arithmetic : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' =>
        Standard_System_and_Solutions_io.get(p_sd,s_sd);
        Call_Standard_Correctors(p_sd.all,s_sd);
      when '2' =>
        DoblDobl_System_and_Solutions_io.get(p_dd,s_dd);
        Call_DoblDobl_Correctors(p_dd.all,s_dd);
      when '3' =>
        QuadDobl_System_and_Solutions_io.get(p_qd,s_qd);
        Call_QuadDobl_Correctors(p_qd.all,s_qd);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_ortocor;
