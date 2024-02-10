with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Floating_Vectors;
with Standard_Homotopy;
with Standard_Laurent_Homotopy;
with DoblDobl_Homotopy;
with DoblDobl_Laurent_Homotopy;
with QuadDobl_Homotopy;
with QuadDobl_Laurent_Homotopy;
with Multprec_Homotopy;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Continuation_Data_io;
with DoblDobl_Continuation_Data_io;
with QuadDobl_Continuation_Data_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Assignments_of_Solutions;          use Assignments_of_Solutions;
with PHCpack_Operations;
with Crude_Path_Trackers;

package body Path_Trackers_Interface is

  function Path_Trackers_Standard_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Homotopy_Random ...");
    end if;
    PHCpack_Operations.Create_Standard_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Homotopy_Random.");
      end if;
     return 152;
  end Path_Trackers_Standard_Homotopy_Random;

  function Path_Trackers_Standard_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    g : Standard_Floating_Vectors.Vector(1..2);
    gamma : Standard_Complex_Numbers.Complex_Number;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    pwt : constant natural32 := natural32(v_a(v_a'first)); -- power of t

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Homotopy_Gamma ...");
    end if;
    Assign(2,c,g);
    gamma := Standard_Complex_Numbers.Create(g(1),g(2));
    if pwt > 0 then -- the power of t must be positive
      PHCpack_Operations.Create_Standard_Homotopy(gamma,pwt);
    else
      if vrblvl > 0 then
        put("Invalid input "); put(pwt,1);
        put_line(" for the power of t.");
        put_line("Will resort to the default power of t ...");
      end if;
      PHCpack_Operations.Create_Standard_Homotopy(gamma);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Homotopy_Gamma.");
      end if;
     return 153;
  end Path_Trackers_Standard_Homotopy_Gamma;

  function Path_Trackers_Standard_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Homotopy_Clear ...");
    end if;
    Standard_Homotopy.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Homotopy_Clear.");
      end if;
     return 154;
  end Path_Trackers_Standard_Homotopy_Clear;

  function Path_Trackers_DoblDobl_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Homotopy ...");
    end if;
    PHCpack_Operations.Create_DoblDobl_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Homotopy_Random.");
      end if;
     return 172;
  end Path_Trackers_DoblDobl_Homotopy_Random;

  function Path_Trackers_DoblDobl_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : double_double;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    pwt : constant natural32 := natural32(v_a(v_a'first)); -- power of t

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Homotopy_Gamma ...");
    end if;
    Assign(2,c,g);
    g_re := create(g(1));
    g_im := create(g(2));
    gamma := DoblDobl_Complex_Numbers.Create(g_re,g_im);
    if pwt > 0 then -- the power of t must be positive
      PHCpack_Operations.Create_DoblDobl_Homotopy(gamma,pwt);
    else
      if vrblvl > 0 then
        put("Invalid input "); put(pwt,1);
        put_line(" for the power of t.");
        put_line("Will resort to the default power of t ...");
      end if;
      PHCpack_Operations.Create_DoblDobl_Homotopy(gamma);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Homotopy_Gamma.");
      end if;
     return 173;
  end Path_Trackers_DoblDobl_Homotopy_Gamma;

  function Path_Trackers_DoblDobl_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Homotopy_Clear ...");
    end if;
    DoblDobl_Homotopy.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Homotopy_Clear.");
      end if;
     return 174;
  end Path_Trackers_DoblDobl_Homotopy_Clear;

  function Path_Trackers_QuadDobl_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Homotopy_Random ...");
    end if;
    PHCpack_Operations.Create_QuadDobl_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Homotopy_Random.");
      end if;
     return 182;
  end Path_Trackers_QuadDobl_Homotopy_Random;

  function Path_Trackers_QuadDobl_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : quad_double;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    pwt : constant natural32 := natural32(v_a(v_a'first)); -- power of t

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Homotopy_Gamma ...");
    end if;
    Assign(2,c,g);
    g_re := create(g(1));
    g_im := create(g(1));
    gamma := QuadDobl_Complex_Numbers.Create(g_re,g_im);
    if pwt > 0 then -- the power of t must be positive
      PHCpack_Operations.Create_QuadDobl_Homotopy(gamma,pwt);
    else
      if vrblvl > 0 then
        put("Invalid input "); put(pwt,1);
        put_line(" for the power of t.");
        put_line("Will resort to the default power of t ...");
      end if;
      PHCpack_Operations.Create_QuadDobl_Homotopy(gamma);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Homotopy_Gamma.");
      end if;
     return 183;
  end Path_Trackers_QuadDobl_Homotopy_Gamma;

  function Path_Trackers_QuadDobl_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Homotopy_Clear ...");
    end if;
    QuadDobl_Homotopy.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Homotopy_Clear.");
      end if;
     return 184;
  end Path_Trackers_QuadDobl_Homotopy_Clear;

  function Path_Trackers_Multprec_Homotopy_Random
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Multprec_Homotopy_Random ...");
    end if;
    PHCpack_Operations.Create_Multprec_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Multprec_Homotopy_Random.");
      end if;
     return 522;
  end Path_Trackers_Multprec_Homotopy_Random;

  function Path_Trackers_Multprec_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : Multprec_Floating_Numbers.Floating_Number;
    gamma : Multprec_Complex_Numbers.Complex_Number;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    pwt : constant natural32 := natural32(v_a(v_a'first)); -- power of t

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Multprec_Homotopy_Gamma ...");
    end if;
    Assign(2,c,g);
    g_re := Multprec_Floating_Numbers.create(g(1));
    g_im := Multprec_Floating_Numbers.create(g(1));
    gamma := Multprec_Complex_Numbers.Create(g_re,g_im);
    Multprec_Floating_Numbers.Clear(g_re);
    Multprec_Floating_Numbers.Clear(g_im);
    if pwt > 0 then -- the power of t must be positive
      PHCpack_Operations.Create_Multprec_Homotopy(gamma,pwt);
    else
      if vrblvl > 0 then
        put("Invalid input "); put(pwt,1);
        put_line(" for the power of t.");
        put_line("Will resort to the default power of t ...");
      end if;
      PHCpack_Operations.Create_Multprec_Homotopy(gamma);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Multprec_Homotopy_Gamma.");
      end if;
     return 523;
  end Path_Trackers_Multprec_Homotopy_Gamma;

  function Path_Trackers_Multprec_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Multprec_Homotopy_Clear ...");
    end if;
    Multprec_Homotopy.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Multprec_Homotopy_Clear.");
      end if;
     return 524;
  end Path_Trackers_Multprec_Homotopy_Clear;

  function Path_Trackers_Standard_Laurent_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Laurent_Homotopy ...");
    end if;
    PHCpack_Operations.Create_Standard_Laurent_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Laurent_Homotopy.");
      end if;
     return 792;
  end Path_Trackers_Standard_Laurent_Homotopy;

  function Path_Trackers_Standard_Laurent_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    g : Standard_Floating_Vectors.Vector(1..2);
    gamma : Standard_Complex_Numbers.Complex_Number;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    pwt : constant natural32 := natural32(v_a(v_a'first)); -- power of t

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Laurent_Homotopy_Gamma ...");
    end if;
    Assign(2,c,g);
    gamma := Standard_Complex_Numbers.Create(g(1),g(2));
    if pwt > 0 then -- the power of t must be positive
      PHCpack_Operations.Create_Standard_Laurent_Homotopy(gamma,pwt);
    else
      if vrblvl > 0 then
        put("Invalid input "); put(pwt,1);
        put_line(" for the power of t.");
        put_line("Will resort to the default power of t ...");
      end if;
      PHCpack_Operations.Create_Standard_Laurent_Homotopy(gamma);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Laurent_Homotopy_Gamma.");
      end if;
     return 921;
  end Path_Trackers_Standard_Laurent_Homotopy_Gamma;

  function Path_Trackers_Standard_Laurent_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Laurent_Homotopy_Clear ...");
    end if;
    Standard_Laurent_Homotopy.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Laurent_Homotopy_Clear.");
      end if;
     return 924;
  end Path_Trackers_Standard_Laurent_Homotopy_Clear;

  function Path_Trackers_DoblDobl_Laurent_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Laurent_Homotopy ...");
    end if;
    PHCpack_Operations.Create_DoblDobl_Laurent_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Laurent_Homotopy.");
      end if;
     return 793;
  end Path_Trackers_DoblDobl_Laurent_Homotopy;

  function Path_Trackers_DoblDobl_Laurent_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : double_double;
    gamma : DoblDobl_Complex_Numbers.Complex_Number;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    pwt : constant natural32 := natural32(v_a(v_a'first)); -- power of t

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Laurent_Homotopy_Gamma ...");
    end if;
    Assign(2,c,g);
    g_re := create(g(1));
    g_im := create(g(2));
    gamma := DoblDobl_Complex_Numbers.Create(g_re,g_im);
    if pwt > 0 then -- the power of t must be positive
      PHCpack_Operations.Create_DoblDobl_Laurent_Homotopy(gamma,pwt);
    else
      if vrblvl > 0 then
        put("Invalid input "); put(pwt,1);
        put_line(" for the power of t.");
        put_line("Will resort to the default power of t ...");
      end if;
      PHCpack_Operations.Create_DoblDobl_Laurent_Homotopy(gamma);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Laurent_Homotopy_Gamma.");
      end if;
     return 922;
  end Path_Trackers_DoblDobl_Laurent_Homotopy_Gamma;

  function Path_Trackers_DoblDobl_Laurent_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Laurent_Homotopy_Clear ...");
    end if;
    DoblDobl_Laurent_Homotopy.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Laurent_Homotopy_Clear.");
      end if;
     return 925;
  end Path_Trackers_DoblDobl_Laurent_Homotopy_Clear;

  function Path_Trackers_QuadDobl_Laurent_Homotopy
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Laurent_Homotopy ...");
    end if;
    PHCpack_Operations.Create_QuadDobl_Laurent_Homotopy;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Laurent_Homotopy.");
      end if;
     return 794;
  end Path_Trackers_QuadDobl_Laurent_Homotopy;

  function Path_Trackers_QuadDobl_Laurent_Homotopy_Gamma
             ( a : C_intarrs.Pointer; c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    g : Standard_Floating_Vectors.Vector(1..2);
    g_re,g_im : quad_double;
    gamma : QuadDobl_Complex_Numbers.Complex_Number;
    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    pwt : constant natural32 := natural32(v_a(v_a'first)); -- power of t

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Laurent_Homotopy_Gamma ...");
    end if;
    Assign(2,c,g);
    g_re := create(g(1));
    g_im := create(g(1));
    gamma := QuadDobl_Complex_Numbers.Create(g_re,g_im);
    if pwt > 0 then -- the power of t must be positive
      PHCpack_Operations.Create_QuadDobl_Homotopy(gamma,pwt);
    else
      if vrblvl > 0 then
        put("Invalid input "); put(pwt,1);
        put_line(" for the power of t.");
        put_line("Will resort to the default power of t ...");
      end if;
      PHCpack_Operations.Create_QuadDobl_Laurent_Homotopy(gamma);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Laurent_Homotopy_Gamma.");
      end if;
     return 923;
  end Path_Trackers_QuadDobl_Laurent_Homotopy_Gamma;

  function Path_Trackers_QuadDobl_Laurent_Homotopy_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Laurent_Homotopy_Clear ...");
    end if;
    QuadDobl_Laurent_Homotopy.Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Laurent_Homotopy_Clear.");
      end if;
     return 926;
  end Path_Trackers_QuadDobl_Laurent_Homotopy_Clear;

  function Path_Trackers_Standard_Silent_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ls : Standard_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Silent_Track ...");
    end if;
    PHCpack_Operations.Silent_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    Standard_Complex_Solutions.Clear(ls);
    if crash
     then return 155;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Silent_Track.");
      end if;
     return 155;
  end Path_Trackers_Standard_Silent_Track;

  function Path_Trackers_DoblDobl_Silent_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ls : DoblDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Silent_Track ...");
    end if;
    PHCpack_Operations.Silent_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    DoblDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 175;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Silent_Track.");
      end if;
     return 175;
  end Path_Trackers_DoblDobl_Silent_Track;

  function Path_Trackers_QuadDobl_Silent_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ls : QuadDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Silent_Track ...");
    end if;
    PHCpack_Operations.Silent_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    QuadDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 185;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Silent_Track.");
      end if;
     return 175;
  end Path_Trackers_QuadDobl_Silent_Track;

  function Path_Trackers_Standard_Report_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ls : Standard_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Report_Track ...");
    end if;
    PHCpack_Operations.Reporting_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    Standard_Complex_Solutions.Clear(ls);
    if crash
     then return 156;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Report_Track.");
      end if;
     return 156;
  end Path_Trackers_Standard_Report_Track;

  function Path_Trackers_DoblDobl_Report_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ls : DoblDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Report_Track ...");
    end if;
    PHCpack_Operations.Reporting_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    DoblDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 176;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Report_Track.");
      end if;
     return 176;
  end Path_Trackers_DoblDobl_Report_Track;

  function Path_Trackers_QuadDobl_Report_Track
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    ls : QuadDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    length : double_float;
    nbstep,nbfail,nbiter,nbsyst : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..4);
    crash : boolean;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Report_Track ...");
    end if;
    PHCpack_Operations.Reporting_Path_Tracker
      (ls,length,nbstep,nbfail,nbiter,nbsyst,crash);
    nbs(1) := nbstep; nbs(2) := nbfail;
    nbs(3) := nbiter; nbs(4) := nbsyst;
    Assign(nbs,a);
    Assign_Solution(ls,b,c);
    QuadDobl_Complex_Solutions.Clear(ls);
    if crash
     then return 186;
     else return 0;
    end if;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Report_Track.");
      end if;
     return 186;
  end Path_Trackers_QuadDobl_Report_Track;

  function Path_Trackers_Standard_Write_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Numbers,Standard_Continuation_Data_io;

    ls : constant Standard_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    nb : Standard_Natural_Vectors.Vector(1..5);
    length_path : constant double_float := IMAG_PART(ls.t);
    nbfail,nbregu,nbsing,kind : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    tol_sing : constant double_float := 1.0E-8;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Write_Solution ...");
    end if;
    Assign(5,a,nb);
    if PHCpack_Operations.Is_File_Defined then
      Write_Statistics(PHCpack_Operations.output_file,
                       integer32(nb(5)),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(PHCpack_Operations.output_file,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    else
      Write_Statistics(standard_output,
                       integer32(nb(5)),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(standard_output,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    end if;
    ls.t := Create(REAL_PART(ls.t),0.0);
    if PHCpack_Operations.Is_File_Defined
     then Write_Solution(PHCpack_Operations.output_file,ls.all,length_path);
     else Write_Solution(standard_output,ls.all,length_path);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Write_Solution.");
      end if;
     return 157;
  end Path_Trackers_Standard_Write_Solution;

  function Path_Trackers_DoblDobl_Write_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Numbers,DoblDobl_Continuation_Data_io;

    ls : constant DoblDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    nb : Standard_Natural_Vectors.Vector(1..5);
    length_path : constant double_double := IMAG_PART(ls.t);
    nbfail,nbregu,nbsing,kind : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    tol_sing : constant double_float := 1.0E-8;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Write_Solution ...");
    end if;
    Assign(5,a,nb);
    if PHCpack_Operations.Is_File_Defined then
      Write_Statistics(PHCpack_Operations.output_file,
                       nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(PHCpack_Operations.output_file,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    else
      Write_Statistics(standard_output,nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(standard_output,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    end if;
    ls.t := Create(REAL_PART(ls.t),Double_Double_Numbers.create(0.0));
    if PHCpack_Operations.Is_File_Defined
     then Write_Solution(PHCpack_Operations.output_file,ls.all,length_path);
     else Write_Solution(standard_output,ls.all,length_path);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Write_Solution.");
      end if;
     return 177;
  end Path_Trackers_DoblDobl_Write_Solution;

  function Path_Trackers_QuadDobl_Write_Solution
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Numbers,QuadDobl_Continuation_Data_io;

    ls : constant QuadDobl_Complex_Solutions.Link_to_Solution
       := Convert_to_Solution(b,c);
    nb : Standard_Natural_Vectors.Vector(1..5);
    length_path : constant quad_double := IMAG_PART(ls.t);
    nbfail,nbregu,nbsing,kind : natural32 := 0;
    tol_zero : constant double_float := 1.0E-12;
    tol_sing : constant double_float := 1.0E-8;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Write_Solution ...");
    end if;
    Assign(5,a,nb);
    if PHCpack_Operations.Is_File_Defined then
      Write_Statistics(PHCpack_Operations.output_file,
                       nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(PHCpack_Operations.output_file,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    else
      Write_Statistics(standard_output,nb(5),nb(1),nb(2),nb(3),nb(4));
      Write_Diagnostics(standard_output,
                        ls.all,tol_zero,tol_sing,nbfail,nbregu,nbsing,kind);
    end if;
    ls.t := Create(REAL_PART(ls.t),Quad_Double_Numbers.create(integer(0)));
    if PHCpack_Operations.Is_File_Defined
     then Write_Solution(PHCpack_Operations.output_file,ls.all,length_path);
     else Write_Solution(standard_output,ls.all,length_path);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Write_Solution.");
      end if;
     return 187;
  end Path_Trackers_QuadDobl_Write_Solution;

  function Path_Trackers_Standard_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks
    rtr : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Polynomial_Solve ...");
    end if;
    rtr := PHCpack_Operations.Solve_by_Standard_Homotopy_Continuation
             (nbt,vrblvl-1);
    return rtr;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Polynomial_Solve.");
      end if;
     return 16;
  end Path_Trackers_Standard_Polynomial_Solve;

  function Path_Trackers_DoblDobl_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks
    rtr : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Polynomial_Solve ...");
    end if;
    rtr := PHCpack_Operations.Solve_by_DoblDobl_Homotopy_Continuation
             (nbt,vrblvl-1);
    return rtr;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Polynomial_Solve.");
      end if;
     return 236;
  end Path_Trackers_DoblDobl_Polynomial_Solve;

  function Path_Trackers_QuadDobl_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbt : constant natural32 := natural32(v_a(v_a'first)); -- #tasks
    rtr : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Polynomial_Solve ...");
    end if;
    rtr := PHCpack_Operations.Solve_by_QuadDobl_Homotopy_Continuation
             (nbt,vrblvl-1);
    return rtr;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Polynomial_Solve.");
      end if;
      return 246;
  end Path_Trackers_QuadDobl_Polynomial_Solve;

  function Path_Trackers_Multprec_Polynomial_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    deci : constant natural32 := natural32(v_a(v_a'first));
    rtr : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Multprec_Polynomial_Solve ...");
    end if;
    rtr := PHCpack_Operations.Solve_by_Multprec_Homotopy_Continuation
             (deci,vrblvl-1);
    return rtr;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Multprec_Polynomial_Solve.");
      end if;
      return 496;
  end Path_Trackers_Multprec_Polynomial_Solve;

  function Path_Trackers_Standard_Laurent_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use PHCpack_Operations;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant natural32 := natural32(v_a(v_a'first));
    rtr : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_Standard_Laurent_Solve ...");
    end if;
    rtr := PHCpack_Operations.Solve_by_Standard_Laurent_Homotopy_Continuation
             (nbr,vrblvl-1);
    return rtr;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Laurent_Solve.");
      end if;
      return 774;
  end Path_Trackers_Standard_Laurent_Solve;

  function Path_Trackers_DoblDobl_Laurent_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use PHCpack_Operations;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant natural32 := natural32(v_a(v_a'first));
    rtr : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_DoblDobl_Laurent_Solve ...");
    end if;
    rtr := PHCpack_Operations.Solve_by_DoblDobl_Laurent_Homotopy_Continuation
             (nbr,vrblvl-1);
    return rtr;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Laurent_Solve.");
      end if;
      return 775;
  end Path_Trackers_DoblDobl_Laurent_Solve;

  function Path_Trackers_QuadDobl_Laurent_Solve
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use PHCpack_Operations;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant natural32 := natural32(v_a(v_a'first));
    rtr : integer32 := 0;

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface.");
      put_line("Path_Trackers_QuadDobl_Laurent_Solve ...");
    end if;
    rtr := PHCpack_Operations.Solve_by_QuadDobl_Laurent_Homotopy_Continuation
             (nbr,vrblvl-1);
    return rtr;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Laurent_Solve.");
      end if;
      return 776;
  end Path_Trackers_QuadDobl_Laurent_Solve;

  function Path_Trackers_Standard_Crude_Track
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    val : constant natural32 := natural32(v_a(v_a'first));
    verbose : constant boolean := (val = 1);

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface");
      put_line("Path_Trackers_Standard_Crude_Track ...");
    end if;
    Crude_Path_Trackers.Standard_Track_Paths(verbose);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_Standard_Crude_Track.");
      end if;
      return 622;
  end Path_Trackers_Standard_Crude_Track;

  function Path_Trackers_DoblDobl_Crude_Track
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    val : constant natural32 := natural32(v_a(v_a'first));
    verbose : constant boolean := (val = 1);

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface");
      put_line("Path_Trackers_DoblDobl_Crude_Track ...");
    end if;
    Crude_Path_Trackers.DoblDobl_Track_Paths(verbose);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_DoblDobl_Crude_Track.");
      end if;
      return 623;
  end Path_Trackers_DoblDobl_Crude_Track;

  function Path_Trackers_QuadDobl_Crude_Track
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    val : constant natural32 := natural32(v_a(v_a'first));
    verbose : constant boolean := (val = 1);

  begin
    if vrblvl > 0 then
      put("-> in path_trackers_interface");
      put_line("Path_Trackers_QuadDobl_Crude_Track ...");
    end if;
    Crude_Path_Trackers.QuadDobl_Track_Paths(verbose);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in path_trackers_interface.");
        put_line("Path_Trackers_QuadDobl_Crude_Track.");
      end if;
      return 624;
  end Path_Trackers_QuadDobl_Crude_Track;

end Path_Trackers_Interface;
