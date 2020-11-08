with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Continuation_Parameters; 
with Continuation_Parameters_io; 
with Pack_Continuation_Parameters;
with Main_Poly_Continuation;            use Main_Poly_Continuation;
with PHCpack_Operations;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;

package body Continuation_Parameters_Interface is

  function Continuation_Parameters_Ask_Values
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Ask_Values ...");
    end if;
    if PHCpack_Operations.Is_File_Defined then
      Driver_for_Continuation_Parameters
        (PHCpack_Operations.output_file);
    else
      Driver_for_Continuation_Parameters;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameters_Ask_Values.");
      end if;
      return 70;
  end Continuation_Parameters_Ask_Values;

  function Continuation_Parameters_Ask_Output_Level
             ( vrblvl : integer32 := 0 ) return integer32 is

    oc : natural32;

  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Ask_Output_Level ...");
    end if;
    if PHCpack_Operations.Is_File_Defined then
      Driver_for_Process_io
        (PHCpack_Operations.output_file,oc);
    else
      Driver_for_Process_io(standard_output,oc);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameters_Ask_Output_Level.");
      end if;
      return 71;
  end Continuation_Parameters_Ask_Output_Level;

  function Continuation_Parameters_Get_All
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : constant Standard_Floating_Vectors.Vector(1..34)
      := Pack_Continuation_Parameters.Get;

  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Get_All ...");
    end if;
    Assign(v,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameters_Get_All.");
      end if;
      return 72;
  end Continuation_Parameters_Get_All;

  function Continuation_Parameters_Set_All
             ( c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v : Standard_Floating_Vectors.Vector(1..34);

  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Set_All ...");
    end if;
    Assign(34,c,v);
    Pack_Continuation_Parameters.Set(v);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameters_Set_All.");
      end if;
      return 73;
  end Continuation_Parameters_Set_All;

  function Continuation_Parameters_Autotune
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    level : constant natural32 := natural32(v_a(v_a'first));
    nbdgt : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Autotune ...");
    end if;
    Continuation_Parameters.Tune(level,nbdgt);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameters_Autotune.");
      end if;
      return 193;
  end Continuation_Parameters_Autotune;

  function Continuation_Parameters_Show
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Show ...");
    end if;
    Continuation_Parameters_io.put;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameter_Show.");
      end if;
      return 194;
  end Continuation_Parameters_Show;

  function Continuation_Parameters_Get_Value
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v_a(v_a'first));
    res : Standard_Floating_Vectors.Vector(1..1);

  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Get_Value ...");
    end if;
    if k = 0 or k > 34 then
      return 189;
    else
      res(1) := Pack_Continuation_Parameters.Get_Value(k);
      Assign(res,c);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameter_Get_Value.");
      end if;
      return 189;
  end Continuation_Parameters_Get_Value;

  function Continuation_Parameters_Set_Value
             ( a : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v_a(v_a'first));
    v_c : constant C_Double_Array
        := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(1));
    v : constant double_float := double_float(v_c(v_c'first));

  begin
    if vrblvl > 0 then
      put("-> in continuation_parameters_interface.");
      put_line("Continuation_Parameters_Set_Value ...");
    end if;
    if k = 0 or k > 34
     then return 190;
     else Pack_Continuation_Parameters.Set_Value(k,v);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in continuation_parameters_interface.");
        put_line("Continuation_Parameter_Set_Value.");
      end if;
      return 190;
  end Continuation_Parameters_Set_Value;

end Continuation_Parameters_Interface;
