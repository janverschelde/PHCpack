with text_io;                           use text_io;
with Interfaces.C;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Assignments_in_Ada_and_C;
with C_to_Ada_Arrays;
with DEMiCs_Output_Data;

function use_outdata ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  function Allocate return integer32 is

  -- DESCRIPTION :
  --   Allocates memory for the lifting values of the supports.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbr : constant integer32 := integer32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nbr-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbr));
    crd : constant Standard_Integer_Vectors.Vector
        := C_to_Ada_Arrays.Convert(v_b);

  begin
    DEMiCs_Output_Data.Initialize_Lifting(crd);
    return 0;
  end Allocate;

  function Assign_Lifting return integer32 is

  -- DESCRIPTION :
  --   Assigns a lifting value.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    idxsup : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    idxpnt : constant integer32 := integer32(v_b(v_b'first));
    v_c : constant C_Double_Array := C_dblarrs.Value(c);
    value : constant double_float := double_float(v_c(v_c'first));

  begin
    DEMiCs_Output_Data.Assign_Lifting(idxsup,idxpnt,value);
    return 0;
  end Assign_Lifting;

  function Retrieve_Lifting return integer32 is

  -- DESCRIPTION :
  --   Retrieves a lifting value.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    idxsup : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    idxpnt : constant integer32 := integer32(v_b(v_b'first));
    value : double_float;

  begin
    value := DEMiCs_Output_Data.Retrieve_Lifting(idxsup,idxpnt);
    Assignments_in_Ada_and_C.Assign(value,c);
    return 0;
  end Retrieve_Lifting;

  function Clear_Lifting return integer32 is
  begin
    DEMiCs_Output_Data.Clear_Lifting;
    return 0;
  end Clear_Lifting;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Allocate;
      when 1 => return Assign_Lifting;
      when 2 => return Retrieve_Lifting;
      when 3 => return Clear_Lifting;
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in use_outdata handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others =>
     put_line("Ignoring the exception, returning job number.");
     return job;
end use_outdata;
