with text_io;                           use text_io;
with Interfaces.C;
with String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with Floating_Mixed_Subdivisions;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with Assignments_in_Ada_and_C;
with C_to_Ada_Arrays;
with Cells_Container;
with DEMiCs_Output_Data;
with Drivers_for_DEMiCs_Algorithm;

function use_outdata ( job : integer32;
                       a : C_intarrs.Pointer;
                       b : C_intarrs.Pointer;
                       c : C_dblarrs.Pointer ) return integer32 is

  function Allocate return integer32 is

  -- DESCRIPTION :
  --   Allocates memory for the lifting values of the supports.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant integer32 := integer32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nbr-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbr));
    crd : constant Standard_Integer_Vectors.Vector
        := C_to_Ada_Arrays.Convert(v_b);
    shifted : Standard_Integer_Vectors.Vector(1..nbr);

  begin
    for i in shifted'range loop
      shifted(i) := crd(i-1);
    end loop;
    DEMiCs_Output_Data.Initialize_Lifting(shifted);
    return 0;
  end Allocate;

  function Assign_Lifting return integer32 is

  -- DESCRIPTION :
  --   Assigns a lifting value.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    idxsup : constant integer32 := integer32(v_a(v_a'first))+1;
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    idxpnt : constant integer32 := integer32(v_b(v_b'first))+1;
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
    idxsup : constant integer32 := integer32(v_a(v_a'first))+1;
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    idxpnt : constant integer32 := integer32(v_b(v_b'first))+1;
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

  function Append_Cell_Indices return integer32 is

  -- DESCRIPTION :
  --   Retrieves the string from b and adds it to DEMiCs_Output_Data.

    use Assignments_in_Ada_and_C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nbr : constant integer32 := integer32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nbr-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbr));
    strcell : constant string
            := C_Integer_Array_to_String(natural32(nbr),v_b);

  begin
    DEMiCs_Output_Data.Add_Cell_Indices(strcell);
    return 0;
  end Append_Cell_Indices;

  function Retrieve_Cell_Indices return integer32 is

  -- DESCRIPTION :
  --   Retrieves the string corresponding to the index in a.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    idx : constant integer32 := integer32(v_a(v_a'first));
    ls : constant String_Splitters.Link_to_String
       := DEMiCs_Output_Data.Get_Cell_Indices(idx);
    nbr : integer32;

    use Assignments_in_Ada_and_C;
    use String_Splitters;

  begin
    if ls = null then
      nbr := 0;
    else
      nbr := integer32(ls'last);
      declare
        sv : constant Standard_Integer_Vectors.Vector
           := String_to_Integer_Vector(ls.all);
      begin
        Assign(sv,b);
      end;
    end if;
    Assign(nbr,a);
    return 0;
  end Retrieve_Cell_Indices;

  function Clear_Cell_Indices return integer32 is

  -- DESCRIPTION :
  --   Clears the list of strings that stores the indices to the cells.

  begin
    DEMiCs_Output_Data.Clear_Cell_Indices;
    return 0;
  end Clear_Cell_Indices;

  function Store_Mixed_Volume return integer32 is

  -- DESCRIPTION :
  --   Stores the mixed volume that is given in a.

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));

  begin
    DEMiCs_Output_Data.mixed_volume := integer32(v_a(v_a'first));
    return 0;
  end Store_Mixed_Volume;

  function Retrieve_Mixed_Volume return integer32 is

  -- DESCRIPTION :
  --   Retrieves the mixed volume and assigns it to a.

    mv : constant integer32 := DEMiCs_Output_Data.mixed_volume;

  begin
    Assignments_in_Ada_and_C.Assign(mv,a);
    return 0;
  end Retrieve_Mixed_Volume;

  function Mixed_Volume_by_DEMiCs return integer32 is 

  -- DESCRIPTION :
  --   Calls DEMiCs to ompute the mixed volume of the system in the
  --   standard systems container, or if that container is empty,
  --   of the system in the Laurent Systems container.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Laur_Systems;
    use Drivers_for_DEMiCs_Algorithm;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    lq : Link_to_Laur_Sys;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv : natural32;


  begin
    DEMiCs_Output_Data.Clear_Lifting;
    DEMiCs_Output_Data.Clear_Cell_Indices;
    DEMiCs_Output_Data.Clear_Allocated_Cells;
    if lp /= null then
      BlackBox_DEMiCs_Algorithm(lp.all,mix,lifsup,mixsub,mv);
    else
      lq := Standard_LaurSys_Container.Retrieve;
      BlackBox_DEMiCs_Algorithm(lq.all,mix,lifsup,mixsub,mv);
    end if;
    Assignments_in_Ada_and_C.Assign(integer32(mv),a);
    Cells_Container.Initialize(mix,lifsup,mixsub);
    return 0;
  end Mixed_Volume_by_DEMiCs;

  function Stable_Mixed_Volume_by_DEMiCs return integer32 is 

  -- DESCRIPTION :
  --   Calls DEMiCs to ompute the stable mixed volume of the system
  --   in the standard systems container, or if that container is empty,
  --   of the system in the Laurent Systems container.

    use Standard_Complex_Poly_Systems;
    use Drivers_for_DEMiCs_Algorithm;

    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    mix : Standard_Integer_Vectors.Link_to_Vector;
    lifsup : Arrays_of_Floating_Vector_Lists.Link_to_Array_of_Lists;
    mixsub : Floating_Mixed_Subdivisions.Mixed_Subdivision;
    mv,smv,tmv : natural32 := 0;

  begin
    DEMiCs_Output_Data.Clear_Lifting;
    DEMiCs_Output_Data.Clear_Cell_Indices;
    DEMiCs_Output_Data.Clear_Allocated_Cells;
    if lp /= null then
      BlackBox_DEMiCs_Algorithm(lp.all,mix,lifsup,mixsub,mv,smv,tmv);
      Cells_Container.Initialize(mix,lifsup,mixsub);
    end if;
    Assignments_in_Ada_and_C.Assign(integer32(mv),a);
    Assignments_in_Ada_and_C.Assign(integer32(smv),b);
    return 0;
  end Stable_Mixed_Volume_by_DEMiCs;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 0 => return Allocate;
      when 1 => return Assign_Lifting;
      when 2 => return Retrieve_Lifting;
      when 3 => return Clear_Lifting;
      when 4 => return Append_Cell_Indices;
      when 5 => return Retrieve_Cell_Indices;
      when 6 => return Clear_Cell_Indices;
      when 7 => return Store_Mixed_Volume;
      when 8 => return Retrieve_Mixed_Volume;
      when 9 => return Mixed_Volume_by_DEMiCs;
      when 10 => return Stable_Mixed_Volume_by_DEMiCs;
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
