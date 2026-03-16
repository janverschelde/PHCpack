with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_IO;       use Standard_Integer_Vectors_IO;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_IO;      use Standard_Floating_Vectors_IO;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Double_VecVecs_Container;

package body Double_VecVecs_Interface is

  function Double_VecVecs_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_vecvecs_interface.Initialize ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      dim : constant integer32 := integer32(v_a(v_a'first));
    begin
      if vrblvl > 0
       then put("  dim : "); put(dim,1); put_line(" ...");
      end if;
      declare
        v_b : constant C_Integer_Array
            := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(dim));
        nbr : Standard_Integer_Vectors.Vector(1..dim);
      begin
        for i in 1..dim loop
          nbr(i) := integer32(v_b(Interfaces.C.size_t(i-1)));
        end loop;
        if vrblvl > 0
         then put("  sizes of the vectors : "); put(nbr); new_line;
        end if;
        Double_VecVecs_Container.Initialize(nbr,vrblvl=>vrblvl-1);
      end;
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_vecvecs_interface.");
        put_line("Initialize.");
      end if;
      return -1;
  end Double_VecVecs_Initialize;

  function Double_VecVecs_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    size : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in double_vecvecs_interface.Get_Dimension ...");
    end if;
    size := Double_VecVecs_Container.Size(vrblvl=>vrblvl-1);
    if vrblvl > 0
     then put("  size : "); put(size,1); new_line;
    end if;
    Assign(size,a);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_vecvecs_interface.");
        put_line("Get_Dimension.");
      end if;
      return -1;
  end Double_VecVecs_Get_Dimension;

  function Double_VecVecs_Get_Size_Array
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    size : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in double_vecvecs_interface.Get_Size_Array ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      idx : constant integer32 := integer32(v_a(v_a'first));
    begin
      if vrblvl > 0 
       then put("  idx : "); put(idx,1); new_line;
      end if;
      size := Double_VecVecs_Container.Size(idx,vrblvl=>vrblvl-1);
    end;
    if vrblvl > 0
     then put("  size : "); put(size,1); new_line;
    end if;
    Assign(size,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_vecvecs_interface.");
        put_line("Get_Size_Array.");
      end if;
      return -1;
  end Double_VecVecs_Get_Size_Array;

  function Double_VecVecs_Get_Size_Vector
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    size : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in double_vecvecs_interface.Get_Size_Vector ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
      idx1 : constant integer32 := integer32(v_a(v_a'first));
      use Interfaces.C;
      idx2 : constant integer32 := integer32(v_a(v_a'first+1));
    begin
      if vrblvl > 0 then
        put("  k idx : "); put(idx1,1);
        put(", i idx : "); put(idx2,1); new_line;
      end if;
      size := Double_VecVecs_Container.Size(idx1,idx2,vrblvl=>vrblvl-1);
    end;
    if vrblvl > 0
     then put("  size : "); put(size,1); new_line;
    end if;
    Assign(size,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_vecvecs_interface.");
        put_line("Get_Size_Vector.");
      end if;
      return -1;
  end Double_VecVecs_Get_Size_Vector;

  function Double_VecVecs_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_vecvecs_interface.Set ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      dim : constant integer32 := integer32(v_a(v_a'first));
      v_b : constant C_Integer_Array
          := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
      i : constant integer32 := integer32(v_b(v_b'first));
      use Interfaces.C;
      j : constant integer32 := integer32(v_b(v_b'first+1));
      v_c : constant C_Double_Array
          := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(dim));
      data : Standard_Floating_Vectors.Vector(1..dim);
    begin
      if vrblvl > 0 then
        put("  dim : "); put(dim,1);
        put(", i : "); put(i,1);
        put(", j : "); put(j,1); new_line;
      end if;
      for i in 1..dim loop
        data(i) := double_float(v_c(Interfaces.C.size_t(i-1)));
      end loop;
      if vrblvl > 0 then
        put_line("the floating-point vector :");
        put_line(data);
      end if;
      Double_VecVecs_Container.Store_Copy(i,j,data,vrblvl=>vrblvl-1);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_vecvecs_interface.");
        put_line("Set.");
      end if;
      return -1;
  end Double_VecVecs_Set;

  function Double_VecVecs_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_vecvecs_interface.Get ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      dim : constant integer32 := integer32(v_a(v_a'first));
      v_b : constant C_Integer_Array
          := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
      i : constant integer32 := integer32(v_b(v_b'first));
      use Interfaces.C;
      j : constant integer32 := integer32(v_b(v_b'first+1));
      data : Standard_Floating_Vectors.Link_to_Vector;
    begin
      if vrblvl > 0 then
        put("  dim : "); put(dim,1);
        put(", i : "); put(i,1);
        put(", j : "); put(j,1); new_line;
      end if;
      data := Double_VecVecs_Container.Get(i,j,vrblvl=>vrblvl-1);
      if vrblvl > 0 then
        put_line("the floating-point vector :");
        put_line(data.all);
      end if;
      Assign(data.all,c);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_vecvecs_interface.");
        put_line("Get.");
      end if;
      return -1;
  end Double_VecVecs_Get;

  function Double_VecVecs_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in double_vecvecs_interface.Clear ...");
    end if;
    Double_VecVecs_Container.Clear(vrblvl=>vrblvl-1);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in double_vecvecs_interface.");
        put_line("Clear.");
      end if;
      return -1;
  end Double_VecVecs_Clear;

end Double_VecVecs_Interface;
