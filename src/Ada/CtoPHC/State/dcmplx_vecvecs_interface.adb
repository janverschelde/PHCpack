with Interfaces.C;
with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_IO;       use Standard_Integer_Numbers_IO;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_IO;       use Standard_Integer_Vectors_IO;
with Standard_Floating_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_IO;       use Standard_Complex_Vectors_IO;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with DCMPLX_VecVecs_Container;

package body DCMPLX_VecVecs_Interface is

  function DCMPLX_VecVecs_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_interface.Initialize ...");
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
        DCMPLX_VecVecs_Container.Initialize(nbr,vrblvl=>vrblvl-1);
      end;
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in dcmplx_vecvecs_interface.");
        put_line("Initialize.");
      end if;
      return -1;
  end DCMPLX_VecVecs_Initialize;

  function DCMPLX_VecVecs_Get_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    size : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_interface.Get_Dimension ...");
    end if;
    size := DCMPLX_VecVecs_Container.Size(vrblvl=>vrblvl-1);
    if vrblvl > 0
     then put("  size : "); put(size,1); new_line;
    end if;
    Assign(size,a);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in dcmplx_vecvecs_interface.");
        put_line("Get_Dimension.");
      end if;
      return -1;
  end DCMPLX_VecVecs_Get_Dimension;

  function DCMPLX_VecVecs_Get_Size
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    size : integer32;

  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_interface.Get_Size ...");
    end if;
    declare
      v_a : constant C_Integer_Array
          := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
      idx : constant integer32 := integer32(v_a(v_a'first));
    begin
      if vrblvl > 0 
       then put("  idx : "); put(idx,1); new_line;
      end if;
      size := DCMPLX_VecVecs_Container.Size(idx,vrblvl=>vrblvl-1);
    end;
    if vrblvl > 0
     then put("  size : "); put(size,1); new_line;
    end if;
    Assign(size,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in dcmplx_vecvecs_interface.");
        put_line("Get_Size.");
      end if;
      return -1;
  end DCMPLX_VecVecs_Get_Size;

  function DCMPLX_VecVecs_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_interface.Set ...");
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
          := C_dblarrs.Value(c,Interfaces.C.ptrdiff_t(2*dim));
      data : Standard_Complex_Vectors.Vector(1..dim);
      idx : Interfaces.C.size_t := 0;
    begin
      if vrblvl > 0 then
        put("  dim : "); put(dim,1);
        put(", i : "); put(i,1);
        put(", j : "); put(j,1); new_line;
      end if;
      for i in 1..dim loop
        data(i) := Create(double_float(v_c(idx)),double_float(v_c(idx+1)));
        if vrblvl > 0
         then put("idx : "); put(integer32(idx),1); new_line;
        end if;
        idx := idx + 2;
      end loop;
      if vrblvl > 0 then
        put_line("the complex vector :");
        put_line(data);
      end if;
      DCMPLX_VecVecs_Container.Store_Copy(i,j,data,vrblvl=>vrblvl-1);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in dcmplx_vecvecs_interface.");
        put_line("Set.");
      end if;
      return -1;
  end DCMPLX_VecVecs_Set;

  function DCMPLX_VecVecs_Get
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_interface.Get ...");
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
      data : Standard_Complex_Vectors.Link_to_Vector;
    begin
      if vrblvl > 0 then
        put("  dim : "); put(dim,1);
        put(", i : "); put(i,1);
        put(", j : "); put(j,1); new_line;
      end if;
      data := DCMPLX_VecVecs_Container.Get(i,j,vrblvl=>vrblvl-1);
      if vrblvl > 0 then
        put_line("the complex vector :");
        put_line(data.all);
        put("assigning to c ...");
      end if;
      declare
        fd : Standard_Floating_Vectors.Vector(1..2*dim);
      begin
        for i in 1..dim loop
          fd(2*i-1) := REAL_PART(data(i));
          fd(2*i)   := IMAG_PART(data(i));
        end loop;
        Assign(fd,c); -- Assign(data.all,c) does not work!
      end;
      if vrblvl > 0
       then put_line(" done assigning");
      end if;
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in dcmplx_vecvecs_interface.");
        put_line("Get.");
      end if;
      return -1;
  end DCMPLX_VecVecs_Get;

  function DCMPLX_VecVecs_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in dcmplx_vecvecs_interface.Clear ...");
    end if;
    DCMPLX_VecVecs_Container.Clear(vrblvl=>vrblvl-1);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in dcmplx_vecvecs_interface.");
        put_line("Clear.");
      end if;
      return -1;
  end DCMPLX_VecVecs_Clear;

end DCMPLX_VecVecs_Interface;
