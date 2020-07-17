with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Double_Double_Vectors;
with Double_Double_vecVecs;
with Quad_Double_Vectors;
with Quad_Double_vecVecs;
with Numerical_Tropisms_Container;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

package body Numerical_Tropisms_Interface is

  function Standard_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nbt : constant integer32 := integer32(v_a(v_a'first));
    dim : constant integer32 := integer32(v_a(v_a'first+1));
    nfl : constant integer32 := nbt*(dim+1);
    wnd : Standard_Integer_Vectors.Vector(1..nbt);
    dir : Standard_Floating_VecVecs.VecVec(1..nbt);
    dre : Standard_Floating_Vectors.Vector(1..nfl);
    err : Standard_Floating_Vectors.Vector(1..nbt);
    idx : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Standard_Initialize ...");
    end if;
    Assign(natural32(nbt),b,wnd);
    Assign(natural32(nfl),c,dre);
    idx := 0;
    for i in dir'range loop
      dir(i) := new Standard_Floating_Vectors.Vector(1..dim);
      for j in 1..dim loop
        idx := idx + 1;
        dir(i)(j) := dre(idx);
      end loop;
    end loop;
    for i in err'range loop
      err(i) := dre(nbt*dim+i);
    end loop;
    Numerical_Tropisms_Container.Standard_Initialize(wnd,dir,err);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.Standard_Initialize.");
      end if;
      return 711;
  end Standard_Initialize;

  function DoblDobl_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nbt : constant integer32 := integer32(v_a(v_a'first));
    dim : constant integer32 := integer32(v_a(v_a'first+1));
    nfl : constant integer32 := 2*nbt*(dim+1);
    wnd : Standard_Integer_Vectors.Vector(1..nbt);
    dre : Standard_Floating_Vectors.Vector(1..nfl);
    dir : Double_Double_VecVecs.VecVec(1..nbt);
    err : Double_Double_Vectors.Vector(1..nbt);
    idx : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.DoblDobl_Initialize ...");
    end if;
    Assign(natural32(nbt),b,wnd);
    Assign(natural32(nfl),c,dre);
    idx := 0;
    for i in dir'range loop
      dir(i) := new Double_Double_Vectors.Vector(1..dim);
      for j in 1..dim loop
        declare
          hi,lo : double_float;
        begin
          idx := idx + 1; hi := dre(idx);
          idx := idx + 1; lo := dre(idx);
          dir(i)(j) := Double_Double_Numbers.create(hi,lo);
        end;
      end loop;
    end loop;
    for i in err'range loop
      declare
        hi,lo : double_float;
      begin
        idx := idx + 1; hi := dre(idx);
        idx := idx + 1; lo := dre(idx);
        err(i) := Double_Double_Numbers.create(hi,lo);
      end;
    end loop;
    Numerical_Tropisms_Container.DoblDobl_Initialize(wnd,dir,err);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.DoblDobl_Initialize.");
      end if;
      return 712;
  end DoblDobl_Initialize;

  function QuadDobl_Initialize
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nbt : constant integer32 := integer32(v_a(v_a'first));
    dim : constant integer32 := integer32(v_a(v_a'first+1));
    nfl : constant integer32 := 4*nbt*(dim+1);
    wnd : Standard_Integer_Vectors.Vector(1..nbt);
    dre : Standard_Floating_Vectors.Vector(1..nfl);
    dir : Quad_Double_VecVecs.VecVec(1..nbt);
    err : Quad_Double_Vectors.Vector(1..nbt);
    idx : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.QuadDobl_Initialize ...");
    end if;
    Assign(natural32(nbt),b,wnd);
    Assign(natural32(nfl),c,dre);
    idx := 0;
    for i in dir'range loop
      dir(i) := new Quad_Double_Vectors.Vector(1..dim);
      for j in 1..dim loop
        declare
          hihi,lohi,hilo,lolo : double_float;
        begin
          idx := idx + 1; hihi := dre(idx);
          idx := idx + 1; lohi := dre(idx);
          idx := idx + 1; hilo := dre(idx);
          idx := idx + 1; lolo := dre(idx);
          dir(i)(j) := Quad_Double_Numbers.create(hihi,lohi,hilo,lolo);
        end;
      end loop;
    end loop;
    for i in err'range loop
      declare
        hihi,lohi,hilo,lolo : double_float;
      begin
        idx := idx + 1; hihi := dre(idx);
        idx := idx + 1; lohi := dre(idx);
        idx := idx + 1; hilo := dre(idx);
        idx := idx + 1; lolo := dre(idx);
        err(i) := Quad_Double_Numbers.create(hihi,lohi,hilo,lolo);
      end;
    end loop;
    Numerical_Tropisms_Container.QuadDobl_Initialize(wnd,dir,err);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.QuadDobl_Initialize.");
      end if;
      return 713;
  end QuadDobl_Initialize;

  function Store_Standard_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    dim : constant integer32 := integer32(v_a(v_a'first));
    idx : constant integer32 := integer32(v_a(v_a'first+1));
    wnd : integer32;
    nfl : constant integer32 := dim + 1;
    dre : Standard_Floating_Vectors.Vector(1..nfl);
    dir : Standard_Floating_Vectors.Vector(1..dim);
    err : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Store_Standard_Tropism ...");
    end if;
    Assign(b,wnd);
    Assign(natural32(nfl),c,dre);
    for i in dir'range loop
      dir(i) := dre(i);
    end loop;
    err := dre(nfl);
    Numerical_Tropisms_Container.Store_Standard_Tropism(idx,wnd,dir,err);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("Store_Standard_Tropism.");
      end if;
      return 714;
  end Store_Standard_Tropism;

  function Store_DoblDobl_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    dim : constant integer32 := integer32(v_a(v_a'first));
    idx : constant integer32 := integer32(v_a(v_a'first+1));
    wnd : integer32;
    nfl : constant integer32 := 2*dim + 2;
    dre : Standard_Floating_Vectors.Vector(1..nfl);
    dir : Double_Double_Vectors.Vector(1..dim);
    err : double_double;
    ind : integer32;
    hi,lo : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Store_DoblDobl_Tropism ...");
    end if;
    Assign(b,wnd);
    Assign(natural32(nfl),c,dre);
    ind := 0;
    for i in dir'range loop
      ind := ind + 1; hi := dre(ind);
      ind := ind + 1; lo := dre(ind);
      dir(i) := Double_Double_Numbers.create(hi,lo);
    end loop;
    ind := ind + 1; hi := dre(ind);
    ind := ind + 1; lo := dre(ind);
    err := Double_Double_Numbers.create(hi,lo);
    Numerical_Tropisms_Container.Store_DoblDobl_Tropism(idx,wnd,dir,err);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("Store_DoblDobl_Tropism.");
      end if;
      return 715;
  end Store_DoblDobl_Tropism;

  function Store_QuadDobl_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    dim : constant integer32 := integer32(v_a(v_a'first));
    idx : constant integer32 := integer32(v_a(v_a'first+1));
    wnd : integer32;
    nfl : constant integer32 := 4*dim + 4;
    dre : Standard_Floating_Vectors.Vector(1..nfl);
    dir : Quad_Double_Vectors.Vector(1..dim);
    err : quad_double;
    ind : integer32;
    hihi,lohi,hilo,lolo : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Store_QuadDobl_Tropism ...");
    end if;
    Assign(b,wnd);
    Assign(natural32(nfl),c,dre);
    ind := 0;
    for i in dir'range loop
      ind := ind + 1; hihi := dre(ind);
      ind := ind + 1; lohi := dre(ind);
      ind := ind + 1; hilo := dre(ind);
      ind := ind + 1; lolo := dre(ind);
      dir(i) := Quad_Double_Numbers.create(hihi,lohi,hilo,lolo);
    end loop;
    ind := ind + 1; hihi := dre(ind);
    ind := ind + 1; lohi := dre(ind);
    ind := ind + 1; hilo := dre(ind);
    ind := ind + 1; lolo := dre(ind);
    err := Quad_Double_Numbers.create(hihi,lohi,hilo,lolo);
    Numerical_Tropisms_Container.Store_QuadDobl_Tropism(idx,wnd,dir,err);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("Store_QuadDobl_Tropism.");
      end if;
      return 716;
  end Store_QuadDobl_Tropism;

  function Standard_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant integer32 := Numerical_Tropisms_Container.Standard_Size;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Standard_Size ...");
    end if;
    Assign(n,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.Standard_Size");
      end if;
      return 720;
  end Standard_Size;

  function DoblDobl_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant integer32 := Numerical_Tropisms_Container.DoblDobl_Size;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.DoblDobl_Size ...");
    end if;
    Assign(n,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.DoblDobl_Size");
      end if;
      return 721;
  end DoblDobl_Size;

  function QuadDobl_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant integer32 := Numerical_Tropisms_Container.QuadDobl_Size;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.QuadDobl_Size ...");
    end if;
    Assign(n,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.QuadDobl_Size");
      end if;
      return 722;
  end QuadDobl_Size;

  function Standard_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant integer32 := Numerical_Tropisms_Container.Standard_Dimension;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Standard_Dimension ...");
    end if;
    Assign(n,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.Standard_Dimension");
      end if;
      return 729;
  end Standard_Dimension;

  function DoblDobl_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant integer32 := Numerical_Tropisms_Container.DoblDobl_Dimension;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.DoblDobl_Dimension ...");
    end if;
    Assign(n,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.DoblDobl_Dimension");
      end if;
      return 730;
  end DoblDobl_Dimension;

  function QuadDobl_Dimension
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    n : constant integer32 := Numerical_Tropisms_Container.QuadDobl_Dimension;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.QuadDobl_Dimension ...");
    end if;
    Assign(n,a);
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.QuadDobl_Dimension");
      end if;
      return 731;
  end QuadDobl_Dimension;

  function Standard_Retrieve_One_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    dim : constant integer32 := integer32(v_a(v_a'first));
    idx : constant integer32 := integer32(v_a(v_a'first+1));
    wnd : integer32;
    dir : Standard_Floating_Vectors.Vector(1..dim);
    err : double_float;
    cff : Standard_Floating_Vectors.Vector(1..dim+1);

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Standard_Retrieve_One_Tropism ...");
    end if;
    Numerical_Tropisms_Container.Standard_Retrieve_Tropism(idx,wnd,dir,err);
    Assign(wnd,b);
    for i in dir'range loop
      cff(i) := dir(i);
    end loop;
    cff(cff'last) := err;
    Assign(cff,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("Standard_Retrieve_One_Tropism.");
      end if;
      return 723;
  end Standard_Retrieve_One_Tropism;

  function DoblDobl_Retrieve_One_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    dim : constant integer32 := integer32(v_a(v_a'first));
    idx : constant integer32 := integer32(v_a(v_a'first+1));
    wnd : integer32;
    dir : Double_Double_Vectors.Vector(1..dim);
    err : double_double;
    cff : Double_Double_Vectors.Vector(1..dim+1);

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.DoblDobl_Retrieve_One_Tropism ...");
    end if;
    Numerical_Tropisms_Container.DoblDobl_Retrieve_Tropism(idx,wnd,dir,err);
    Assign(wnd,b);
    for i in dir'range loop
      cff(i) := dir(i);
    end loop;
    cff(cff'last) := err;
    Assign(cff,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("DoblDobl_Retrieve_One_Tropism.");
      end if;
      return 724;
  end DoblDobl_Retrieve_One_Tropism;

  function QuadDobl_Retrieve_One_Tropism
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    dim : constant integer32 := integer32(v_a(v_a'first));
    idx : constant integer32 := integer32(v_a(v_a'first+1));
    wnd : integer32;
    dir : Quad_Double_Vectors.Vector(1..dim);
    err : quad_double;
    cff : Quad_Double_Vectors.Vector(1..dim+1);

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.QuadDobl_Retrieve_One_Tropism ...");
    end if;
    Numerical_Tropisms_Container.QuadDobl_Retrieve_Tropism(idx,wnd,dir,err);
    Assign(wnd,b);
    for i in dir'range loop
      cff(i) := dir(i);
    end loop;
    cff(cff'last) := err;
    Assign(cff,c);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("QuadDobl_Retrieve_One_Tropism.");
      end if;
      return 725;
  end QuadDobl_Retrieve_One_Tropism; 

  function Standard_Retrieve_All_Tropisms
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    w : Standard_Integer_Vectors.Link_to_Vector;
    v : Standard_Floating_VecVecs.Link_to_VecVec;
    e : Standard_Floating_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Standard_Retrieve_All_Tropisms ...");
    end if;
    Numerical_Tropisms_Container.Standard_Retrieve(w,v,e);
    declare
      nbt : constant integer32 := w'last;
      flv : constant Standard_Floating_Vectors.Link_to_Vector := v(v'first);
      dim : constant integer32 := flv'last;
      nbd : Standard_Integer_Vectors.Vector(1..2);
      csz : constant integer32 := nbt*(dim+1);
      cff : Standard_Floating_Vectors.Vector(1..csz);
      idx : integer32 := 0;
    begin
      nbd(1) := nbt;
      nbd(2) := dim;
      Assign(nbd,a);
      Assign(w.all,b);
      for i in v'range loop
        for j in v(i)'range loop
          idx := idx + 1;
          cff(idx) := v(i)(j);
        end loop;
      end loop;
      for i in e'range loop
        idx := idx + 1;
        cff(idx) := e(i);
      end loop;
      Assign(cff,c);
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("Standard_Retrieve_All_Tropisms.");
      end if;
      return 717;
  end Standard_Retrieve_All_Tropisms;

  function DoblDobl_Retrieve_All_Tropisms
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    w : Standard_Integer_Vectors.Link_to_Vector;
    v : Double_Double_VecVecs.Link_to_VecVec;
    e : Double_Double_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.DoblDobl_Retrieve_All_Tropisms ...");
    end if;
    Numerical_Tropisms_Container.DoblDobl_Retrieve(w,v,e);
    declare
      nbt : constant integer32 := w'last;
      flv : constant Double_Double_Vectors.Link_to_Vector := v(v'first);
      dim : constant integer32 := flv'last;
      nbd : Standard_Integer_Vectors.Vector(1..2);
      csz : constant integer32 := 2*nbt*(dim+1);
      cff : Standard_Floating_Vectors.Vector(1..csz);
      idx : integer32 := 0;
      hi,lo : double_float;
    begin
      nbd(1) := nbt;
      nbd(2) := dim;
      Assign(nbd,a);
      Assign(w.all,b);
      for i in v'range loop
        for j in v(i)'range loop
          hi := hi_part(v(i)(j)); idx := idx + 1; cff(idx) := hi;
          lo := lo_part(v(i)(j)); idx := idx + 1; cff(idx) := lo;
        end loop;
      end loop;
      for i in e'range loop
        hi := hi_part(e(i)); idx := idx + 1; cff(idx) := hi;
        lo := lo_part(e(i)); idx := idx + 1; cff(idx) := lo;
      end loop;
      Assign(cff,c);
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("DoblDobl_Retrieve_All_Tropisms.");
      end if;
      return 718;
  end DoblDobl_Retrieve_All_Tropisms;

  function QuadDobl_Retrieve_All_Tropisms
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               c : C_dblarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    w : Standard_Integer_Vectors.Link_to_Vector;
    v : Quad_Double_VecVecs.Link_to_VecVec;
    e : Quad_Double_Vectors.Link_to_Vector;

  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.QuadDobl_Retrieve_All_Tropisms ...");
    end if;
    Numerical_Tropisms_Container.QuadDobl_Retrieve(w,v,e);
    declare
      nbt : constant integer32 := w'last;
      flv : constant Quad_Double_Vectors.Link_to_Vector := v(v'first);
      dim : constant integer32 := flv'last;
      nbd : Standard_Integer_Vectors.Vector(1..2);
      csz : constant integer32 := 4*nbt*(dim+1);
      cff : Standard_Floating_Vectors.Vector(1..csz);
      idx : integer32 := 0;
      dfl : double_float;
    begin
      nbd(1) := nbt;
      nbd(2) := dim;
      Assign(nbd,a);
      Assign(w.all,b);
      for i in v'range loop
        for j in v(i)'range loop
          dfl := hihi_part(v(i)(j)); idx := idx + 1; cff(idx) := dfl;
          dfl := lohi_part(v(i)(j)); idx := idx + 1; cff(idx) := dfl;
          dfl := hilo_part(v(i)(j)); idx := idx + 1; cff(idx) := dfl;
          dfl := lolo_part(v(i)(j)); idx := idx + 1; cff(idx) := dfl;
        end loop;
      end loop;
      for i in e'range loop
        dfl := hihi_part(e(i)); idx := idx + 1; cff(idx) := dfl;
        dfl := lohi_part(e(i)); idx := idx + 1; cff(idx) := dfl;
        dfl := hilo_part(e(i)); idx := idx + 1; cff(idx) := dfl;
        dfl := lolo_part(e(i)); idx := idx + 1; cff(idx) := dfl;
      end loop;
      Assign(cff,c);
    end;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in use_numbtrop.");
        put_line("QuadDobl_Retrieve_All_Tropisms.");
      end if;
      return 719;
  end QuadDobl_Retrieve_All_Tropisms;

  function Standard_Clear ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.Standard_Clear ...");
    end if;
    Numerical_Tropisms_Container.Standard_Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.Standard_Clear.");
      end if;
      return 726;
  end Standard_Clear;

  function DoblDobl_Clear ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.DoblDobl_Clear ...");
    end if;
    Numerical_Tropisms_Container.DoblDobl_Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.DoblDobl_Clear.");
      end if;
      return 727;
  end DoblDobl_Clear;

  function QuadDobl_Clear ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in use_numbtrop.QuadDobl_Clear ...");
    end if;
    Numerical_Tropisms_Container.QuadDobl_Clear;
    return 0;
  exception
    when others => 
      if vrblvl > 0
       then put_line("Exception raised in use_numbtrop.QuadDobl_Clear.");
      end if;
      return 728;
  end QuadDobl_Clear;

end Numerical_Tropisms_Interface;
