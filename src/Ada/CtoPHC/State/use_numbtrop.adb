with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
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
with Numerical_Tropisms_Container;       use Numerical_Tropisms_Container;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

function use_numbtrop ( job : integer32;
                        a : C_intarrs.Pointer;
                        b : C_intarrs.Pointer;
                        c : C_dblarrs.Pointer ) return integer32 is

  function Job1 return integer32 is -- standard initialize

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
  end Job1;

  function Job2 return integer32 is -- dobldobl initialize

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
  end Job2;

  function Job3 return integer32 is -- quaddobl initialize

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
  end Job3;

  function Job4 return integer32 is -- store standard tropisms

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
    Assign(b,wnd);
    Assign(natural32(nfl),c,dre);
    for i in dir'range loop
      dir(i) := dre(i);
    end loop;
    err := dre(nfl);
    Numerical_Tropisms_Container.Store_Standard_Tropism(idx,wnd,dir,err);
    return 0;
  end Job4;

  function Job5 return integer32 is -- store dobldobl tropisms

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
  end Job5;

  function Job6 return integer32 is -- store quaddobl tropisms

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
  end Job6;

  function Job7 return integer32 is -- standard retrieve

    w : Standard_Integer_Vectors.Link_to_Vector;
    v : Standard_Floating_VecVecs.Link_to_VecVec;
    e : Standard_Floating_Vectors.Link_to_Vector;

  begin
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
  end Job7;

  function Job8 return integer32 is -- dobldobl retrieve

    w : Standard_Integer_Vectors.Link_to_Vector;
    v : Double_Double_VecVecs.Link_to_VecVec;
    e : Double_Double_Vectors.Link_to_Vector;

  begin
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
  end Job8;

  function Job9 return integer32 is -- quaddobl retrieve

    w : Standard_Integer_Vectors.Link_to_Vector;
    v : Quad_Double_VecVecs.Link_to_VecVec;
    e : Quad_Double_Vectors.Link_to_Vector;

  begin
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
  end Job9;

  function Job10 return integer32 is -- number of standard tropisms

    n : constant integer32 := Numerical_Tropisms_Container.Standard_Size;

  begin
    Assign(n,a);
    return 0;
  end Job10;

  function Job11 return integer32 is -- number of dobldobl tropisms

    n : constant integer32 := Numerical_Tropisms_Container.DoblDobl_Size;

  begin
    Assign(n,a);
    return 0;
  end Job11;

  function Job12 return integer32 is -- number of quaddobl tropisms

    n : constant integer32 := Numerical_Tropisms_Container.QuadDobl_Size;

  begin
    Assign(n,a);
    return 0;
  end Job12;

  function Job19 return integer32 is -- dimension of standard tropisms

    n : constant integer32 := Numerical_Tropisms_Container.Standard_Dimension;

  begin
    Assign(n,a);
    return 0;
  end Job19;

  function Job20 return integer32 is -- dimension of dobldobl tropisms

    n : constant integer32 := Numerical_Tropisms_Container.DoblDobl_Dimension;

  begin
    Assign(n,a);
    return 0;
  end Job20;

  function Job21 return integer32 is -- dimension of quaddobl tropisms

    n : constant integer32 := Numerical_Tropisms_Container.QuadDobl_Dimension;

  begin
    Assign(n,a);
    return 0;
  end Job21;

  function Job13 return integer32 is -- standard retrieve tropism

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
    Numerical_Tropisms_Container.Standard_Retrieve_Tropism(idx,wnd,dir,err);
    Assign(wnd,b);
    for i in dir'range loop
      cff(i) := dir(i);
    end loop;
    cff(cff'last) := err;
    Assign(cff,c);
    return 0;
  end Job13;

  function Job14 return integer32 is -- dobldobl retrieve tropism

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
    Numerical_Tropisms_Container.DoblDobl_Retrieve_Tropism(idx,wnd,dir,err);
    Assign(wnd,b);
    for i in dir'range loop
      cff(i) := dir(i);
    end loop;
    cff(cff'last) := err;
    Assign(cff,c);
    return 0;
  end Job14;

  function Job15 return integer32 is -- quaddobl retrieve tropism

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
    Numerical_Tropisms_Container.QuadDobl_Retrieve_Tropism(idx,wnd,dir,err);
    Assign(wnd,b);
    for i in dir'range loop
      cff(i) := dir(i);
    end loop;
    cff(cff'last) := err;
    Assign(cff,c);
    return 0;
  end Job15; 

  function Job16 return integer32 is -- standard clear
  begin
    Numerical_Tropisms_Container.Standard_Clear;
    return 0;
  end Job16;

  function Job17 return integer32 is -- dobldobl clear
  begin
    Numerical_Tropisms_Container.DoblDobl_Clear;
    return 0;
  end Job17;

  function Job18 return integer32 is -- quaddobl clear
  begin
    Numerical_Tropisms_Container.QuadDobl_Clear;
    return 0;
  end Job18;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 1 => return Job1; -- standard initialize
      when 2 => return Job2; -- dobldobl initialize
      when 3 => return Job3; -- quaddobl initialize
      when 4 => return Job4; -- store standard tropism
      when 5 => return Job5; -- store dobldobl tropism
      when 6 => return Job6; -- store quaddobl tropism
      when 7 => return Job7; -- standard retrieve
      when 8 => return Job8; -- dobldobl retrieve
      when 9 => return Job9; -- quaddobl retrieve
      when 10 => return Job10; -- number of standard double tropisms
      when 11 => return Job11; -- number of double double tropisms
      when 12 => return Job12; -- number of quad double tropisms
      when 13 => return Job13; -- standard retrieve tropism
      when 14 => return Job14; -- dobldobl retrieve tropism
      when 15 => return Job15; -- quaddobl retrieve tropism
      when 16 => return Job16; -- standard clear
      when 17 => return Job17; -- dobldobl clear
      when 18 => return Job18; -- quaddobl clear
      when 19 => return Job19; -- dimension of standard double tropisms
      when 20 => return Job20; -- dimension of double double tropisms
      when 21 => return Job21; -- dimension of quad double tropisms
      when others => put_line("  Sorry.  Invalid operation."); return -1;
    end case;
  exception
    when others => put("Exception raised in use_numbtrop handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_numbtrop;
