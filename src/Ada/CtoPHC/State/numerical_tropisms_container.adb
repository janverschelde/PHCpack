package body Numerical_Tropisms_Container is

-- INTERNAL DATA :

  wnd : Standard_Integer_Vectors.Link_to_Vector;
  st_dir : Standard_Floating_VecVecs.Link_to_VecVec;
  dd_dir : Double_Double_VecVecs.Link_to_VecVec;
  qd_dir : Quad_Double_VecVecs.Link_to_VecVec;
  st_err : Standard_Floating_Vectors.Link_to_Vector;
  dd_err : Double_Double_Vectors.Link_to_Vector;
  qd_err : Quad_Double_Vectors.Link_to_Vector;

-- CONSTRUTORS :

  procedure Standard_Initialize 
              ( w : in Standard_Integer_Vectors.Vector;
                v : in Standard_Floating_VecVecs.VecVec;
                e : in Standard_Floating_Vectors.Vector ) is
  begin
    wnd := new Standard_Integer_Vectors.Vector'(w);
    st_dir := new Standard_Floating_VecVecs.VecVec(v'range);
    for i in v'range loop
      st_dir(i) := new Standard_Floating_Vectors.Vector'(v(i).all);
    end loop;
    st_err := new Standard_Floating_Vectors.Vector'(e);
  end Standard_Initialize;

  procedure DoblDobl_Initialize 
              ( w : in Standard_Integer_Vectors.Vector;
                v : in Double_Double_VecVecs.VecVec;
                e : in Double_Double_Vectors.Vector ) is
  begin
    wnd := new Standard_Integer_Vectors.Vector'(w);
    dd_dir := new Double_Double_VecVecs.VecVec(v'range);
    for i in v'range loop
      dd_dir(i) := new Double_Double_Vectors.Vector'(v(i).all);
    end loop;
    dd_err := new Double_Double_Vectors.Vector'(e);
  end DoblDobl_Initialize;

  procedure QuadDobl_Initialize 
              ( w : in Standard_Integer_Vectors.Vector;
                v : in Quad_Double_VecVecs.VecVec;
                e : in Quad_Double_Vectors.Vector ) is
  begin
    wnd := new Standard_Integer_Vectors.Vector'(w);
    qd_dir := new Quad_Double_VecVecs.VecVec(v'range);
    for i in v'range loop
      qd_dir(i) := new Quad_Double_Vectors.Vector'(v(i).all);
    end loop;
    qd_err := new Quad_Double_Vectors.Vector'(e);
  end QuadDobl_Initialize;

  procedure Store_Standard_Tropism
              ( k : integer32; w : in integer32;
                v : in Standard_Floating_Vectors.Vector;
                e : in double_float ) is
  begin
    wnd(k) := w;
    for i in v'range loop
      st_dir(k)(i) := v(i);
    end loop;
    st_err(k) := e;
  end Store_Standard_Tropism;

  procedure Store_DoblDobl_Tropism
              ( k : integer32; w : in integer32;
                v : in Double_Double_Vectors.Vector;
                e : in double_double ) is
  begin
    wnd(k) := w;
    for i in v'range loop
      dd_dir(k)(i) := v(i);
    end loop;
    dd_err(k) := e;
  end Store_DoblDobl_Tropism;

  procedure Store_QuadDobl_Tropism
              ( k : integer32; w : in integer32;
                v : in Quad_Double_Vectors.Vector;
                e : in quad_double ) is
  begin
    wnd(k) := w;
    for i in v'range loop
      qd_dir(k)(i) := v(i);
    end loop;
    qd_err(k) := e;
  end Store_QuadDobl_Tropism;

-- SELECTORS :

  procedure Standard_Retrieve
              ( w : out Standard_Integer_Vectors.Link_to_Vector;
                v : out Standard_Floating_VecVecs.Link_to_VecVec;
                e : out Standard_Floating_Vectors.Link_to_Vector ) is
  begin
    w := wnd;
    v := st_dir;
    e := st_err;
  end Standard_Retrieve;

  procedure DoblDobl_Retrieve
              ( w : out Standard_Integer_Vectors.Link_to_Vector;
                v : out Double_Double_VecVecs.Link_to_VecVec;
                e : out Double_Double_Vectors.Link_to_Vector ) is
  begin
    w := wnd;
    v := dd_dir;
    e := dd_err;
  end DoblDobl_Retrieve;

  procedure QuadDobl_Retrieve
              ( w : out Standard_Integer_Vectors.Link_to_Vector;
                v : out Quad_Double_VecVecs.Link_to_VecVec;
                e : out Quad_Double_Vectors.Link_to_Vector ) is
  begin
    w := wnd;
    v := qd_dir;
    e := qd_err;
  end QuadDobl_Retrieve;

  function Standard_Size return integer32 is

    use Standard_Floating_VecVecs;

  begin
    if st_dir = null
     then return 0;
     else return st_dir'last;
    end if;
  end Standard_Size;

  function DoblDobl_Size return integer32 is

    use Double_Double_VecVecs;

  begin
    if dd_dir = null
     then return 0;
     else return dd_dir'last;
    end if;
  end DoblDobl_Size;

  function QuadDobl_Size return integer32 is

    use Quad_Double_VecVecs;

  begin
    if qd_dir = null
     then return 0;
     else return qd_dir'last;
    end if;
  end QuadDobl_Size;

  function Standard_Dimension return integer32 is

    use Standard_Floating_Vectors;
    use Standard_Floating_VecVecs;

  begin
    if st_dir = null
     then return 0;
     else return st_dir(st_dir'first)'last;
    end if;
  end Standard_Dimension;

  function DoblDobl_Dimension return integer32 is

    use Double_Double_Vectors;
    use Double_Double_VecVecs;

  begin
    if dd_dir = null
     then return 0;
     else return dd_dir(dd_dir'first)'last;
    end if;
  end DoblDobl_Dimension;

  function QuadDobl_Dimension return integer32 is

    use Quad_Double_Vectors;
    use Quad_Double_VecVecs;

  begin
    if qd_dir = null
     then return 0;
     else return qd_dir(qd_dir'first)'last;
    end if;
  end QuadDobl_Dimension;

  procedure Standard_Retrieve_Tropism
              ( k : in integer32; w : out integer32;
                v : out Standard_Floating_Vectors.Vector;
                e : out double_float ) is
  begin
    w := wnd(k);
    for i in v'range loop
      v(i) := st_dir(k)(i);
    end loop;
    e := st_err(k);
  end Standard_Retrieve_Tropism;

  procedure DoblDobl_Retrieve_Tropism
              ( k : in integer32; w : out integer32;
                v : out Double_Double_Vectors.Vector;
                e : out double_double ) is
  begin
    w := wnd(k);
    for i in v'range loop
      v(i) := dd_dir(k)(i);
    end loop;
    e := dd_err(k);
  end DoblDobl_Retrieve_Tropism;

  procedure QuadDobl_Retrieve_Tropism
              ( k : in integer32; w : out integer32;
                v : out Quad_Double_Vectors.Vector;
                e : out quad_double ) is
  begin
    w := wnd(k);
    for i in v'range loop
      v(i) := qd_dir(k)(i);
    end loop;
    e := qd_err(k);
  end QuadDobl_Retrieve_Tropism;

-- DESTRUCTORS :

  procedure Standard_Clear is
  begin
    Standard_Integer_Vectors.Clear(wnd);
    Standard_Floating_VecVecs.Deep_Clear(st_dir);
    Standard_Floating_Vectors.Clear(st_err);
  end Standard_Clear;

  procedure DoblDobl_Clear is
  begin
    Standard_Integer_Vectors.Clear(wnd);
    Double_Double_VecVecs.Deep_Clear(dd_dir);
    Double_Double_Vectors.Clear(dd_err);
  end DoblDobl_Clear;

  procedure QuadDobl_Clear is
  begin
    Standard_Integer_Vectors.Clear(wnd);
    Quad_Double_VecVecs.Deep_Clear(qd_dir);
    Quad_Double_Vectors.Clear(qd_err);
  end QuadDobl_Clear;

end Numerical_Tropisms_Container;
