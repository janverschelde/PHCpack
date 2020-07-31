with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Symbol_Table; -- ,Symbol_Table_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Laurentials;
with Standard_Complex_Poly_Strings;
with Standard_Complex_Laur_Strings;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Laurentials;
with DoblDobl_Complex_Poly_Strings;
with DoblDobl_Complex_Laur_Strings;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Laurentials;
with QuadDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Laur_Strings;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Witness_Sets_io;
with Standard_Hypersurface_Witdrivers;
with DoblDobl_Hypersurface_Witdrivers;
with QuadDobl_Hypersurface_Witdrivers;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with QuadDobl_Solutions_Container;
with PHCpack_Operations;

package body Diagonal_Homotopy_Interface is

  function Diagonal_Homotopy_Standard_Polynomial_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(v_a(v_a'first));
    nc : constant integer := integer(v_a(v_a'first+1));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : Standard_Complex_Polynomials.Poly;
    eps : constant double_float := 1.0E-12;
    fail : boolean;
    e : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    esols : Standard_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Standard_Polynomial_Set ...");
    end if;
    if Symbol_Table.Empty
     then Symbol_Table.Init(nv);
    end if;
    p := Standard_Complex_Poly_Strings.Parse(nv,s);
    Silent_Root_Finder(p,eps,fail,e,esols);
   -- if fail
   --  then put_line("a failure occurred");
   --  else put_line("no failure occurred");
   -- end if;
    Witness_Sets_io.Add_Embed_Symbols(nv-1);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(e.all);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(esols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Standard_Polynomial_Set.");
      end if;
      return 270;
  end Diagonal_Homotopy_Standard_Polynomial_Set;

  function Diagonal_Homotopy_Standard_Laurential_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(v_a(v_a'first));
    nc : constant integer := integer(v_a(v_a'first+1));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : Standard_Complex_Laurentials.Poly;
    eps : constant double_float := 1.0E-12;
    fail : boolean;
    e : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
    esols : Standard_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Standard_Laurential_Set ...");
    end if;
    if Symbol_Table.Empty
     then Symbol_Table.Init(nv);
    end if;
    p := Standard_Complex_Laur_Strings.Parse(nv,s);
    Silent_Root_Finder(p,eps,fail,e,esols);
   -- if fail
   --  then put_line("a failure occurred");
   --  else put_line("no failure occurred");
   -- end if;
    Witness_Sets_io.Add_Embed_Symbols(nv-1);
    Standard_LaurSys_Container.Clear;
    Standard_LaurSys_Container.Initialize(e.all);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(esols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Standard_Laurential_Set.");
      end if;
      return 814;
  end Diagonal_Homotopy_Standard_Laurential_Set;

  function Diagonal_Homotopy_DoblDobl_Polynomial_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(v_a(v_a'first));
    nc : constant integer := integer(v_a(v_a'first+1));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : DoblDobl_Complex_Polynomials.Poly;
    eps : constant double_double := create(1.0E-12);
    fail : boolean;
    e : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    esols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_DoblDobl_Polynomial_Set ...");
    end if;
    if Symbol_Table.Empty
     then Symbol_Table.Init(nv);
    end if;
    p := DoblDobl_Complex_Poly_Strings.Parse(nv,s);
    Silent_Root_Finder(p,eps,fail,e,esols);
   -- if fail
   --  then put_line("a failure occurred");
   --  else put_line("no failure occurred");
   -- end if;
    Witness_Sets_io.Add_Embed_Symbols(nv-1);
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(e.all);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(esols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_DoblDobl_Polynomial_Set.");
      end if;
      return 259;
  end Diagonal_Homotopy_DoblDobl_Polynomial_Set;

  function Diagonal_Homotopy_DoblDobl_Laurential_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(v_a(v_a'first));
    nc : constant integer := integer(v_a(v_a'first+1));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : DoblDobl_Complex_Laurentials.Poly;
    eps : constant double_double := create(1.0E-12);
    fail : boolean;
    e : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    esols : DoblDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_DoblDobl_Laurential_Set ...");
    end if;
    if Symbol_Table.Empty
     then Symbol_Table.Init(nv);
    end if;
    p := DoblDobl_Complex_Laur_Strings.Parse(nv,s);
    Silent_Root_Finder(p,eps,fail,e,esols);
   -- if fail
   --  then put_line("a failure occurred");
   --  else put_line("no failure occurred");
   -- end if;
    Witness_Sets_io.Add_Embed_Symbols(nv-1);
    DoblDobl_LaurSys_Container.Clear;
    DoblDobl_LaurSys_Container.Initialize(e.all);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(esols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_DoblDobl_Laurential_Set.");
      end if;
      return 814;
  end Diagonal_Homotopy_DoblDobl_Laurential_Set;

  function Diagonal_Homotopy_QuadDobl_Polynomial_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(v_a(v_a'first));
    nc : constant integer := integer(v_a(v_a'first+1));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : QuadDobl_Complex_Polynomials.Poly;
    eps : constant quad_double := create(1.0E-12);
    fail : boolean;
    e : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    esols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_QuadDobl_Polynomial_Set ...");
    end if;
    if Symbol_Table.Empty
     then Symbol_Table.Init(nv);
    end if;
    p := QuadDobl_Complex_Poly_Strings.Parse(nv,s);
    Silent_Root_Finder(p,eps,fail,e,esols);
   -- if fail
   --  then put_line("a failure occurred");
   --  else put_line("no failure occurred");
   -- end if;
    Witness_Sets_io.Add_Embed_Symbols(nv-1);
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(e.all);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(esols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_QuadDobl_Polynomial_Set.");
      end if;
      return 269;
  end Diagonal_Homotopy_QuadDobl_Polynomial_Set;

  function Diagonal_Homotopy_QuadDobl_Laurential_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    nv : constant natural32 := natural32(v_a(v_a'first));
    nc : constant integer := integer(v_a(v_a'first+1));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    p : QuadDobl_Complex_Laurentials.Poly;
    eps : constant quad_double := create(1.0E-12);
    fail : boolean;
    e : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
    esols : QuadDobl_Complex_Solutions.Solution_List;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_QuadDobl_Laurential_Set ...");
    end if;
    if Symbol_Table.Empty
     then Symbol_Table.Init(nv);
    end if;
    p := QuadDobl_Complex_Laur_Strings.Parse(nv,s);
    Silent_Root_Finder(p,eps,fail,e,esols);
   -- if fail
   --  then put_line("a failure occurred");
   --  else put_line("no failure occurred");
   -- end if;
    Witness_Sets_io.Add_Embed_Symbols(nv-1);
    QuadDobl_LaurSys_Container.Clear;
    QuadDobl_LaurSys_Container.Initialize(e.all);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(esols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_QuadDobl_Laurential_Set.");
      end if;
      return 814;
  end Diagonal_Homotopy_QuadDobl_Laurential_Set;

  function Diagonal_Homotopy_Standard_Polynomial_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Standard_Make ...");
    end if;
    PHCpack_Operations.Standard_Diagonal_Homotopy(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Standard_Polynomial_Make.");
      end if;
      return 165;
  end Diagonal_Homotopy_Standard_Polynomial_Make;

  function Diagonal_Homotopy_DoblDobl_Polynomial_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_DoblDobl_Polynomial_Make ...");
    end if;
    PHCpack_Operations.DoblDobl_Diagonal_Homotopy(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_DoblDobl_Polynomial_Make.");
      end if;
      return 289;
  end Diagonal_Homotopy_DoblDobl_Polynomial_Make;

  function Diagonal_Homotopy_QuadDobl_Polynomial_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_QuadDobl_Polynomial_Make ...");
    end if;
    PHCpack_Operations.QuadDobl_Diagonal_Homotopy(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_QuadDobl_Polynomial_Make.");
      end if;
      return 290;
  end Diagonal_Homotopy_QuadDobl_Polynomial_Make;

  function Diagonal_Homotopy_Standard_Laurent_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Standard_Laurent_Make ...");
    end if;
    PHCpack_Operations.Standard_Diagonal_Laurent_Homotopy(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Standard_Laurent_Make.");
      end if;
      return 810;
  end Diagonal_Homotopy_Standard_Laurent_Make;

  function Diagonal_Homotopy_DoblDobl_Laurent_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_DoblDobl_Laurent_Make ...");
    end if;
    PHCpack_Operations.DoblDobl_Diagonal_Laurent_Homotopy(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_DoblDobl_Laurent_Make.");
      end if;
      return 811;
  end Diagonal_Homotopy_DoblDobl_Laurent_Make;

  function Diagonal_Homotopy_QuadDobl_Laurent_Make 
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    a_dim : constant natural32 := natural32(v_a(v_a'first));
    b_dim : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_QuadDobl_Laurent_Make ...");
    end if;
    PHCpack_Operations.QuadDobl_Diagonal_Laurent_Homotopy(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_QuadDobl_Laurent_Make.");
      end if;
      return 812;
  end Diagonal_Homotopy_QuadDobl_Laurent_Make;

end Diagonal_Homotopy_Interface;
