with text_io;                           use text_io;
with Interfaces.C;                      use Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Natural_Vectors;
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
with Extrinsic_Diagonal_Homotopies;
with Extrinsic_Diagonal_Homotopies_io;
with Extrinsic_Diagonal_Solvers;
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
with PHCpack_Operations_io;

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

  function Diagonal_Homotopy_Symbols_Doubler
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Symbol_Table;
    use Extrinsic_Diagonal_Homotopies_io;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(3));
    n : constant integer32 := integer32(v_a(v_a'first));
    d : constant natural32 := natural32(v_a(v_a'first+1));
    nc : constant integer := integer(v_a(v_a'first+2));
    nc1 : constant Interfaces.C.size_t := Interfaces.C.size_t(nc-1);
    v_b : constant C_Integer_Array(0..nc1)
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nc));
    s : constant String(1..nc) := C_Integer_Array_to_String(natural32(nc),v_b);
    sa1 : Array_of_Symbols(1..n);
    nb2 : constant natural32 := Symbol_Table.Number;
    sa2e : constant Array_of_Symbols(1..integer32(nb2)) := Symbol_Table.Content;
    sa2 : constant Array_of_Symbols := Remove_Embed_Symbols(sa2e);
    s11 : Array_of_Symbols(sa1'range);
    s22 : constant Array_of_Symbols(sa2'range) := Add_Suffix(sa2,'2');
    ind : integer := 0;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Symbols_Doubler ...");
    end if;
   -- put_line("The string of names in Ada : " & s); 
    for i in 1..n loop
      declare
        sb : Symbol;
        ksb : integer;
      begin
        sb := (sb'range => ' ');
        ind := ind + 1;
        ksb := sb'first-1;
        while ind <= s'last loop
          exit when s(ind) = ' ';
          ksb := ksb + 1;
          sb(ksb) := s(ind);
          ind := ind + 1;
        end loop;
        sa1(i) := sb;
      end;
    end loop;
    s11 := Add_Suffix(sa1,'1');
   -- put("The symbols in the first array of symbols :");
   -- Write_Symbols(sa1);
   -- put("The symbols in the second array of symbols :");
   -- Write_Symbols(sa2);
   -- put("The first suffixed symbols :"); Write_Symbols(s11);
   -- put("The second suffixed symbols :"); Write_Symbols(s22);
    Symbol_Table.Clear;
    Assign_Symbol_Table(s11,s22);
    Witness_Sets_io.Add_Embed_Symbols(d);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Symbols_Doubler.");
      end if;
      return 230;
  end Diagonal_Homotopy_Symbols_Doubler;

  function Diagonal_Homotopy_Standard_Start_Solutions
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
      put_line("Diagonal_Homotopy_Standard_Start_Solutions ...");
    end if;
    PHCpack_Operations.Standard_Diagonal_Cascade_Solutions(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Standard_Start_Solutions.");
      end if;
      return 271;
  end Diagonal_Homotopy_Standard_Start_Solutions;

  function Diagonal_Homotopy_DoblDobl_Start_Solutions
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
      put_line("Diagonal_Homotopy_DoblDobl_Start_Solutions ...");
    end if;
    PHCpack_Operations.DoblDobl_Diagonal_Cascade_Solutions(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_DoblDobl_Start_Solutions.");
      end if;
      return 297;
  end Diagonal_Homotopy_DoblDobl_Start_Solutions;

  function Diagonal_Homotopy_QuadDobl_Start_Solutions
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
      put_line("Diagonal_Homotopy_QuadDobl_Start_Solutions ...");
    end if;
    PHCpack_Operations.QuadDobl_Diagonal_Cascade_Solutions(a_dim,b_dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_QuadDobl_Start_Solutions.");
      end if;
      return 298;
  end Diagonal_Homotopy_QuadDobl_Start_Solutions;

  function Diagonal_Homotopy_Prompt_Set
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v_a(v_a'first));
    n,dim,deg : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..2);
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Prompt_Set ...");
    end if;
    if k = 1 or k = 2 then
      PHCpack_Operations_io.Read_Witness_Set_for_Diagonal_Homotopy
        (k,n,dim,deg,fail);
      if fail then
        return 166;
      else
        Assign(integer32(n),a);
        nbs(1) := dim; nbs(2) := deg;
        Assign(nbs,b);
      end if;
    else
      put("Wrong value on input : "); put(k,1); new_line;
      return 166;
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Prompt_Set.");
      end if;
      return 166;
  end Diagonal_Homotopy_Prompt_Set;

  function Diagonal_Homotopy_Reset_Input
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    k : constant natural32 := natural32(v_a(v_a'first));
    deg,dim : natural32;
    nbs : Standard_Natural_Vectors.Vector(1..2);
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Reset_Input ...");
    end if;
   -- put("resetting the input file for witness set "); put(k,1); new_line;
    PHCpack_Operations_io.Reset_Witness_Input_File(k,deg,dim,fail);
    if fail  then
      return 167;
    else
      -- put("  degree : "); put(deg,1);
      -- put("  n : "); put(dim,1); new_line;
      nbs(1) := deg; nbs(2) := dim;
      Assign(nbs,b);
    end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Reset_Input.");
      end if;
      return 167;
  end Diagonal_Homotopy_Reset_Input;

  function Diagonal_Homotopy_Cascade_Dimension
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    n1 : constant natural32 := natural32(v_a(v_a'first));
    n2 : constant natural32 := natural32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(2));
    a_d : constant natural32 := natural32(v_b(v_b'first));
    b_d : constant natural32 := natural32(v_b(v_b'first+1));
    cd : natural32;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Cascade_Dimension ...");
    end if;
   -- put("  n1 = "); put(n1,1);
   -- put("  n2 = "); put(n2,1);
   -- put("  a = "); put(a_d,1);
   -- put("  b = "); put(b_d,1); new_line;
    if a_d >= b_d then
      cd := Extrinsic_Diagonal_Homotopies.Cascade_Dimension(n1,n2,a_d,b_d);
    else
      cd := Extrinsic_Diagonal_Homotopies.Cascade_Dimension(n2,n1,b_d,a_d);
    end if;
   -- put("cascade dimension cd = "); put(cd,1); new_line;
    Assign(integer32(cd),a);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Cascade_Dimension.");
      end if;
      return 168;
  end Diagonal_Homotopy_Cascade_Dimension;

  function Diagonal_Homotopy_Standard_Hyperset
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Hypersurface_Witdrivers;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant integer32 := integer32(v_a(v_a'first));
    n : constant integer := integer(v_a(v_a'first+1));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..n)
             := C_Integer_Array_to_String(natural32(n),v_b);
    file : file_type;
    fail : boolean;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Standard_Hyperset ...");
    end if;
   -- put("The number of the equation : "); put(k,1); new_line;
   -- if lp = null then
   --   put_line("But the systems container is empty!");
   -- elsif lp'last < k then
   --   put("But there are only "); put(lp'last,1);
   --   put_line(" polynomials in container!");
   -- else
   --   put("Polynomial : "); put(lp(k)); new_line;
   --   put_line("Creating the output file with name " & filename);
      Create(file,out_file,filename);
      Call_Root_Finder(file,lp(k),false,1.0E-10,fail);
      Close(file);
   -- end if;
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Standard_Hyperset.");
      end if;
      return 169;
  end Diagonal_Homotopy_Standard_Hyperset;

  function Diagonal_Homotopy_Standard_Collapse
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    d : constant natural32 := natural32(v_a(v_a'first+1));
    lp : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
       := Standard_PolySys_Container.Retrieve;
    sols : constant Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;
    clps : Standard_Complex_Solutions.Solution_List;
    cp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_Standard_Collapse ...");
    end if;
   -- put("#equations in systems container : "); put(lp'last,1); new_line;
   -- put("#solutions in solutions container : ");
   -- put(Standard_Complex_Solutions.Length_Of(sols),1); new_line;
    Standard_Complex_Solutions.Copy(sols,clps);
   -- put("Collapse system with k = ");
   -- put(k,1); put(" and d = "); put(d,1); new_line;
    Extrinsic_Diagonal_Solvers.Collapse_System(lp.all,clps,k,d,cp);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(cp.all);
    Standard_Solutions_Container.Clear;
    Standard_Solutions_Container.Initialize(clps);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_Standard_Collapse.");
      end if;
      return 170;
  end Diagonal_Homotopy_Standard_Collapse;

  function Diagonal_Homotopy_DoblDobl_Collapse
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    d : constant natural32 := natural32(v_a(v_a'first+1));
    lp : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := DoblDobl_PolySys_Container.Retrieve;
    sols : constant DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;
    clps : DoblDobl_Complex_Solutions.Solution_List;
    cp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_DoblDobl_Collapse ...");
    end if;
   -- put("#equations in systems container : "); put(lp'last,1); new_line;
   -- put("#solutions in solutions container : ");
   -- put(Length_Of(sols),1); new_line;
    DoblDobl_Complex_Solutions.Copy(sols,clps);
    Extrinsic_Diagonal_Solvers.Collapse_System(lp.all,clps,k,d,cp);
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(cp.all);
    DoblDobl_Solutions_Container.Clear;
    DoblDobl_Solutions_Container.Initialize(clps);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_DoblDobl_Collapse.");
      end if;
      return 299;
  end Diagonal_Homotopy_DoblDobl_Collapse;

  function Diagonal_Homotopy_QuadDobl_Collapse
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    k : constant natural32 := natural32(v_a(v_a'first));
    d : constant natural32 := natural32(v_a(v_a'first+1));
    lp : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
       := QuadDobl_PolySys_Container.Retrieve;
    sols : constant QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Solutions_Container.Retrieve;
    clps : QuadDobl_Complex_Solutions.Solution_List;
    cp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in diagonal_homotopy_interface.");
      put_line("Diagonal_Homotopy_QuadDobl_Collapse ...");
    end if;
   -- put("#equations in systems container : "); put(lp'last,1); new_line;
   -- put("#solutions in solutions container : ");
   -- put(Length_Of(sols),1); new_line;
    QuadDobl_Complex_Solutions.Copy(sols,clps);
    Extrinsic_Diagonal_Solvers.Collapse_System(lp.all,clps,k,d,cp);
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(cp.all);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(clps);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in diagonal_homotopy_interface.");
        put_line("Diagonal_Homotopy_QuadDobl_Collapse.");
      end if;
      return 312;
  end Diagonal_Homotopy_QuadDobl_Collapse;

end Diagonal_Homotopy_Interface;
