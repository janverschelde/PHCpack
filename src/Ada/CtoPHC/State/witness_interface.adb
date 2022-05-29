with text_io;                           use text_io;
with Interfaces.C;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Laur_Systems;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Laur_Systems;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_System_and_Solutions_io;
with DoblDobl_System_and_Solutions_io;
with QuadDobl_System_and_Solutions_io;
with Square_and_Embed_Systems;
with Witness_Sets_io;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Standard_PolySys_Container;
with Standard_LaurSys_Container;
with DoblDobl_PolySys_Container;
with DoblDobl_LaurSys_Container;
with QuadDobl_PolySys_Container;
with QuadDobl_LaurSys_Container;
with Standard_Solutions_Container;
with DoblDobl_Solutions_Container;
with QuadDobl_Solutions_Container;

package body Witness_Interface is

  function Witness_Standard_Polynomial_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Polynomial_Prompt ...");
    end if;
    Witness_Sets_io.Standard_Read_Embedding(p,sols,dim);
    Standard_PolySys_Container.Initialize(p.all);
    Standard_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Polynomial_Prompt.");
      end if;
      return 41;
  end Witness_Standard_Polynomial_Prompt;

  function Witness_DoblDobl_Polynomial_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Polynomial_Prompt ...");
    end if;
    Witness_Sets_io.DoblDobl_Read_Embedding(p,sols,dim);
    DoblDobl_PolySys_Container.Initialize(p.all);
    DoblDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_DoblDobl_Polynomial_Prompt.");
      end if;
      return 631;
  end Witness_DoblDobl_Polynomial_Prompt;

  function Witness_QuadDobl_Polynomial_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Polynomial_Prompt ...");
    end if;
    Witness_Sets_io.QuadDobl_Read_Embedding(p,sols,dim);
    QuadDobl_PolySys_Container.Initialize(p.all);
    QuadDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Polynomial_Prompt.");
      end if;
      return 661;
  end Witness_QuadDobl_Polynomial_Prompt;

  function Witness_Standard_Laurent_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    p : Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Laurent_Prompt ...");
    end if;
    Witness_Sets_io.Standard_Read_Embedding(p,sols,dim);
    Standard_LaurSys_Container.Initialize(p.all);
    Standard_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Laurent_Prompt.");
      end if;
      return 798;
  end Witness_Standard_Laurent_Prompt;

  function Witness_DoblDobl_Laurent_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    p : Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Laurent_Prompt ...");
    end if;
    Witness_Sets_io.DoblDobl_Read_Embedding(p,sols,dim);
    DoblDobl_LaurSys_Container.Initialize(p.all);
    DoblDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_DoblDobl_Laurent_Prompt.");
      end if;
      return 799;
  end Witness_DoblDobl_Laurent_Prompt;

  function Witness_QuadDobl_Laurent_Prompt
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    p : Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Laurent_Prompt ...");
    end if;
    Witness_Sets_io.QuadDobl_Read_Embedding(p,sols,dim);
    QuadDobl_LaurSys_Container.Initialize(p.all);
    QuadDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Laurent_Prompt.");
      end if;
      return 800;
  end Witness_QuadDobl_Laurent_Prompt;

  function Witness_Standard_Polynomial_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Polynomial_Read ...");
    end if;
    Open(file,in_file,filename);
    Witness_Sets_io.Standard_Read_Embedding(file,p,sols,dim);
    Standard_PolySys_Container.Initialize(p.all);
    Standard_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Polynomial_Read.");
      end if;
      return 64;
  end Witness_Standard_Polynomial_Read;

  function Witness_DoblDobl_Polynomial_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Polynomial_Read ...");
    end if;
    Open(file,in_file,filename);
    Witness_Sets_io.DoblDobl_Read_Embedding(file,p,sols,dim);
    DoblDobl_PolySys_Container.Initialize(p.all);
    DoblDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Polynomial_Read.");
      end if;
      return 654;
  end Witness_DoblDobl_Polynomial_Read;

  function Witness_QuadDobl_Polynomial_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    p : Link_to_Poly_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Polynomial_Read ...");
    end if;
    Open(file,in_file,filename);
    Witness_Sets_io.QuadDobl_Read_Embedding(file,p,sols,dim);
    QuadDobl_PolySys_Container.Initialize(p.all);
    QuadDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Polynomial_Read.");
      end if;
      return 684;
  end Witness_QuadDobl_Polynomial_Read;

  function Witness_Standard_Laurent_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    p : Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Laurent_Read ...");
    end if;
    Open(file,in_file,filename);
    Witness_Sets_io.Standard_Read_Embedding(file,p,sols,dim);
    Standard_LaurSys_Container.Initialize(p.all);
    Standard_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Laurent_Read.");
      end if;
      return 801;
  end Witness_Standard_Laurent_Read;

  function Witness_DoblDobl_Laurent_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    p : Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Laurent_Read ...");
    end if;
    Open(file,in_file,filename);
    Witness_Sets_io.DoblDobl_Read_Embedding(file,p,sols,dim);
    DoblDobl_LaurSys_Container.Initialize(p.all);
    DoblDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_DoblDobl_Laurent_Read.");
      end if;
      return 802;
  end Witness_DoblDobl_Laurent_Read;

  function Witness_QuadDobl_Laurent_Read
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    p : Link_to_Laur_Sys;
    sols : Solution_List;
    dim : natural32;
    data : Standard_Natural_Vectors.Vector(1..2);

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Laurent_Read ...");
    end if;
    Open(file,in_file,filename);
    Witness_Sets_io.QuadDobl_Read_Embedding(file,p,sols,dim);
    QuadDobl_LaurSys_Container.Initialize(p.all);
    QuadDobl_Solutions_Container.Initialize(sols);
    data(1) := dim;
    data(2) := Length_Of(sols);
    Assign(p'last,a);
    Assign(data,b);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Laurent_Read.");
      end if;
      return 803;
  end Witness_QuadDobl_Laurent_Read;

  function Witness_Standard_Polynomial_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_b(v_b'first));
    p : constant Standard_Complex_Poly_Systems.Link_to_Poly_Sys
      := Standard_PolySys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Polynomal_Swap ...");
    end if;
    Witness_Sets_io.Swap_Symbols_to_End(nvr,dim,"zz",p.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Polynomial_Swap.");
      end if;
      return 816;
  end Witness_Standard_Polynomial_Swap;

  function Witness_DoblDobl_Polynomial_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_b(v_b'first));
    p : constant DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := DoblDobl_PolySys_Container.Retrieve;
    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Polynomal_Swap ...");
    end if;
    Witness_Sets_io.Swap_Symbols_to_End(nvr,dim,"zz",p.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_DoblDobl_Polynomial_Swap.");
      end if;
      return 817;
  end Witness_DoblDobl_Polynomial_Swap;

  function Witness_QuadDobl_Polynomial_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_b(v_b'first));
    p : constant QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
      := QuadDobl_PolySys_Container.Retrieve;
    sols : QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Polynomal_Swap ...");
    end if;
    Witness_Sets_io.Swap_Symbols_to_End(nvr,dim,"zz",p.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Polynomial_Swap.");
      end if;
      return 818;
  end Witness_QuadDobl_Polynomial_Swap;

  function Witness_Standard_Laurent_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_b(v_b'first));
    p : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
      := Standard_LaurSys_Container.Retrieve;
    sols : Standard_Complex_Solutions.Solution_List
         := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Laurent_Swap ...");
    end if;
    Witness_Sets_io.Swap_Symbols_to_End(nvr,dim,"zz",p.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Laurent_Swap.");
      end if;
      return 819;
  end Witness_Standard_Laurent_Swap;

  function Witness_DoblDobl_Laurent_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_b(v_b'first));
    p : constant DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys
      := DoblDobl_LaurSys_Container.Retrieve;
    sols : DoblDobl_Complex_Solutions.Solution_List
         := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Laurent_Swap ...");
    end if;
    Witness_Sets_io.Swap_Symbols_to_End(nvr,dim,"zz",p.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_DoblDobl_Laurent_Swap.");
      end if;
      return 820;
  end Witness_DoblDobl_Laurent_Swap;

  function Witness_QuadDobl_Laurent_Swap
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    nvr : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_b(v_b'first));
    p : constant QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys
      := QuadDobl_LaurSys_Container.Retrieve;
    sols : QuadDobl_Complex_Solutions.Solution_List
         := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Laurent_Swap ...");
    end if;
    Witness_Sets_io.Swap_Symbols_to_End(nvr,dim,"zz",p.all,sols);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Laurent_Swap.");
      end if;
      return 821;
  end Witness_QuadDobl_Laurent_Swap;

  function Witness_Standard_Polynomial_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Poly_Systems;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    ep : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Polynomial_Embed ...");
    end if;
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    Standard_PolySys_Container.Clear;
    Standard_PolySys_Container.Initialize(ep.all);
   -- Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Polynomial_Embed.");
      end if;
      return 66;
  end Witness_Standard_Polynomial_Embed;

  function Witness_DoblDobl_Polynomial_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Poly_Systems;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    ep : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Polynomial_Embed ...");
    end if;
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    DoblDobl_PolySys_Container.Clear;
    DoblDobl_PolySys_Container.Initialize(ep.all);
   -- Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_DoblDobl_Polynomial_Embed.");
      end if;
      return 129;
  end Witness_DoblDobl_Polynomial_Embed;

  function Witness_QuadDobl_Polynomial_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Poly_Systems;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    ep : Link_to_Poly_Sys;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Polynomial_Embed ...");
    end if;
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    QuadDobl_PolySys_Container.Clear;
    QuadDobl_PolySys_Container.Initialize(ep.all);
   -- Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Polynomial_Embed.");
      end if;
      return 260;
  end Witness_QuadDobl_Polynomial_Embed;

  function Witness_Standard_Laurent_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Standard_Complex_Laur_Systems;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys := Standard_LaurSys_Container.Retrieve;
    ep : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Laurent_Embed ...");
    end if;
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    Standard_LaurSys_Container.Clear;
    Standard_LaurSys_Container.Initialize(ep.all);
   -- Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Polynomial_Embed.");
      end if;
      return 625;
  end Witness_Standard_Laurent_Embed;

  function Witness_DoblDobl_Laurent_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use DoblDobl_Complex_Laur_Systems;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys := DoblDobl_LaurSys_Container.Retrieve;
    ep : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Laurent_Embed ...");
    end if;
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    DoblDobl_LaurSys_Container.Clear;
    DoblDobl_LaurSys_Container.Initialize(ep.all);
   -- Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_DoblDobl_Laurent_Embed.");
      end if;
      return 626;
  end Witness_DoblDobl_Laurent_Embed;

  function Witness_QuadDobl_Laurent_Embed
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use QuadDobl_Complex_Laur_Systems;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    dim : constant natural32 := natural32(v_a(v_a'first));
    lp : constant Link_to_Laur_Sys := QuadDobl_LaurSys_Container.Retrieve;
    ep : Link_to_Laur_Sys;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Laurent_Embed ...");
    end if;
    Square_and_Embed_Systems.Square_and_Embed(lp.all,dim,ep);
    QuadDobl_LaurSys_Container.Clear;
    QuadDobl_LaurSys_Container.Initialize(ep.all);
   -- Witness_Sets_io.Add_Embed_Symbols(dim);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Laurent_Embed.");
      end if;
      return 627;
  end Witness_QuadDobl_Laurent_Embed;

  function Witness_Standard_Polynomial_Write
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    lp : constant Link_to_Poly_Sys := Standard_PolySys_Container.Retrieve;
    sols : constant Solution_List := Standard_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_Standard_Polynomial_Write ...");
    end if;
    Create(file,out_file,filename);
    Standard_System_and_Solutions_io.put(file,lp.all,sols);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Polynomial_Write.");
      end if;
      return 65;
  end Witness_Standard_Polynomial_Write;

  function Witness_DoblDobl_Polynomial_Write
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    lp : constant Link_to_Poly_Sys := DoblDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := DoblDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_DoblDobl_Polynomial_Write ...");
    end if;
    Create(file,out_file,filename);
    DoblDobl_System_and_Solutions_io.put(file,lp.all,sols);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_Standard_Polynomial_Write.");
      end if;
      return 655;
  end Witness_DoblDobl_Polynomial_Write;

  function Witness_QuadDobl_Polynomial_Write
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is 

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Solutions;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(1));
    n : constant natural32 := natural32(v_a(v_a'first));
    n1 : constant Interfaces.C.size_t := Interfaces.C.size_t(n-1);
    v_b : constant C_Integer_Array(0..n1)
        := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(n));
    filename : constant String(1..integer(n))
             := C_Integer_Array_to_String(n,v_b);
    file : file_type;
    lp : constant Link_to_Poly_Sys := QuadDobl_PolySys_Container.Retrieve;
    sols : constant Solution_List := QuadDobl_Solutions_Container.Retrieve;

  begin
    if vrblvl > 0 then
      put("-> in witness_interface.");
      put_line("Witness_QuadDobl_Polynomial_Write ...");
    end if;
    Create(file,out_file,filename);
    QuadDobl_System_and_Solutions_io.put(file,lp.all,sols);
    Close(file);
    return 0;
  exception
    when others => 
      if vrblvl > 0 then
        put("Exception raised in witness_interface.");
        put_line("Witness_QuadDobl_Polynomial_Write.");
      end if;
      return 685;
  end Witness_QuadDobl_Polynomial_Write;

end Witness_Interface;
