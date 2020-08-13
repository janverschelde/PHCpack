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

end Witness_Interface;
