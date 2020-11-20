with Interfaces.C;
with text_io;                            use text_io;
--with String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
--with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
--with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
--with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Multprec_Integer_Matrices;
--with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;
with Multprec_Lattice_Polygons;
with Multprec_Lattice_3d_Facets;
with Multprec_Lattice_4d_Facets;
with Standard_Complex_Laur_Systems;
--with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Standard_LaurSys_Container;
with Lists_of_Integer_Vectors;
with Supports_of_Polynomial_Systems;
with Standard_Initial_Forms;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;
with Multprec_Giftwrap_Container;
with Point_Lists_and_Strings;
with Facets_and_Strings;

package body Giftwrap_Interface is

  function String_of_Point_Configuration
             ( a : C_intarrs.Pointer; b : C_intarrs.Pointer ) return string is

  -- DESCRIPTION :
  --   Returns the string representation of the point configuration
  --   stored in b, where the number of characters is the value of a.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbc : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array(0..Interfaces.C.size_t(nbc))
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc+1));
    res : String(1..integer(nbc));

  begin
    res := C_Integer_Array_to_String(natural32(nbc),v_b);
    return res;
  end String_of_Point_Configuration;

  function Giftwrap_Planar_Hull
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    sv : constant string := String_of_Point_Configuration(a,b);
    r,c : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in giftwrap_interface.Giftwrap_Planar_Hull ...");
    end if;
   -- put_line("The string representation : "); put_line(sv);
    Point_Lists_and_Strings.Extract_Dimensions(sv,r,c);
    declare
      M : Multprec_Integer_Matrices.Matrix(1..r,1..c)
        := Point_Lists_and_Strings.parse(sv,r,c);
    begin
     -- put_line("The matrix : "); put(M);
      Multprec_Lattice_Polygons.Lexicographic_Decreasing_Sort(M);
      declare
        V : constant Multprec_Integer_Matrices.Matrix
          := Multprec_Lattice_Polygons.Convex_Hull_2D(M);
        vertices : constant string := Point_Lists_and_Strings.Write(V);
        N : constant Multprec_Integer_Matrices.Matrix
          := Multprec_Lattice_Polygons.Inner_Normals(V);
        normals : constant string := Point_Lists_and_Strings.write(N);
        result : constant string := "(" & vertices & ", " & normals & ")";
        resvec : constant Standard_Integer_Vectors.Vector
               := String_to_Integer_Vector(result);
      begin
       -- put_line("The vertices : "); put(V);
       -- put_line("The inner normals : "); put(N);
       -- put_line("The result : "); put_line(result);
        Assign(integer32(resvec'last),a);
        Assign(resvec,b);
      end;
      Multprec_Integer_Matrices.Clear(M);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_Planar_Hull.");
      end if;
      return 580;
  end Giftwrap_Planar_Hull;

  function Giftwrap_Spatial_Hull
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    sv : constant string := String_of_Point_Configuration(a,b);
    r,c : integer32;

  begin
    if vrblvl > 0 then
      put_line("-> in giftwrap_interface.Giftwrap_Spatial_Hull ...");
    end if;
   -- put_line("The string representation : "); put_line(sv);
    Point_Lists_and_Strings.Extract_Dimensions(sv,r,c);
    declare
      M : constant Multprec_Integer_Matrices.Matrix(1..r,1..c)
        := Point_Lists_and_Strings.parse(sv,r,c);
    begin
     -- put_line("The matrix : "); put(M);
      Multprec_Giftwrap_Container.Create(M);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_Spatial_Hull.");
      end if;
      return 581;
  end Giftwrap_Spatial_Hull;

  function Giftwrap_Number_of_Facets
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is
 
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    res : natural32 := 0;

  begin
    if vrblvl > 0 then
      put_line("-> in giftwrap_interface.Giftwrap_Number_of_Facets ...");
    end if;
    if dim = 3
     then res := Multprec_Giftwrap_Container.Number_of_3d_Facets;
     else res := Multprec_Giftwrap_Container.Number_of_4d_Facets;
    end if;
    Assign(integer32(res),b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_Number_of_Facets.");
      end if;
      return 582;
  end Giftwrap_Number_of_Facets;

  function Giftwrap_String_of_Facet
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    fcn : constant natural32 := natural32(v_b(v_b'first));

  begin
    if vrblvl > 0 then
      put_line("-> in giftwrap_interface.Giftwrap_String_of_Facet ...");
    end if;
    if dim = 3 then
      declare
        lft : constant Multprec_Lattice_3d_Facets.Link_to_3d_Facet
            := Multprec_Giftwrap_Container.Facet_3d_Data(fcn);
        pts : constant Multprec_Integer_Matrices.Link_to_Matrix
            := Multprec_Giftwrap_Container.Point_Configuration_in_3d;
        sfc : constant string := Facets_and_Strings.write(pts.all,lft.all);
        resvec : constant Standard_Integer_Vectors.Vector
               := String_to_Integer_Vector(sfc);
      begin
        Assign(integer32(resvec'last),a);
        Assign(resvec,b);
      end;
    elsif dim = 4 then
      declare
        lft : constant Multprec_Lattice_4d_Facets.Link_to_4d_Facet
            := Multprec_Giftwrap_Container.Facet_4d_Data(fcn);
        pts : constant Multprec_Integer_Matrices.Link_to_Matrix
            := Multprec_Giftwrap_Container.Point_Configuration_in_4d;
        sfc : constant string := Facets_and_Strings.write(pts.all,lft.all);
        resvec : constant Standard_Integer_Vectors.Vector
               := String_to_Integer_Vector(sfc);
      begin
        Assign(integer32(resvec'last),a);
        Assign(resvec,b);
      end;
    end if;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_String_of_Facet.");
      end if;
      return 583;
  end Giftwrap_String_of_Facet;

  function Giftwrap_String_Size
             ( a : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    idx : constant integer32 := integer32(v_a(v_a'first));
    lp : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
       := Standard_LaurSys_Container.retrieve;
    sup : Lists_of_Integer_Vectors.List;

  begin
    if vrblvl > 0
     then put_line("-> in giftwrap_interface.Giftwrap_String_Size ...");
    end if;
    sup := Supports_of_Polynomial_Systems.Create(lp(idx));
   -- put("The number of elements in the support : ");
   -- put(Lists_of_Integer_Vectors.Length_Of(sup),1); new_line;
    declare
      support : constant string := Point_Lists_and_Strings.Write(sup);
    begin
     -- put_line("The support as string : "); put_line(support);
      Multprec_Giftwrap_Container.Store_String(support);
      Assign(integer32(support'last),a);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_String_Size.");
      end if;
      return 586;
  end Giftwrap_String_Size;

  function Giftwrap_String_of_Support
             ( b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    support : constant string := Multprec_Giftwrap_Container.Retrieve_String;
    sv : constant Standard_Integer_Vectors.Vector
       := String_to_Integer_Vector(support);

  begin
    if vrblvl > 0 then
      put_line("-> in giftwrap_interface.Giftwrap_String_of_Support ...");
    end if;
   -- put_line("In job 8, the support string :"); put_line(support);
    Assign(sv,b);
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_String_of_Support.");
      end if;
      return 587;
  end GiftWrap_String_of_Support;

  function Giftwrap_String_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is 
  begin
    if vrblvl > 0 then
      put_line("-> in giftwrap_interface.Giftwrap_String_Clear ...");
    end if;
    Multprec_Giftwrap_Container.Clear_String;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_String_Clear.");
      end if;
      return 588;
  end Giftwrap_String_Clear;

  function Giftwrap_Laurent_Initial_Form
             ( a : C_intarrs.Pointer;
               b : C_intarrs.Pointer;
               vrblvl : integer32 := 0 ) return integer32 is

    use Interfaces.C;

    v_a : constant C_Integer_Array
        := C_intarrs.Value(a,Interfaces.C.ptrdiff_t(2));
    dim : constant integer32 := integer32(v_a(v_a'first));
    nbc : constant integer32 := integer32(v_a(v_a'first+1));
    v_b : constant C_Integer_Array(0..Interfaces.C.size_t(nbc))
        := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc+1));
    strnrm : String(1..integer(nbc));
    normal : Standard_Integer_Vectors.Vector(1..dim);
    lp : constant Standard_Complex_Laur_Systems.Link_to_Laur_Sys
       := Standard_LaurSys_Container.retrieve;

  begin
    if vrblvl > 0 then
      put_line("-> in giftwrap_interface.Giftwrap_Laurent_Initial_Form ...");
    end if;
    strnrm := C_Integer_Array_to_String(natural32(nbc),v_b);
   -- put("The dimension : "); put(dim,1); put_line(".");
   -- put("The normal in Ada : "); put(strnrm); put_line(".");
    normal := Point_Lists_and_Strings.parse(dim,strnrm);
   -- put("The parsed normal : "); put(normal); put_line(".");
    declare
      q : Standard_Complex_Laur_Systems.Laur_Sys(lp'range);
    begin
      q := Standard_Initial_Forms.Initial(lp.all,normal);
     -- put_line("The initial form : "); put_line(q);
      Standard_LaurSys_Container.Clear;
      Standard_LaurSys_Container.Initialize(q);
      Standard_Complex_Laur_Systems.Clear(q);
    end;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_Laurent_Initial_Form.");
      end if;
      return 589;
  end Giftwrap_Laurent_Initial_Form;

  function Giftwrap_3d_Clear
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in giftwrap_interface.Giftwrap_3d_Clear ...");
    end if;
    Multprec_Giftwrap_Container.Clear_3d;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_3d_Clear.");
      end if;
      return 584;
  end Giftwrap_3d_Clear;

  function Giftwrap_4d_Clear 
             ( vrblvl : integer32 := 0 ) return integer32 is
  begin
    if vrblvl > 0
     then put_line("-> in giftwrap_interface.Giftwrap_4d_Clear ...");
    end if;
    Multprec_Giftwrap_Container.Clear_4d;
    return 0;
  exception
    when others =>
      if vrblvl > 0 then
        put("Exception raised in giftwrap_interface.");
        put_line("Giftwrap_4d_Clear.");
      end if;
      return 585;
  end Giftwrap_4d_Clear;

end Giftwrap_Interface;
