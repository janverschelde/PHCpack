with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Multprec_Integer_Matrices;
-- with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;
with Multprec_Lattice_Polygons;
with Multprec_Giftwrap_Container;
with Point_Lists_and_Strings;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

function use_giftwrap ( job : integer32;
                        a : C_intarrs.Pointer;
                        b : C_intarrs.Pointer;
                        c : C_dblarrs.Pointer ) return integer32 is

  function String_of_Point_Configuration return string is

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

  function Job1 return integer32 is -- planar convex hull

  -- DESCRIPTION :
  --   Parses the input data of a[0] integers in b as a string,
  --   constructs the planar convex hull and returns the vertices
  --   and inner normals in their string representation.

    sv : constant string := String_of_Point_Configuration;
    r,c : integer32;

  begin
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
        Assign(integer32(result'last),a);
        Assign(resvec,b);
      end;
      Multprec_Integer_Matrices.Clear(M);
    end;
    return 0;
  end Job1;

  function Job2 return integer32 is -- create convex hull in 3d or 4d

    sv : constant string := String_of_Point_Configuration;
    r,c : integer32;

  begin
   -- put_line("The string representation : "); put_line(sv);
    Point_Lists_and_Strings.Extract_Dimensions(sv,r,c);
    declare
      M : Multprec_Integer_Matrices.Matrix(1..r,1..c)
        := Point_Lists_and_Strings.parse(sv,r,c);
    begin
     -- put_line("The matrix : "); put(M);
      Multprec_Giftwrap_Container.Create(M);
    end;
    return 0;
  end Job2;

  function Job3 return integer32 is -- returns number of facets

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant integer32 := integer32(v_a(v_a'first));
    res : natural32 := 0;

  begin
    if dim = 3
     then res := Multprec_Giftwrap_Container.Number_of_3d_Facets;
     else res := Multprec_Giftwrap_Container.Number_of_4d_Facets;
    end if;
    Assign(integer32(res),b);
    return 0;
  end Job3;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 1 => return Job1; -- planar convex hull
      when 2 => return Job2; -- convex hull in 3d or 4d
      when 3 => return Job3; -- returns the number of facets
      when others => put_line("  Sorry.  Invalid operation."); return 1;
    end case;
  exception
    when others => put("Exception raised in unisolve handling job ");
                   put(job,1); put_line(".  Will not ignore."); raise;
  end Handle_Jobs;

begin
  return Handle_Jobs;
exception
  when others => put_line("Ignoring the exception, returning job number.");
                 return job;
end use_giftwrap;
