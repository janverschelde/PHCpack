with Interfaces.C;
with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Multprec_Integer_Matrices;
-- with Multprec_Integer_Matrices_io;       use Multprec_Integer_Matrices_io;
with Multprec_Lattice_Polygons;
with Point_Lists_and_Strings;
with Assignments_in_Ada_and_C;           use Assignments_in_Ada_and_C;

function use_giftwrap ( job : integer32;
                        a : C_intarrs.Pointer;
                        b : C_intarrs.Pointer;
                        c : C_dblarrs.Pointer ) return integer32 is

  function Job1 return integer32 is -- planar convex hull

  -- DESCRIPTION :
  --   Parses the input data of a[0] integers in b as a string.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbc : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array(0..Interfaces.C.size_t(nbc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc+1));
    sv : String(1..integer(nbc));
    r,c : integer32;

  begin
   -- put("The number of characters : "); put(nbc,1); new_line;
    sv := C_Integer_Array_to_String(natural32(nbc),v_b);
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
    end;
    return 0;
  end Job1;

  function Handle_Jobs return integer32 is
  begin
    case job is
      when 1 => return Job1; -- planar convex hull
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
