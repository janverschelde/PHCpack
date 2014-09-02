with Interfaces.C;
with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;    use Standard_Integer64_Matrices_io;
with Convex_Hull_Methods;               use Convex_Hull_Methods;
with Point_Lists_and_Strings;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with use_giftwrap;

procedure ts_use_giftwrap is

-- DESCRIPTION :
--   Test on the gateway interface to the giftwrapping method.

  procedure Show_Result ( a,b : in C_IntArrs.Pointer ) is

  -- DESCRIPTION :
  --   Extracts the string information of b.

    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    nbc : constant integer32 := integer32(v_a(v_a'first));
    v_b : constant C_Integer_Array(0..Interfaces.C.size_t(nbc))
       := C_Intarrs.Value(b,Interfaces.C.ptrdiff_t(nbc+1));
    sv : String(1..integer(nbc));

  begin
    put("The number of characters : "); put(nbc,1); new_line;
    sv := C_Integer_Array_to_String(natural32(nbc),v_b);
    put_line("The string representation : "); put_line(sv);
  end Show_Result;

  procedure Planar_Convex_Hull ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random point configuration of m points.

    A : Standard_Integer64_Matrices.Matrix(1..2,1..m) := Random_Data(2,m);
    s : constant string := Point_Lists_and_Strings.write(A);
    ar : C_Integer_Array(0..0);
    br : C_Integer_Array(0..Interfaces.C.size_t(s'last-1))
       := String_to_C_Integer_Array(natural32(s'last),s);
    ap,bp : C_IntArrs.Pointer;
    cp : C_DblArrs.Pointer;
    r : integer32;

  begin
    put_line("The point configuration : "); put(A);
    put_line("The string representation : "); put_line(s);
    ar(0) := Interfaces.C.int(s'last);
    ap := ar(0)'unchecked_access;
    bp := br(0)'unchecked_access;
    r := use_giftwrap(1,ap,bp,cp);
    Show_Result(ap,bp);
  end Planar_Convex_Hull;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of points
  --   and then calls the planar convex hull method.

    m : integer32 := 0;

  begin
    new_line;
    put("Give the number of points : "); get(m);
    Planar_Convex_Hull(m);
  end Main;

begin
  Main;
end ts_use_giftwrap;
