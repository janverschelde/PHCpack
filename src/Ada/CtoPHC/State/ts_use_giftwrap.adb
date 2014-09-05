with Interfaces.C;
with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Multprec_Integer_Matrices;
with Multprec_Integer_Matrices_io;      use Multprec_Integer_Matrices_io;
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

  procedure Show_Number_of_Facets ( a,b : in C_IntArrs.Pointer ) is

  -- DESCRIPTION :
  --   Shows the number of facets as stored in b, of the convex
  --   hull in dimension equal to the value of a.
  --
    v_a : constant C_Integer_Array := C_intarrs.Value(a);
    dim : constant natural32 := natural32(v_a(v_a'first));
    v_b : constant C_Integer_Array := C_intarrs.Value(b);
    nfc : constant natural32 := natural32(v_b(v_b'first));

  begin
    put("The number of facets in ");
    put(dim,1); put("-space : "); put(nfc,1); new_line;
  end Show_Number_of_Facets;

  procedure Planar_Convex_Hull ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random point configuration of m points
  --   in the plane and then computes the convex hull.

    A : Multprec_Integer_Matrices.Matrix(1..2,1..m) := Random_Data(2,m);
    s : constant string := Point_Lists_and_Strings.write(A);
    ar : C_Integer_Array(0..0);
    cr : C_Double_Array(0..0);
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
    cp := cr(0)'unchecked_access;
    r := use_giftwrap(1,ap,bp,cp);
    Show_Result(ap,bp);
  end Planar_Convex_Hull;

  procedure Convex_Hull_in_3d ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random point configuration of m points
  --   in 3-space and computes its convex hull.

    A : Multprec_Integer_Matrices.Matrix(1..3,1..m) := Random_Data(3,m);
    s : constant string := Point_Lists_and_Strings.write(A);
    ar : C_Integer_Array(0..0);
    br : C_Integer_Array(0..Interfaces.C.size_t(s'last-1))
       := String_to_C_Integer_Array(natural32(s'last),s);
    cr : C_Double_Array(0..0);
    ap,bp : C_IntArrs.Pointer;
    cp : C_DblArrs.Pointer;
    r : integer32;

  begin
    put_line("The point configuration : "); put(A);
    put_line("The string representation : "); put_line(s);
    ar(0) := Interfaces.C.int(s'last);
    cr(0) := Interfaces.C.double(0.0);
    ap := ar(0)'unchecked_access;
    bp := br(0)'unchecked_access;
    cp := cr(0)'unchecked_access;
    put_line("Constructing the list of facets ...");
    r := use_giftwrap(2,ap,bp,cp);
    ar(0) := Interfaces.C.int(3);
    r := use_giftwrap(3,ap,bp,cp);
    Show_Number_of_Facets(ap,bp);
  end Convex_Hull_in_3d;

  procedure Convex_Hull_in_4d ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random point configuration of m points
  --   in 4-space and computes its convex hull.

    A : Multprec_Integer_Matrices.Matrix(1..4,1..m) := Random_Data(4,m);
    s : constant string := Point_Lists_and_Strings.write(A);
    ar : C_Integer_Array(0..0);
    br : C_Integer_Array(0..Interfaces.C.size_t(s'last-1))
       := String_to_C_Integer_Array(natural32(s'last),s);
    cr : C_Double_Array(0..0);
    ap,bp : C_IntArrs.Pointer;
    cp : C_DblArrs.Pointer;
    r : integer32;

  begin
    put_line("The point configuration : "); put(A);
    put_line("The string representation : "); put_line(s);
    ar(0) := Interfaces.C.int(s'last);
    cr(0) := Interfaces.C.double(0.0);
    ap := ar(0)'unchecked_access;
    bp := br(0)'unchecked_access;
    cp := cr(0)'unchecked_access;
    put_line("Constructing the list of facets ...");
    r := use_giftwrap(2,ap,bp,cp);
    ar(0) := Interfaces.C.int(4);
    r := use_giftwrap(3,ap,bp,cp);
    Show_Number_of_Facets(ap,bp);
  end Convex_Hull_in_4d;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of points
  --   and then calls the planar convex hull method.

    m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("What is the ambient dimension ? (2, 3, or 4) ");
    Ask_Alternative(ans,"234");
    new_line;
    put("Give the number of points : "); get(m);
    case ans is
      when '2' => Planar_Convex_Hull(m);
      when '3' => Convex_Hull_in_3d(m);
      when '4' => Convex_Hull_in_4d(m);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_use_giftwrap;
