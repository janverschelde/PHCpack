with Interfaces.C;
with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with String_Splitters;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with C_Integer_Arrays;                  use C_Integer_Arrays;
with C_Double_Arrays;                   use C_Double_Arrays;
with Assignments_in_Ada_and_C;          use Assignments_in_Ada_and_C;
with Parameter_Homotopy_State;
with use_sweep;

procedure ts_use_sweep is

-- DESCRIPTION :
--   Test on the gateway interface to the sweep homotopy.

  procedure Test_Numeric_Definition ( k,n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Interactive test on the numerical definition of the indices
  --   of the parameters.

    ar : C_Integer_Array(0..Interfaces.C.size_t(3));
    br : C_Integer_Array(0..Interfaces.C.size_t(m));
    ap,bp : C_IntArrs.Pointer;
    cp : C_DblArrs.Pointer := null;
    r,p : integer32 := 0;

  begin
    ar(0) := Interfaces.C.int(k);
    ar(1) := Interfaces.C.int(n);
    ar(2) := Interfaces.C.int(m);
    ap := ar(0)'unchecked_access;
    for i in 1..m loop
      put("Give index of parameter "); put(i,1); put(" : ");
      get(p);
      br(Interfaces.C.size_t(i-1)) := Interfaces.C.int(p);
    end loop;
    bp := br(0)'unchecked_access;
    r := use_sweep(0,ap,bp,cp);
    new_line;
    put_line("-> after storing the data via the use_sweep routine ...");
    put("The number of equations : ");
    put(Parameter_Homotopy_State.Get_Number_of_Equations,1); new_line;
    put("The number of variables : ");
    put(Parameter_Homotopy_State.Get_Number_of_Variables,1); new_line;
    put("The number of parameters : ");
    put(Parameter_Homotopy_State.Get_Number_of_Parameters,1); new_line;
    put("The indices of parameters : ");
    put(Parameter_Homotopy_State.Get_Indices.all); new_line;
  end Test_Numeric_Definition;

  procedure Test_Numeric_Retrieval is

  -- DESCRIPTION :
  --   Tests the numerical retrieval of the definitions via the use_sweep.

    ar1 : C_Integer_Array(0..Interfaces.C.size_t(1));
    ap,bp : C_IntArrs.Pointer := null;
    cp : C_DblArrs.Pointer := null;
    r,nbpar : integer32 := 0;

  begin
    new_line;
    put_line("-> retrieving definitions numerically via use_sweep ...");
    ap := ar1(0)'unchecked_access;
    r := use_sweep(2,ap,bp,cp);
    declare
      v_a : constant C_Integer_Array := C_intarrs.Value(ap);
      val : constant integer32 := integer32(v_a(v_a'first));
    begin
      put("The number of equations : "); put(val,1); new_line;
    end;
    ap := ar1(0)'unchecked_access;
    r := use_sweep(3,ap,bp,cp);
    declare
      v_a : constant C_Integer_Array := C_intarrs.Value(ap);
      val : constant integer32 := integer32(v_a(v_a'first));
    begin
      put("The number of variables : "); put(val,1); new_line;
    end;
    ap := ar1(0)'unchecked_access;
    r := use_sweep(4,ap,bp,cp);
    declare
      v_a : constant C_Integer_Array := C_intarrs.Value(ap);
      val : constant integer32 := integer32(v_a(v_a'first));
    begin
      put("The number of parameters : "); put(val,1); new_line;
      nbpar := val;
    end;
    declare
      ar2 : C_Integer_Array(0..Interfaces.C.size_t(nbpar));
    begin
      ap := ar2(0)'unchecked_access;
      r := use_sweep(5,ap,bp,cp);
      declare
        v_a : constant C_Integer_Array
            := C_intarrs.Value(ap,Interfaces.C.ptrdiff_t(nbpar));
        val : integer32;
      begin
        put("The indices :");
        for i in 1..nbpar loop
          val := integer32(v_a(Interfaces.C.size_t(i-1)));     
          put(" "); put(val,1);
        end loop;
        new_line;
      end;
    end;
  end Test_Numeric_Retrieval;

  procedure Test_Clear is

  -- DESCRIPTION :
  --   Performs the clear of the definitions via use_sweep.

    a,b : C_IntArrs.Pointer := null;
    c : C_DblArrs.Pointer := null;
    r : integer32 := 0;

  begin
    r := use_sweep(7,a,b,c);
  end Test_Clear;

  procedure Test_Symbolic_Definition ( k,n,m : in integer32 ) is

  -- DESCRIPTION :
  --   Initializes the symbol table first by asking for a polynomial
  --   and then initializes the indices to the parameters.

    p : Poly;

  begin
    Symbol_Table.Init(natural32(n));
    put("Give a polynomial in "); put(n,1);
    put_line(" variables ...");
    get(p);
    put("-> your polynomial : "); put(p); new_line;
    new_line;
    put("Reading the names of the "); put(m,1);
    put_line(" parameters ...");
    skip_line; -- skip newline from get(p)
    declare
      s : constant string := String_Splitters.Read_String;
      nc : constant integer32 := integer32(s'last);
      sv : constant Standard_Integer_Vectors.Vector
         := String_to_Integer_Vector(s);
      ar : C_Integer_Array(0..Interfaces.C.size_t(4));
      br : C_Integer_Array(0..Interfaces.C.size_t(nc));
      a,b : C_IntArrs.Pointer;
      c : C_DblArrs.Pointer := null;
      r : integer32 := 0;
    begin
      put_line("The string of parameter names : " & s);
      put("The number of characters : "); put(nc,1); new_line;
      ar(0) := Interfaces.C.int(k);
      ar(1) := Interfaces.C.int(n);
      ar(2) := Interfaces.C.int(m);
      ar(3) := Interfaces.C.int(nc);
      a := ar(0)'unchecked_access;
      b := br(0)'unchecked_access;
      Assign(sv,b);
      r := use_sweep(1,a,b,c);
    end;
  end Test_Symbolic_Definition;

  procedure Test_Symbolic_Retrieval is

    nbpar : constant integer32
          := Parameter_Homotopy_State.Get_Number_of_Parameters;
    ar : C_Integer_Array(0..Interfaces.C.size_t(1));
    br : C_Integer_Array(0..Interfaces.C.size_t(nbpar*10));
    a,b : C_IntArrs.Pointer;
    c : C_DblArrs.Pointer := null;
    r : integer32;

  begin
    new_line;
    put_line("-> retrieving definitions symbolically via use_sweep ...");
    a := ar(0)'unchecked_access;
    b := br(0)'unchecked_access;
    r := use_sweep(6,a,b,c);
    declare
      v_a : constant C_Integer_Array := C_intarrs.Value(a);
      val : constant integer32 := integer32(v_a(v_a'first));
      v_b : constant C_Integer_Array
          := C_intarrs.Value(b,Interfaces.C.ptrdiff_t(val));
      s : constant string := C_Integer_Array_to_String(natural32(val),v_b);
    begin
      put("The number of characters : "); put(val,1); new_line;
      put_line("The names of the parameters : " & s);
    end;
  end Test_Symbolic_Retrieval;

  procedure Standard_Test_Set_Start_Target_Values ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random values in standard double precision for start and
  --   target and passes those values to the use_sweep function.

    stv : constant Standard_Complex_Vectors.Vector
        := Standard_Random_Vectors.Random_Vector(1,m);
    lnk_stv : Standard_Complex_Vectors.Link_to_Vector;
    ar : C_Integer_Array(0..Interfaces.C.size_t(2));
    br : C_Integer_Array(0..Interfaces.C.size_t(1));
    cr : C_Double_Array(0..Interfaces.C.size_t(2*m));
    a,b : C_IntArrs.Pointer := null;
    c : C_DblArrs.Pointer := null;
    r : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Assigning the random complex values :");
    put_line(stv);
    new_line;
    put("Start or target values ? (0/1) ");
    Ask_Alternative(ans,"01");
    ar(0) := Interfaces.C.int(0); -- standard double precision
    if ans = '0'
     then ar(1) := Interfaces.C.int(0); -- start values
     else ar(1) := Interfaces.C.int(1); -- target values
    end if;
    a := ar(0)'unchecked_access;
    br(0) := Interfaces.C.int(m);
    b := br(0)'unchecked_access;
    c := cr(0)'unchecked_access;
    Assign(stv,c);
    r := use_sweep(8,a,b,c);
    new_line;
    put_line("After assigning via use_sweep :");
    if ans = '0'
     then lnk_stv := Parameter_Homotopy_State.Get_Start;
     else lnk_stv := Parameter_Homotopy_State.Get_Target;
    end if;
    put_line(lnk_stv.all);
  end Standard_Test_Set_Start_Target_Values;

  procedure DoblDobl_Test_Set_Start_Target_Values ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random values in double double precision for start and
  --   target and passes those values to the use_sweep function.

    ddv : constant DoblDobl_Complex_Vectors.Vector
        := DoblDobl_Random_Vectors.Random_Vector(1,m);
    lnk_ddv : DoblDobl_Complex_Vectors.Link_to_Vector;
    ar : C_Integer_Array(0..Interfaces.C.size_t(2));
    br : C_Integer_Array(0..Interfaces.C.size_t(1));
    cr : C_Double_Array(0..Interfaces.C.size_t(4*m));
    a,b : C_IntArrs.Pointer := null;
    c : C_DblArrs.Pointer := null;
    r : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Assigning the random complex values :");
    put_line(ddv);
    new_line;
    put("Start or target values ? (0/1) ");
    Ask_Alternative(ans,"01");
    ar(0) := Interfaces.C.int(1); -- double double precision
    if ans = '0'
     then ar(1) := Interfaces.C.int(0); -- start values
     else ar(1) := Interfaces.C.int(1); -- target values
    end if;
    a := ar(0)'unchecked_access;
    br(0) := Interfaces.C.int(m);
    b := br(0)'unchecked_access;
    c := cr(0)'unchecked_access;
    Assign(ddv,c);
    r := use_sweep(8,a,b,c);
    new_line;
    put_line("After assigning via use_sweep :");
    if ans = '0'
     then lnk_ddv := Parameter_Homotopy_State.Get_Start;
     else lnk_ddv := Parameter_Homotopy_State.Get_Target;
    end if;
    put_line(lnk_ddv.all);
  end DoblDobl_Test_Set_Start_Target_Values;

  procedure QuadDobl_Test_Set_Start_Target_Values ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Generates m random values in quad double precision for start and
  --   target and passes those values to the use_sweep function.

    qdv : constant QuadDobl_Complex_Vectors.Vector
        := QuadDobl_Random_Vectors.Random_Vector(1,m);
    lnk_qdv : QuadDobl_Complex_Vectors.Link_to_Vector;
    ar : C_Integer_Array(0..Interfaces.C.size_t(2));
    br : C_Integer_Array(0..Interfaces.C.size_t(1));
    cr : C_Double_Array(0..Interfaces.C.size_t(8*m));
    a,b : C_IntArrs.Pointer := null;
    c : C_DblArrs.Pointer := null;
    r : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Assigning the random complex values :");
    put_line(qdv);
    new_line;
    put("Start or target values ? (0/1) ");
    Ask_Alternative(ans,"01");
    ar(0) := Interfaces.C.int(2); -- quad double precision
    if ans = '0'
     then ar(1) := Interfaces.C.int(0); -- start values
     else ar(1) := Interfaces.C.int(1); -- target values
    end if;
    a := ar(0)'unchecked_access;
    br(0) := Interfaces.C.int(m);
    b := br(0)'unchecked_access;
    c := cr(0)'unchecked_access;
    Assign(qdv,c);
    r := use_sweep(8,a,b,c);
    new_line;
    put_line("After assigning via use_sweep :");
    if ans = '0'
     then lnk_qdv := Parameter_Homotopy_State.Get_Start;
     else lnk_qdv := Parameter_Homotopy_State.Get_Target;
    end if;
    put_line(lnk_qdv.all);
  end QuadDobl_Test_Set_Start_Target_Values;

  procedure Standard_Test_Get_Start_Target_Values ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Test the retrieval of the values in standard double precision
  --   for the parameters, start or target, via use_sweep.

    stv : Standard_Complex_Vectors.Vector(1..m);
    ar : C_Integer_Array(0..Interfaces.C.size_t(2));
    br : C_Integer_Array(0..Interfaces.C.size_t(1));
    cr : C_Double_Array(0..Interfaces.C.size_t(2*m));
    a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
    b : constant C_IntArrs.Pointer := br(0)'unchecked_access;
    c : constant C_DblArrs.Pointer := cr(0)'unchecked_access;
    r : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Start or target values ? (0/1) ");
    Ask_Alternative(ans,"01");
    ar(0) := Interfaces.C.int(0); -- standard double precision
    if ans = '0'
     then ar(1) := Interfaces.C.int(0); -- start values
     else ar(1) := Interfaces.C.int(1); -- target values
    end if;
    br(0) := Interfaces.C.int(m);
    r := use_sweep(9,a,b,c);
    Assign(natural32(2*m),c,stv);
    new_line;
    put_line("After retrieval vie use_sweep :");
    put_line(stv);
  end Standard_Test_Get_Start_Target_Values;

  procedure DoblDobl_Test_Get_Start_Target_Values ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Test the retrieval of the values in double double precision
  --   for the parameters, start or target, via use_sweep.

    stv : DoblDobl_Complex_Vectors.Vector(1..m);
    ar : C_Integer_Array(0..Interfaces.C.size_t(2));
    br : C_Integer_Array(0..Interfaces.C.size_t(1));
    cr : C_Double_Array(0..Interfaces.C.size_t(4*m));
    a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
    b : constant C_IntArrs.Pointer := br(0)'unchecked_access;
    c : constant C_DblArrs.Pointer := cr(0)'unchecked_access;
    r : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Start or target values ? (0/1) ");
    Ask_Alternative(ans,"01");
    ar(0) := Interfaces.C.int(1); -- double double precision
    if ans = '0'
     then ar(1) := Interfaces.C.int(0); -- start values
     else ar(1) := Interfaces.C.int(1); -- target values
    end if;
    br(0) := Interfaces.C.int(m);
    r := use_sweep(9,a,b,c);
    Assign(natural32(4*m),c,stv);
    new_line;
    put_line("After retrieval vie use_sweep :");
    put_line(stv);
  end DoblDobl_Test_Get_Start_Target_Values;

  procedure QuadDobl_Test_Get_Start_Target_Values ( m : in integer32 ) is

  -- DESCRIPTION :
  --   Test the retrieval of the values in quad double precision
  --   for the parameters, start or target, via use_sweep.

    stv : QuadDobl_Complex_Vectors.Vector(1..m);
    ar : C_Integer_Array(0..Interfaces.C.size_t(2));
    br : C_Integer_Array(0..Interfaces.C.size_t(1));
    cr : C_Double_Array(0..Interfaces.C.size_t(8*m));
    a : constant C_IntArrs.Pointer := ar(0)'unchecked_access;
    b : constant C_IntArrs.Pointer := br(0)'unchecked_access;
    c : constant C_DblArrs.Pointer := cr(0)'unchecked_access;
    r : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Start or target values ? (0/1) ");
    Ask_Alternative(ans,"01");
    ar(0) := Interfaces.C.int(2); -- quad double precision
    if ans = '0'
     then ar(1) := Interfaces.C.int(0); -- start values
     else ar(1) := Interfaces.C.int(1); -- target values
    end if;
    br(0) := Interfaces.C.int(m);
    r := use_sweep(9,a,b,c);
    Assign(natural32(8*m),c,stv);
    new_line;
    put_line("After retrieval vie use_sweep :");
    put_line(stv);
  end QuadDobl_Test_Get_Start_Target_Values;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the total number of variables.

    k,n,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the number of equations : "); get(k);
    put("Give the number of variables : "); get(n);
    put("Give the number of parameters : "); get(m);
    Test_Numeric_Definition(k,n,m);
    Test_Numeric_Retrieval;
    Test_Clear;
    new_line;
    put_line("-> after clearing the definitions ...");
    Test_Numeric_Retrieval;
    Test_Symbolic_Definition(k,n,m);
    Test_Numeric_Retrieval;
    Test_Symbolic_Retrieval;
    new_line;
    put("Standard, Double Double, or Quad Double ? (0/1/2) ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' =>
        Standard_Test_Set_Start_Target_Values(m);
        Standard_Test_Get_Start_Target_Values(m);
      when '1' =>
        DoblDobl_Test_Set_Start_Target_Values(m);
        DoblDobl_Test_Get_Start_Target_Values(m);
      when '2' =>
        QuadDobl_Test_Set_Start_Target_Values(m);
        QuadDobl_Test_Get_Start_Target_Values(m);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_use_sweep;
