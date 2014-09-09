with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;      use Standard_Integer_Matrices_io;
with Standard_Integer64_Vectors;
with Standard_Integer64_Vectors_io;     use Standard_Integer64_Vectors_io;
with Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;    use Standard_Integer64_Matrices_io;
with Multprec_Integer_Vectors;
with Multprec_Integer_Vectors_io;       use Multprec_Integer_Vectors_io;
with Multprec_Integer_Matrices;
with Multprec_Integer_Matrices_io;      use Multprec_Integer_Matrices_io;
with Standard_Power_Transformations;
with Multprec_Power_Transformations;

procedure ts_powtrans is

-- DESCRIPTION :
--   Interactive testing of the operations in Power_Transformations.

  procedure Test_Operations ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Tests n-dimensional unimodular operations in 32-bit integer arithmetic.

    use Standard_Integer_Vectors,Standard_Integer_Matrices;

    idm : constant Matrix(1..n,1..n)
        := Standard_Power_Transformations.Identity_Matrix(natural32(n));
    nrv : Vector(1..n);
    piv : integer32;
    elm,rot : Matrix(1..n,1..n);

  begin
    put_line("The identity matrix : "); put(idm);
    put("Give "); put(n,1);
    put(" components of a normal vector : "); get(nrv);
    piv := Standard_Power_Transformations.Pivot(nrv);
    put("The pivot : "); put(piv,1); new_line;
    if piv = 0 then
      put_line("A zero vector is not a normal vector!");
    else
      elm := Standard_Power_Transformations.Eliminate(nrv,piv);
      put_line("An eliminating power transformation :"); put(elm);
      rot := Standard_Power_Transformations.Rotate(nrv,piv);
      put_line("A rotating power transformation :"); put(rot);
    end if;
  end Test_Operations;

  procedure Test64_Operations ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Tests n-dimensional unimodular operations in 64-bit integer arithmetic.

    use Standard_Integer64_Vectors,Standard_Integer64_Matrices;

    idm : constant Matrix(1..n,1..n)
        := Standard_Power_Transformations.Identity_Matrix(natural32(n));
    nrv : Vector(1..n);
    piv : integer32;
    elm,rot : Matrix(1..n,1..n);

  begin
    put_line("The identity matrix : "); put(idm);
    put("Give "); put(n,1);
    put(" components of a normal vector : "); get(nrv);
    piv := Standard_Power_Transformations.Pivot(nrv);
    put("The pivot : "); put(piv,1); new_line;
    if piv = 0 then
      put_line("A zero vector is not a normal vector!");
    else
      elm := Standard_Power_Transformations.Eliminate(nrv,piv);
      put_line("An eliminating power transformation :"); put(elm);
      rot := Standard_Power_Transformations.Rotate(nrv,piv);
      put_line("A rotating power transformation :"); put(rot);
    end if;
  end Test64_Operations;

  procedure Test_Multprec_Operations ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Tests n-dimensional unimodular operations in multiprecision
  --   integer arithmetic.

    use Multprec_Integer_Vectors,Multprec_Integer_Matrices;

    idm : constant Matrix(1..n,1..n)
        := Multprec_Power_Transformations.Identity_Matrix(natural32(n));
    nrv : Vector(1..n);
    piv : integer32;
    elm,rot : Matrix(1..n,1..n);

  begin
    put_line("The identity matrix : "); put(idm);
    put("Give "); put(n,1);
    put(" components of a normal vector : "); get(nrv);
    piv := Multprec_Power_Transformations.Pivot(nrv);
    put("The pivot : "); put(piv,1); new_line;
    if piv = 0 then
      put_line("A zero vector is not a normal vector!");
    else
      elm := Multprec_Power_Transformations.Eliminate(nrv,piv);
      put_line("An eliminating power transformation :"); put(elm);
      rot := Multprec_Power_Transformations.Rotate(nrv,piv);
      put_line("A rotating power transformation :"); put(rot);
    end if;
  end Test_Multprec_Operations;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the ambient dimension and for the type
  --   of integer arithmetic.  Then launches the corresponding test.

    n : integer32 := 0;
    ans : character;
 
  begin
    new_line;
    put_line("Testing power transformations ...");
    new_line;
    put("Give the ambient dimension : "); get(n);
    new_line;
    put_line("MENU for the type of arithmetic : ");
    put_line("  0. integer 32-bit arithmetic;");
    put_line("  1. integer 64-bit arithmetic;");
    put_line("  2. multiprecision integer arithmetic.");
    put("Type 0, 1, or 2 to select the arithmetic : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Test_Operations(n);
      when '1' => Test64_Operations(n);
      when '2' => Test_Multprec_Operations(n);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_powtrans;

