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
with Standard_Power_Transformations;

procedure ts_powtrans is

-- DESCRIPTION :
--   Interactive testing of the operations in Power_Transformations.

  procedure Test_Operations ( n : in integer32 ) is

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

  procedure Main is

    n : integer32 := 0;
    ans : character;
 
  begin
    new_line;
    put_line("Testing power transformations ...");
    new_line;
    put("Give the ambient dimension : "); get(n);
    put("Use 64-bit arithmetic ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Test64_Operations(n);
     else Test_Operations(n);
    end if;
  end Main;

begin
  Main;
end ts_powtrans;

