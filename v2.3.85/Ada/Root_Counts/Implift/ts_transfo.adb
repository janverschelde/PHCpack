with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer32_Transformations;
 use Standard_Integer32_Transformations;
with Standard_Integer_Transformations_io;
 use Standard_Integer_Transformations_io;

procedure ts_transfo is

-- DESCRIPTION :
--   Interactive test on some of the operations in Transformations.

  procedure Test_Building_Transformations ( n : in integer32 ) is

  -- DESCRIPION :
  --   Prompts the user for a vector with n components 
  --   and shows the results of various creators.

    v : Standard_Integer_Vectors.Vector(1..n);
    p : integer32;
    t1,t2,t3 : Transfo;
    
  begin
    put("Give "); put(n,1); put(" integers : "); get(v);
    put("-> your vector is "); put(v); new_line;
    put("Give the pivot : "); get(p);
    t1 := Create(v,p);
    put_line("The transformation after create : "); put(t1);
    t2 := Rotate(v,p);
    put_line("The transformation after rotate : "); put(t2);
    t3 := Build_Transfo(v,p);
    put_line("The transformation after build : "); put(t3);
  end Test_Building_Transformations;

  procedure Main is

    n : integer32 := 0;
    it : Transfo;

  begin
    new_line;
    put_line("Testing unimodular transformations ...");
    new_line;
    put("Give the number of variables : "); get(n);
    it := Id(natural32(n));
    put_line("The identity transformation : "); put(it);
    Test_Building_Transformations(n);
  end Main;

begin
  Main;
end ts_transfo;
