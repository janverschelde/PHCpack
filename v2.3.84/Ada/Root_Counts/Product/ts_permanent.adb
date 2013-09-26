with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Matrices_io;       use Standard_Integer_Matrices_io;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Sets_of_Unknowns;                   use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with Degrees_in_Sets_of_Unknowns;        use Degrees_in_Sets_of_Unknowns;

procedure ts_permanent is

-- DESCRIPTION :
--   Test on computation of a permanent for the degree matrix for
--   a partition of the variables in a polynomial system.

  procedure Read_Partition_and_Compute_Permanent ( p : in Poly_Sys ) is

    n : constant integer32 := p'last;
    m : integer32 := 0;
 
  begin
    new_line;
    put("The number of variables : "); put(n,1); new_line;
    put("Give the number of sets in the partition : "); get(m);
    declare
      z : Partition(1..natural32(m));
      b : integer32;
      k : Vector(1..m);
      d : Matrix(1..n,1..m);
    begin
      put("Give "); put(m,1); put(" sets : ");
      Create(z,natural32(n)); get(z(1..natural32(m)));
      put("Your partition : "); put(z); new_line;
      for i in 1..m loop
        k(i) := integer32(Extent_Of(z(natural32(i))));
      end loop;
      put("The cardinalities in the partition : "); put(k); new_line;
      d := Degree_Table(p,z);
      put_line("The degree table : ");
      put(d);
      b := Per(d,k); -- b := Bezout_Number(p,z);
      put("The "); put(m,1); put("-homogeneous Bezout number is ");
      put(b,1); put_line(".");
    end;
  end Read_Partition_and_Compute_Permanent;

  procedure Main is

    lp : Link_to_Poly_Sys;

  begin
    new_line;
    put_line("Permanent computation of a degree matrix for a partition.");
    new_line;
    get(lp);
    Read_Partition_and_Compute_Permanent(lp.all);
  end Main;

begin
  Main;
end ts_permanent;
