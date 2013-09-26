with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;          use Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with Witness_Sets;                      use Witness_Sets;
with Witness_Sets_io;                   use Witness_Sets_io;

procedure ts_coords is

-- DESCRIPTION :
--   Interactive development of moving to intrinsic coordinates.

  function Coefficients ( h : VecVec; n : integer32 ) return Matrix is

  -- DESCRIPTION :
  --   Returns the coefficients of the hyperplanes in the original
  --   space, without slack variables, before the embedding.

    res : Matrix(h'range,0..n);

  begin 
    for i in h'range loop
      for j in 0..n loop
        res(i,j) := h(i)(j);
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Max_Norm_of_Difference 
             ( n : in integer32; v1,v2 : in Vector ) return double_float is

  -- DESCRIPTION :
  --   Returns the max norm of the difference between the first n
  --   components of the two vectors v1 and v2.

    d : constant Vector := v1(1..n) - v2(1..n);

  begin
    return Max_Norm(d);
  end Max_Norm_of_Difference;

  procedure Compare ( n : in integer32; a,b : in Solution_List ) is

  -- DESCRIPTION :
  --   Compares the first n coordinates of the two solution lists.

    aa,bb : Link_to_Solution;
    a_tmp : Solution_List := a;
    b_tmp : Solution_List := b;
    cnt : integer32 := 0;

  begin
    put_line("Difference between the solution lists :");
    while not Is_Null(a_tmp) and not Is_Null(b_tmp) loop
      aa := Head_Of(a_tmp); bb := Head_Of(b_tmp);
      cnt := cnt + 1;
      put("  Solution "); put(cnt,1); put(" : ");
      put(Max_Norm_of_Difference(n,aa.v,bb.v)); new_line;
      a_tmp := Tail_Of(a_tmp); b_tmp := Tail_Of(b_tmp);
    end loop;
  end Compare;

  procedure Test_Project_and_Expand
              ( n : in integer32; esols : in Solution_List; 
                m : in Matrix; b : in Vector; v : in VecVec ) is

  -- DESCRIPTION :
  --   The projection and expansion of solution lists is tested
  --   for linear spaces in vector and matrix format.

    isols : Solution_List := Project(esols,b,v);
    asols : Solution_List := Expand(isols,b,v);

  begin
    put_line("Testing project and expand of solutions...");
   -- put_line("The list of intrinsic solutions : ");
   -- put(Standard_Output,Length_Of(isols),Head_Of(isols).n,isols);
   -- put_line("The solutions in the ambient space : ");
   -- put(Standard_Output,Length_Of(asols),Head_Of(asols).n,asols);
    Compare(n,esols,asols);
    Clear(isols); Clear(asols);
    put_line("Doing the project and expand using matrices...");
    isols := Project(esols,m);
    asols := Expand(isols,m);
   -- put_line("The list of intrinsic solutions : ");
   -- put(Standard_Output,Length_Of(isols),Head_Of(isols).n,isols);
   -- put_line("The solutions in the ambient space : ");
   -- put(Standard_Output,Length_Of(asols),Head_Of(asols).n,asols);
    Compare(n,esols,asols);
    Clear(isols); Clear(asols);
  end Test_Project_and_Expand;

  procedure Test_Sampling ( n,d,k : in integer32;
                            ep : in Poly_Sys; esols : in Solution_List ) is

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    p : constant Poly_Sys := Remove_Embedding1(ep,natural32(d));
    s : constant VecVec(1..d) := Slices(ep,natural32(d));
    b : Vector(1..n);
    v,w : VecVec(1..k);
    eqs : constant Matrix(1..d,0..n) := Coefficients(s,n);
    gen : constant Matrix(1..n,0..k) := Generators(eqs);
    ortgen : constant Matrix(1..n,0..k) := Orthogonalize(gen);

  begin
    new_line;
    put_line("The original polynomial system :");
    put(p);
    put_line("The coefficients of the slices :"); put(eqs,3);
    Generators(n,d,s,b,v);
    w := Orthogonalize(v);
    put_line("The offset vector :"); put_line(b);
    put_line("The directions : "); put(w);
    put_line("The plane in matrix format : "); put(ortgen,3);
    Test_Project_and_Expand(n,esols,ortgen,b,w);
  end Test_Sampling;

  procedure Main is

    ep : Link_to_Poly_Sys;
    sols : Solution_List;
    n,d,k : integer32 := 0;

  begin
    Standard_Read_Embedding(ep,sols,natural32(d));
    n := ep'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    Test_Sampling(n,d,k,ep.all,sols);
  end Main;

begin
  new_line;
  put_line("Moving to intrinsic coordinates...");
  Main;
end ts_coords;
