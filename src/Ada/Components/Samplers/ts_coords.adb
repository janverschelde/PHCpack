with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_VecVecs;
with Standard_Complex_VecVecs_io;       use Standard_Complex_VecVecs_io;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;      use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecVecs_io;       use DoblDobl_Complex_VecVecs_io;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Matrices_io;      use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs_io;       use QuadDobl_Complex_VecVecs_io;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices_io;      use QuadDobl_Complex_Matrices_io;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;
with Standard_Plane_Representations;    use Standard_Plane_Representations;
with Standard_Intrinsic_Solutions;      use Standard_Intrinsic_Solutions;
with DoblDobl_Plane_Representations;    use DoblDobl_Plane_Representations;
with DoblDobl_Intrinsic_Solutions;      use DoblDobl_Intrinsic_Solutions;
with QuadDobl_Plane_Representations;    use QuadDobl_Plane_Representations;
with QuadDobl_Intrinsic_Solutions;      use QuadDobl_Intrinsic_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;

procedure ts_coords is

-- DESCRIPTION :
--   Interactive development of moving to intrinsic coordinates.

  function Coefficients
             ( h : Standard_Complex_VecVecs.VecVec; n : integer32 ) 
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coefficients of the hyperplanes in the original
  --   space, without slack variables, before the embedding.

    res : Standard_Complex_Matrices.Matrix(h'range,0..n);

  begin 
    for i in h'range loop
      for j in 0..n loop
        res(i,j) := h(i)(j);
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Coefficients
             ( h : DoblDobl_Complex_VecVecs.VecVec; n : integer32 ) 
             return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coefficients of the hyperplanes in the original
  --   space, without slack variables, before the embedding.

    res : DoblDobl_Complex_Matrices.Matrix(h'range,0..n);

  begin 
    for i in h'range loop
      for j in 0..n loop
        res(i,j) := h(i)(j);
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Coefficients
             ( h : QuadDobl_Complex_VecVecs.VecVec; n : integer32 ) 
             return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coefficients of the hyperplanes in the original
  --   space, without slack variables, before the embedding.

    res : QuadDobl_Complex_Matrices.Matrix(h'range,0..n);

  begin 
    for i in h'range loop
      for j in 0..n loop
        res(i,j) := h(i)(j);
      end loop;
    end loop;
    return res;
  end Coefficients;

  function Max_Norm_of_Difference 
             ( n : in integer32;
               v1,v2 : in Standard_Complex_Vectors.Vector )
             return double_float is

  -- DESCRIPTION :
  --   Returns the max norm of the difference between the first n
  --   components of the two vectors v1 and v2.

    use Standard_Complex_Vectors;

    d : constant Vector := v1(1..n) - v2(1..n);

  begin
    return Max_Norm(d);
  end Max_Norm_of_Difference;

  function Max_Norm_of_Difference 
             ( n : in integer32;
               v1,v2 : in DoblDobl_Complex_Vectors.Vector )
             return double_double is

  -- DESCRIPTION :
  --   Returns the max norm of the difference between the first n
  --   components of the two vectors v1 and v2.

    use DoblDobl_Complex_Vectors;

    d : constant Vector := v1(1..n) - v2(1..n);

  begin
    return Max_Norm(d);
  end Max_Norm_of_Difference;

  function Max_Norm_of_Difference 
             ( n : in integer32;
               v1,v2 : in QuadDobl_Complex_Vectors.Vector )
             return quad_double is

  -- DESCRIPTION :
  --   Returns the max norm of the difference between the first n
  --   components of the two vectors v1 and v2.

    use QuadDobl_Complex_Vectors;

    d : constant Vector := v1(1..n) - v2(1..n);

  begin
    return Max_Norm(d);
  end Max_Norm_of_Difference;

  procedure Compare ( n : in integer32;
                      a,b : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Compares the first n coordinates of the two solution lists.

    use Standard_Complex_Solutions;

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

  procedure Compare ( n : in integer32;
                      a,b : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Compares the first n coordinates of the two solution lists.

    use DoblDobl_Complex_Solutions;

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

  procedure Compare ( n : in integer32;
                      a,b : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Compares the first n coordinates of the two solution lists.

    use QuadDobl_Complex_Solutions;

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
              ( n : in integer32;
                esols : in Standard_Complex_Solutions.Solution_List; 
                m : in Standard_Complex_Matrices.Matrix;
                b : in Standard_Complex_Vectors.Vector;
                v : in Standard_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The projection and expansion of solution lists is tested
  --   for linear spaces in vector and matrix format,
  --   in standard double precision.

    use Standard_Complex_Solutions;

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

  procedure Test_Project_and_Expand
              ( n : in integer32;
                esols : in DoblDobl_Complex_Solutions.Solution_List; 
                m : in DoblDobl_Complex_Matrices.Matrix;
                b : in DoblDobl_Complex_Vectors.Vector;
                v : in DoblDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The projection and expansion of solution lists is tested
  --   for linear spaces in vector and matrix format,
  --   in double double precision.

    use DoblDobl_Complex_Solutions;

    isols : Solution_List := Project(esols,b,v);
    asols : Solution_List := Expand(isols,b,v);

  begin
    put_line("Testing project and expand of solutions...");
   -- put_line("The list of intrinsic solutions : ");
   -- put(DoblDobl_Output,Length_Of(isols),Head_Of(isols).n,isols);
   -- put_line("The solutions in the ambient space : ");
   -- put(DoblDobl_Output,Length_Of(asols),Head_Of(asols).n,asols);
    Compare(n,esols,asols);
    Clear(isols); Clear(asols);
    put_line("Doing the project and expand using matrices...");
    isols := Project(esols,m);
    asols := Expand(isols,m);
   -- put_line("The list of intrinsic solutions : ");
   -- put(DoblDobl_Output,Length_Of(isols),Head_Of(isols).n,isols);
   -- put_line("The solutions in the ambient space : ");
   -- put(DoblDobl_Output,Length_Of(asols),Head_Of(asols).n,asols);
    Compare(n,esols,asols);
    Clear(isols); Clear(asols);
  end Test_Project_and_Expand;

  procedure Test_Project_and_Expand
              ( n : in integer32;
                esols : in QuadDobl_Complex_Solutions.Solution_List; 
                m : in QuadDobl_Complex_Matrices.Matrix;
                b : in QuadDobl_Complex_Vectors.Vector;
                v : in QuadDobl_Complex_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   The projection and expansion of solution lists is tested
  --   for linear spaces in vector and matrix format,
  --   in quad double precision.

    use QuadDobl_Complex_Solutions;

    isols : Solution_List := Project(esols,b,v);
    asols : Solution_List := Expand(isols,b,v);

  begin
    put_line("Testing project and expand of solutions...");
   -- put_line("The list of intrinsic solutions : ");
   -- put(QuadDobl_Output,Length_Of(isols),Head_Of(isols).n,isols);
   -- put_line("The solutions in the ambient space : ");
   -- put(QuadDobl_Output,Length_Of(asols),Head_Of(asols).n,asols);
    Compare(n,esols,asols);
    Clear(isols); Clear(asols);
    put_line("Doing the project and expand using matrices...");
    isols := Project(esols,m);
    asols := Expand(isols,m);
   -- put_line("The list of intrinsic solutions : ");
   -- put(QuadDobl_Output,Length_Of(isols),Head_Of(isols).n,isols);
   -- put_line("The solutions in the ambient space : ");
   -- put(QuadDobl_Output,Length_Of(asols),Head_Of(asols).n,asols);
    Compare(n,esols,asols);
    Clear(isols); Clear(asols);
  end Test_Project_and_Expand;

  procedure Test_Sampling
              ( n,d,k : in integer32;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                esols : in Standard_Complex_Solutions.Solution_List ) is

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    use Standard_Complex_Vectors;
    use Standard_Complex_VecVecs;
    use Standard_Complex_Matrices;
    use Standard_Complex_Poly_Systems;

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

  procedure Test_Sampling
              ( n,d,k : in integer32;
                ep : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                esols : in DoblDobl_Complex_Solutions.Solution_List ) is

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_VecVecs;
    use DoblDobl_Complex_Matrices;
    use DoblDobl_Complex_Poly_Systems;

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

  procedure Test_Sampling
              ( n,d,k : in integer32;
                ep : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                esols : in QuadDobl_Complex_Solutions.Solution_List ) is

  -- ON ENTRY :
  --   n        number of variables before embedding;
  --   d        dimension of the solution component;
  --   k        co-dimension of the solution component;
  --   ep       embedded polynomial system.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_VecVecs;
    use QuadDobl_Complex_Matrices;
    use QuadDobl_Complex_Poly_Systems;

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

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Reads a witness set in standard double precision
  --   and tests the conversion to intrinsic coordinates.

    ep : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
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
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --    Reads a witness set in double double precision
  --    and tests the conversion to intrinsic coordinates.

    ep : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : DoblDobl_Complex_Solutions.Solution_List;
    n,d,k : integer32 := 0;

  begin
    DoblDobl_Read_Embedding(ep,sols,natural32(d));
    n := ep'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    Test_Sampling(n,d,k,ep.all,sols);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Reads a witness set in quad double precision
  --   and tests the conversion to intrisic coordinates.

    ep : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : QuadDobl_Complex_Solutions.Solution_List;
    n,d,k : integer32 := 0;

  begin
    QuadDobl_Read_Embedding(ep,sols,natural32(d));
    n := ep'last;
    new_line;
    put("Dimension of the solution component : "); put(d,1); new_line;
    put("Ambient dimension after embedding   : "); put(n,1); new_line;
    put("Ambient dimension before embedding  : "); put(n-d,1); new_line;
    n := n-d;
    k := n-d;
    put("Co-dimension of the component : "); put(k,1); new_line;
    Test_Sampling(n,d,k,ep.all,sols);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision and calls
  --   then the corresponding test routine.

    ans : character;

  begin
    new_line;
    put_line("Moving to intrinsic coordinates...");
    new_line;
    put_line("MENU to select precision level : ");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the level of precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_coords;
