with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Characters_and_Numbers;             use Characters_and_Numbers;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Functions;
with Standard_Random_Polynomials;        use Standard_Random_Polynomials;
with DoblDobl_Random_Polynomials;        use DoblDobl_Random_Polynomials;
with QuadDobl_Random_Polynomials;        use QuadDobl_Random_Polynomials;
with Sets_of_Unknowns_io;                use Sets_of_Unknowns_io;
with Partitions_of_Sets_of_Unknowns;     use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io;  use Partitions_of_Sets_of_Unknowns_io;
with Set_Structure,Set_Structure_io;
with Random_Product_Start_Systems;       use Random_Product_Start_Systems;
with Standard_Lined_Hypersurfaces;
with DoblDobl_Lined_Hypersurfaces;
with QuadDobl_Lined_Hypersurfaces;
with Hypersurface_Samplers;
with DoblDobl_Gridded_Hypersurfaces;
with QuadDobl_Gridded_Hypersurfaces;

procedure ts_hypsam is

-- DESCRIPTION :
--   Interactive environment to sample points on hypersurfaces.

  procedure Test_Sampler 
              ( n : in integer32; p : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Samples generic points from the surface defined by p,
  --   using standard double precision arithmetic.

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Functions;
    use Standard_Lined_Hypersurfaces;

    d : constant integer32 := Degree(p);
    ep : constant Eval_Poly := Create(p);
    b : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    v : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    t : Standard_Complex_Vectors.Vector(1..d);
    m : Standard_Natural_Vectors.Vector(1..d);
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 10*natural32(d);
    fail : boolean;

  begin
    for i in b'range loop
      b(i) := b(i)/2.0;
      v(i) := v(i)/2.0;
    end loop;
    put_line("Computing generic points on a random line ...");
    Generic_Points(Standard_Output,p,ep,natural32(d),b,v,eps,maxit,t,fail,m);
    if fail
     then put_line("Calculation of generic points failed.");
     else put_line("Calculation of generic points succeeded.");
    end if;
    put("The multiplicities of the roots : "); put(m); new_line;
  end Test_Sampler;

  procedure Test_Sampler 
              ( n : in integer32; p : in DoblDobl_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Samples generic points from the surface defined by p,
  --   using double double precision arithmetic.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Functions;
    use DoblDobl_Lined_Hypersurfaces;

    d : constant integer32 := Degree(p);
    ep : constant Eval_Poly := Create(p);
    b : DoblDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    v : DoblDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    t : DoblDobl_Complex_Vectors.Vector(1..d);
    m : Standard_Natural_Vectors.Vector(1..d);
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 10*natural32(d);
    fail : boolean;
    two : constant double_double := create(2.0);

  begin
    for i in b'range loop
      b(i) := b(i)/two;
      v(i) := v(i)/two;
    end loop;
    put_line("Computing generic points on a random line ...");
    Generic_Points(Standard_Output,p,ep,natural32(d),b,v,eps,maxit,t,fail,m);
    if fail
     then put_line("Calculation of generic points failed.");
     else put_line("Calculation of generic points succeeded.");
    end if;
    put("The multiplicities of the roots : "); put(m); new_line;
  end Test_Sampler;

  procedure Test_Sampler 
              ( n : in integer32; p : in QuadDobl_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Samples generic points from the surface defined by p,
  --   using double double precision arithmetic.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Functions;
    use QuadDobl_Lined_Hypersurfaces;

    d : constant integer32 := Degree(p);
    ep : constant Eval_Poly := Create(p);
    b : QuadDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    v : QuadDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    t : QuadDobl_Complex_Vectors.Vector(1..d);
    m : Standard_Natural_Vectors.Vector(1..d);
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 10*natural32(d);
    fail : boolean;
    two : constant quad_double := create(2.0);

  begin
    for i in b'range loop
      b(i) := b(i)/two;
      v(i) := v(i)/two;
    end loop;
    put_line("Computing generic points on a random line ...");
    Generic_Points(Standard_Output,p,ep,natural32(d),b,v,eps,maxit,t,fail,m);
    if fail
     then put_line("Calculation of generic points failed.");
     else put_line("Calculation of generic points succeeded.");
    end if;
    put("The multiplicities of the roots : "); put(m); new_line;
  end Test_Sampler;

  procedure Test_Multihomogeneous_Sampler
              ( n,m : in integer32;
                p : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Samples generic points of an m-homogeneous polynomial p in n variables.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Functions;
    use Hypersurface_Samplers;

    z : Partition(1..natural32(m));
    b : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    v,t : Standard_Complex_VecVecs.VecVec(1..m);
    dp : constant integer32 := Degree(p);
    ep : constant Eval_Poly := Create(p);
    eps : constant double_float := 1.0E-13;
    maxit : constant natural32 := 10*natural32(dp);
    fail : boolean;

  begin
    Create(z,natural32(n));
    put_line("Reading the partition of the set of variables...");
    put_line("Every set in the partition is enclosed by curly braces.");
    for i in 1..m loop
      put("Give set "); put(i,1); put(" : "); get(z(natural32(i)));
    end loop;
    put("The partition : "); put(z); new_line;
    v := Random_Multihomogeneous_Directions(natural32(n),z);
    Generic_Points(Standard_Output,p,ep,z,b,v,eps,maxit,t,fail);
    put("Calculation of generic points for all sets");
    if fail
     then put_line(" failed.");
     else put_line(" succeeded.");
    end if;
  end Test_Multihomogeneous_Sampler;

  procedure Test_Set_Structure_Sampler 
              ( n : in integer32;
                p : in Standard_Complex_Polynomials.Poly ) is

  -- DESCRIPTION :
  --   Computes witness points on the surface defined by p(x) = 0,
  --   in n-space, exploiting its product structure.

    ns : Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);
    d : constant integer32 := Standard_Complex_Polynomials.Degree(p);

  begin
    ns(1) := natural32(d);
    Set_Structure.Init(ns);
    Build_Set_Structure(1,natural32(d),p);
    put("The set structure : ");
    Set_Structure_io.put(1); new_line;
  end Test_Set_Structure_Sampler;

  procedure Add_Default_Symbols ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Adds n default symbols to denote the variables in a polynomial,
  --   as needed when working with random polynomials.

    use Symbol_Table;

  begin
    for i in 1..n loop
      declare
        sb : Symbol;
        sn : constant string := Convert(i);
      begin
        sb := (sb'range => ' ');
        sb(sb'first) := 'x';
        for j in sn'range loop
          sb(sb'first+j) := sn(j);
        end loop;
        Symbol_Table.Add(sb);
      end;
    end loop;
  end Add_Default_Symbols;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Tests the sampling in standard double precision.

    use Standard_Complex_Polynomials;

    ans : character;
    n,d,m : natural32 := 0;
    p : Poly;
 
  begin
    new_line;
    put_line("Testing sampling of points on a hypersurface.");
    new_line;
    put_line("Choose one of the following options : ");
    put_line("  1. sample a randomly generated complex polynomial; or");
    put_line("  2. give your own polynomial to sample.");
    put("Type 1 or 2 to choose : "); 
    Ask_Alternative(ans,"12");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    case ans is
      when '1' =>
        Add_Default_Symbols(integer32(n));
        put("Give the degree bound : "); get(d);
        put("Give the number of terms : "); get(m);
        p := Random_Sparse_Poly(n,d,m,0);
        put_line("A random polynomial : ");
        put_line(p);
      when '2' =>
        put("Give your polynomial : "); get(p);
        put("-> Your polynomial : "); put(p); new_line;
      when others => null;
    end case;
    new_line;
    put_line("MENU for exploiting structure :");
    put_line("  0. exploit no structure, just plain degree;");
    put_line("  1. exploit given multi-homogeneous structure;");
    put_line("  2. let automatically generate a set structure.");
    put("Type 0, 1, or 2 to make your choice : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Test_Sampler(integer32(n),p);
      when '1' =>
        put("Give number of sets in the partition : "); 
        get(m);
        Test_Multihomogeneous_Sampler(integer32(n),integer32(m),p);
      when '2' => 
        Test_Set_Structure_Sampler(integer32(n),p);
      when others => null;
    end case;
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Tests the sampling in double double precision.

    use DoblDobl_Complex_Polynomials;

    ans : character;
    n,d,m : natural32 := 0;
    p : Poly;
 
  begin
    new_line;
    put_line("Testing sampling of points on a hypersurface.");
    new_line;
    put_line("Choose one of the following options : ");
    put_line("  1. sample a randomly generated complex polynomial; or");
    put_line("  2. give your own polynomial to sample.");
    put("Type 1 or 2 to choose : "); 
    Ask_Alternative(ans,"12");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    case ans is
      when '1' =>
        Add_Default_Symbols(integer32(n));
        put("Give the degree bound : "); get(d);
        put("Give the number of terms : "); get(m);
        p := Random_Sparse_Poly(n,d,m,0);
        put_line("A random polynomial : ");
        put_line(p);
      when '2' =>
        put("Give your polynomial : "); get(p);
        put("-> Your polynomial : "); put(p); new_line;
      when others => null;
    end case;
    Test_Sampler(integer32(n),p);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Tests the sampling in quad double precision.

    use QuadDobl_Complex_Polynomials;

    ans : character;
    n,d,m : natural32 := 0;
    p : Poly;
 
  begin
    new_line;
    put_line("Testing sampling of points on a hypersurface.");
    new_line;
    put_line("Choose one of the following options : ");
    put_line("  1. sample a randomly generated complex polynomial; or");
    put_line("  2. give your own polynomial to sample.");
    put("Type 1 or 2 to choose : "); 
    Ask_Alternative(ans,"12");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    case ans is
      when '1' =>
        Add_Default_Symbols(integer32(n));
        put("Give the degree bound : "); get(d);
        put("Give the number of terms : "); get(m);
        p := Random_Sparse_Poly(n,d,m,0);
        put_line("A random polynomial : ");
        put_line(p);
      when '2' =>
        put("Give your polynomial : "); get(p);
        put("-> Your polynomial : "); put(p); new_line;
      when others => null;
    end case;
    Test_Sampler(integer32(n),p);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the level of precision
  --   and then calls the appropriate test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the working precision :");
    put_line("  0. standard double precision; or");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
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
end ts_hypsam;
