with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;          use Standard_Complex_Vectors;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;
--with Random_Complex_Polynomials;        use Random_Complex_Polynomials;
with Polynomial_Roots;                  use Polynomial_Roots;
with Hypersurface_Roots;                use Hypersurface_Roots;
with Monodromy_Partitions;              use Monodromy_Partitions;

procedure ts_hypfac is

-- DESCRIPTION :
--   Interactive development of factorization of multivariate polynomials,
--   or equivalently, the decomposition of a hypersurface into irreducible
--   components.

  tol : constant double_float := 1.0E-8;

  procedure Random_Irreducible_Poly ( n : out natural32; p : out Poly ) is

  -- DESCRIPTION :
  --   Interactive generation of a random irreducible polynomial.

    d,m : natural32 := 0;

  begin
    n := 0;
    put_line("Generating a random multivariate complex polynomial...");
    put("  Give the number of variables : "); get(n);
    put("  Give the degree : "); get(d);
    put("  Give the number of terms : "); get(m);
   -- p := Random(n,d,m);
    p := Random_Sparse_Poly(n,d,m,0);
    put_line("The random polynomial : ");
    put_line(p);
  end Random_Irreducible_Poly;

  procedure Random_Reducible_Poly ( n,m : out natural32; p : out Poly ) is

  -- DESCRIPTION :
  --   Generates a random reducible polynomial in n variables consisting
  --   m irreducible factors.  The degrees of the factors are determined
  --   interactively.

    f : Poly;
    d,k,mu : natural32 := 0; 

  begin
    n := 0; m := 0;
    put_line("Generating a random multivariate complex polynomial...");
    put("  Give the number of variables : "); get(n);
    put("  Give the number of factors : "); get(m);
    for i in 1..m loop
      put("Generating factor "); put(i,1); put_line(" :");
      put("  Give its degree : "); get(d);
      put("  Give its number of terms : "); get(k);
      put("  Give its multiplicity : "); get(mu);
     -- f := Random(n,d,k);
      f := Random_Sparse_Poly(n,d,k,0);
      if i = 1
       then Copy(f,p);
       else Mul(p,f);
      end if;
      for j in 1..mu-1 loop
        Mul(p,f);
      end loop;
      Clear(f);
    end loop;
  end Random_Reducible_Poly;

  procedure Build_Affine_Loop
                ( deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                  nbdc : in out natural32; n : in natural32; p : in Poly;
                  cp,v0,sA : in Vector; fail : out boolean ) is

  -- DESCRIPTION :
  --   Makes one loop v0 -> v1 -> v2 -> v0 and updates the decomposition.
  --   Failure is reported when the paths failed to converge.
  --   One retry is allowed.

    v1,v2 : Vector(1..integer32(n));
    sB : Vector(sA'range);
    mp : Standard_Natural_Vectors.Vector(sA'range);
    res : double_float;

  begin
    for i in 1..2 loop
      v1 := Random_Vector(1,integer32(n));
      v2 := Random_Vector(1,integer32(n));
      sB := sA;
      Affine_Track_Moving_Line(p,v0,v1,sB);
      Affine_Track_Moving_Line(p,v1,v2,sB);
      Affine_Track_Moving_Line(p,v2,v0,sB);
      Test_Affine_Roots(Standard_Output,cp,sB,res);
      fail := (res > tol);
      if not fail then
        mp := Map(sA,sB,tol);
        Write_Map(Standard_Output,mp);
        Add_Map(deco,nbdc,mp);
        Write_Factors(Standard_Output,deco.all);
      end if;
      exit when not fail;
    end loop;
  end Build_Affine_Loop;

  procedure Build_Projective_Loop
                ( deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                  nbdc : in out natural32; n : in natural32; p : in Poly;
                  cp,v0,s0A,s1A : in Vector; fail : out boolean ) is

  -- DESCRIPTION :
  --   Makes one loop v0 -> v1 -> v2 -> v0 and updates the decomposition.
  --   Failure is reported when the paths failed to converge.

    v1,v2 : Vector(1..integer32(n));
    s0B,s1B,sA,sB : Vector(s0A'range);
    mp : Standard_Natural_Vectors.Vector(s0A'range);
    res : double_float;

  begin
    for i in 1..2 loop
      v1 := Random_Vector(1,integer32(n));
      v2 := Random_Vector(1,integer32(n));
      s0B := s0A; s1B := s1A;
      Projective_Track_Moving_Line(p,v0,v1,s0B,s1B);
      Projective_Track_Moving_Line(p,v1,v2,s0B,s1B);
      Projective_Track_Moving_Line(p,v2,v0,s0B,s1B);
      Test_Projective_Roots(Standard_Output,cp,s0B,s1B,res);
      fail := (res > tol);
      if not fail then
        for i in sA'range loop
          sA(i) := s1A(i)/s0A(i);
          sB(i) := s1B(i)/s0B(i);
        end loop;
        mp := Map(sA,sB,tol);
        Write_Map(Standard_Output,mp);
        Add_Map(deco,nbdc,mp);
        Write_Factors(Standard_Output,deco.all);
      end if;
      exit when not fail;
    end loop;
  end Build_Projective_Loop;

  procedure Test_Affine_Monodromy_Loop
                ( n : in natural32; p : in Poly; threshold : in natural32;
                  fail : out boolean ) is

  -- DESCRIPTION :
  --   Tests the monodromy breakup algorithm to factor p, stops when
  --   either the polynomial is found to be irreducible, 
  --   or there is no change in the partition after #threshold loops,
  --   or numerical failure prevents the formation of loops.

    v : constant Vector(1..integer32(n)) := Random_Vector(1,integer32(n));
    c : constant Vector(0..Degree(p)) := Substitute(p,v);
    s : Vector(1..c'last);
    nbfac : natural32 := natural32(c'last);
    prev_nbfac,cnt : natural32;
    res : double_float;
    deco : Standard_Natural_VecVecs.Link_to_VecVec := Init_Factors(nbfac);

  begin
    for i in 1..2 loop
      Affine_Solve(Standard_Output,c,s,res);
      fail := (res > tol);
      exit when not fail;
    end loop;
    prev_nbfac := nbfac;
    cnt := 0;
    while ((nbfac > 1) and (cnt < threshold) and (not fail)) loop
      Build_Affine_Loop(deco,nbfac,n,p,c,v,s,fail);
      nbfac := Number_of_Factors(deco.all);
      if prev_nbfac = nbfac then
        cnt := cnt + 1;
      else
        cnt := 0;
        prev_nbfac := nbfac;
      end if;
    end loop;
  end Test_Affine_Monodromy_Loop;

  procedure Test_Projective_Monodromy_Loop
                ( n : in natural32; p : in Poly; threshold : in natural32;
                  fail : out boolean ) is

  -- DESCRIPTION :
  --   Tests the monodromy breakup algorithm in projective coordinates.

    v : constant Vector(1..integer32(n)) := Random_Vector(1,integer32(n));
    c : constant Vector(0..Degree(p)) := Substitute(p,v);
    s0,s1 : Vector(1..c'last);
    nbfac : natural32 := natural32(c'last);
    prev_nbfac,cnt : natural32;
    res : double_float;
    deco : Standard_Natural_VecVecs.Link_to_VecVec := Init_Factors(nbfac);

  begin
    for i in 1..2 loop
      Projective_Solve(Standard_Output,c,s0,s1,res);
      fail := (res > tol);
      exit when not fail;
    end loop;
    prev_nbfac := nbfac;
    cnt := 0;
    while ((nbfac > 1) and (cnt < threshold) and (not fail)) loop
      Build_Projective_Loop(deco,nbfac,n,p,c,v,s0,s1,fail);
      nbfac := Number_of_Factors(deco.all);
      if prev_nbfac = nbfac then
        cnt := cnt + 1;
      else
        cnt := 0;
        prev_nbfac := nbfac;
      end if;
    end loop;
  end Test_Projective_Monodromy_Loop;

  procedure Random_Tester ( nb : in natural32; ap : character ) is

  -- DESCRIPTION :
  --   Performs nb random tests on the monodromy breakup algorithm.

    n,d,m : natural32 := 0;
    p : Poly;
    fail : boolean;

  begin
    put_line("Generation of random polynomials...");
    put("  Give the number of variables : "); get(n);
    put("  Give the degree : "); get(d);
    put("  Give the number of terms : "); get(m);
    for i in 1..nb loop
     -- p := Random(n,d,m);
      p := Random_Sparse_Poly(n,d,m,0);
      put_line(p);
      if ap = 'a'
       then Test_Affine_Monodromy_Loop(n,p,10,fail);
       else Test_Projective_Monodromy_Loop(n,p,10,fail);
      end if;
      exit when fail;
    end loop;
    if fail
     then put_line("Ended in failure.");
     else put("Tested "); put(nb,1); put_line(" cases successfully.");
    end if;
  end Random_Tester;

  procedure Factor_Given_Polynomial is

    ans : character;
    file : file_type;
    n : natural32 := 0;
    p : Poly;
    fail : boolean;

  begin
    put("Is the polynomial on file? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      new_line;
      put_line("Reading the name of the file for the polynomial...");
      Read_Name_and_Open_File(file);
      get(file,n);
      Symbol_Table.Init(n);
      get(file,p);
    else
      new_line;
      put("Give the number of variables : ");
      get(n);
      Symbol_Table.Init(n);
      put_line("Give your polynomial : ");
      get(p);
    end if;
    new_line;
    put_line("Your polynomial is "); put(p);
    new_line;
    put_line("Starting the monodromy algorithm...");
    Test_Affine_Monodromy_Loop(n,p,10,fail);
  end Factor_Given_Polynomial;

  procedure Factor_Random_Polynomial is

    n,m : natural32 := 0;
    p : Poly;
    fail : boolean;

  begin
    Random_Reducible_Poly(n,m,p);
    new_line;
    put_line("The random polynomial is ");
    put_line(p);
    new_line;
    put_line("Starting the monodromy algorithm...");
    Test_Affine_Monodromy_Loop(n,p,10,fail);
  end Factor_Random_Polynomial;  

  procedure Main is

    n : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Decomposing a hypersurface into irreducible components.");
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. Perform sequence of random tests.");
    put_line("  2. Factor a given polynomial.");
    put_line("  3. Factor a random reducible polynomial.");
    put("Type 1, 2, or 3 to select your choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    case ans is
      when '1' => put("Give the number of random tests : "); get(n);
                  put("Affine or projective coordinates ? (a/p) ");
                  Ask_Alternative(ans,"ap");
                  Random_Tester(n,ans);
      when '2' => Factor_Given_Polynomial;
      when '3' => Factor_Random_Polynomial;
      when others => put_line("Invalid option.");
    end case;
  end Main;

begin
  Main;
end ts_hypfac;
