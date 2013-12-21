with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer64_Matrices;
with Standard_Lattice_Polytopes;
with Multprec_Integer_Matrices;
with Multprec_Lattice_Polytopes;
with Convex_Hull_Methods;                use Convex_Hull_Methods;

procedure ts_convhull is

-- DESCRIPTION :
--   Test on convex hull algorithms.

  procedure Standard_Convex_Hull
              ( A : in Standard_Integer64_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes the convex hull of the set spanned by the points
  --   in the columns of A.

    n : constant integer32 := A'last(1);

  begin
    case n is
      when 2 => Standard_Planar_Hull(A);
      when 3 => Standard_3D_Hull(A);
      when 4 => Standard_4D_Hull(A);
      when others => Standard_General_Hull(A);
    end case;
  end Standard_Convex_Hull;

  procedure Multprec_Convex_Hull
              ( A : in Multprec_Integer_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Computes the convex hull of the set spanned by the points
  --   in the columns of A.

    n : constant integer32 := A'last(1);

  begin
    case n is
      when 2 => Multprec_Planar_Hull(A);
      when 3 => Multprec_3D_Hull(A);
      when 4 => Multprec_4D_Hull(A);
      when others => Multprec_General_Hull(A);
    end case;
  end Multprec_Convex_Hull;

  procedure Standard_Test_Convex_Hull ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Test on computing the convex hull of a lattice polygon
  --   or 3D polytope, using standard 64-bit arithmetic.

    m : integer32 := 0;
    ans : character;

  begin
    loop
      put("Give the number of points : "); get(m);
      declare
        A : Standard_Integer64_Matrices.Matrix(1..n,1..m);
      begin
        put("-> generate random coordinates ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          A := Random_Data(n,m);
        else
          put("-> work with a cyclic polytope ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'y'
           then A := Cyclic_Polytope(n,m);
           else A := User_Input(n,m);
          end if;
        end if;
        declare
          B : constant Standard_Integer64_Matrices.Matrix
            := Filter_Duplicates(A);
          r : constant natural32 := Standard_Lattice_Polytopes.Rank(B);
        begin
          put("The rank of the support : "); put(r,1); 
          if integer32(r) = n then
            new_line;
            Standard_Convex_Hull(B);
          else
            put("  Continue ? (y/n) ");
            Ask_Yes_or_No(ans);
            if ans = 'y'
             then Standard_Convex_Hull(B);
            end if;
          end if;
        end;
      end;
      put("Test another case ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Standard_Test_Convex_Hull;

  procedure Multprec_Test_Convex_Hull ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Computes the convex hull of a lattice polygon or 3D polytope,
  --   using exact multiprecision integer arithmetic.

    m : integer32 := 0;
    ans : character;

  begin
    loop
      put("Give the number of points : "); get(m);
      declare
        A : Multprec_Integer_Matrices.Matrix(1..n,1..m);
      begin
        put("-> generate random coordinates ? (y/n) ");
        Ask_Yes_or_No(ans);
        if ans = 'y' then
          A := Random_Data(n,m);
        else
          put("-> work with a cyclic polytope ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'y' 
           then A := Cyclic_Polytope(n,m);
           else A := User_Input(n,m);
          end if;
        end if;
        declare
          B : constant Multprec_Integer_Matrices.Matrix
            := Filter_Duplicates(A);
          r : constant natural32 := Multprec_Lattice_Polytopes.Rank(B);
        begin
          put("The rank of the support : "); put(r,1);
          if integer32(r) = n then
            new_line;
            Multprec_Convex_Hull(B);
          else
            put("  Continue ? (y/n) ");
            Ask_Yes_or_No(ans);
            if ans = 'y'
             then Multprec_Convex_Hull(B);
            end if;
          end if;
        end;
      end;
      put("Test another case ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Multprec_Test_Convex_Hull;

  procedure Main is

    ans : character;
    n : integer32 := 0;

  begin
    new_line;
    put_line("MENU to test convex hull computations :");
    put_line("  1. planar convex hull in standard arithmetic;");
    put_line("  2. use multiprecision to compute the planar convex hull");
    put_line("  3. convex hull in 3D using standard arithmetic;");
    put_line("  4. use multiprecision for 3D convex hulls;");
    put_line("  5. convex hull in 4D using standard arithmetic;");
    put_line("  6. use multiprecision for 4D convex hulls;");
    put_line("  7. convex hull in any dimension with standard arithmetic;");
    put_line("  8. use multiprecision for any dimensional convex hulls.");
    put("Type 1, 2, 3, 4, 5, 6, 7, or 8 to select test : ");
    Ask_Alternative(ans,"12345678");
    new_line;
    if ans = '7' or ans = '8'
     then put("Give the dimension : "); get(n);
    end if;
    case ans is
      when '1' => Standard_Test_Convex_Hull(2);
      when '2' => Multprec_Test_Convex_Hull(2);
      when '3' => Standard_Test_Convex_Hull(3);
      when '4' => Multprec_Test_Convex_Hull(3);
      when '5' => Standard_Test_Convex_Hull(4);
      when '6' => Multprec_Test_Convex_Hull(4);
      when '7' => Standard_Test_Convex_Hull(n);
      when '8' => Multprec_Test_Convex_Hull(n);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_convhull;
