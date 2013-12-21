with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Integer64_Vectors;
with Standard_Integer64_Vectors_io;      use Standard_Integer64_Vectors_io;
with Standard_Integer64_Matrices;        use Standard_Integer64_Matrices;
with Standard_Integer64_Matrices_io;     use Standard_Integer64_Matrices_io;
with Standard_Integer64_VecMats;         use Standard_Integer64_VecMats;
with Standard_Integer64_VecMats_io;      use Standard_Integer64_VecMats_io;
with Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Lists_of_Integer64_Vectors;
with Lists_of_Integer64_Vectors_io;      use Lists_of_Integer64_Vectors_io;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laurentials_io;    use Standard_Complex_Laurentials_io;
with Standard_Random_Laurentials;        use Standard_Random_Laurentials;
with Symbol_Table;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_Systems_io;   use Standard_Complex_Laur_Systems_io;
with Supports_of_Polynomial_Systems;
with Standard_Lattice_Supports;
with Standard_Lattice_3d_Facets;         use Standard_Lattice_3d_Facets;
with Standard_Lattice_3d_Facets_io;      use Standard_Lattice_3d_Facets_io;
with Standard_Pretropical_Facets;        use Standard_Pretropical_Facets;

procedure ts_pretrop is

-- DESCRIPTION :
--   Interactive development of the computation of pretropisms.

  procedure Write_Facet_Tuples ( f : in Array_of_Facet_3d_Lists ) is

  -- DESCRIPTION :
  --   Given in f an array of three equally long lists,
  --   the tuples are written.

    tmp : Array_of_Facet_3d_Lists(f'range) := f;
    lft : Link_to_3d_Facet;

  begin
    for i in 1..Length_Of(f(f'first)) loop
      put("tuple "); put(i,1); put(" with normal");
      lft := Head_Of(tmp(tmp'first));
      put(lft.normal);
      put(" spanned by (");
      put(lft.points);
      for j in tmp'first+1..tmp'last loop
        lft := Head_Of(tmp(j));
        put(","); put(lft.points);
      end loop;
      put_line(")");
      for j in tmp'range loop
        tmp(j) := Tail_Of(tmp(j));
      end loop;
    end loop;
  end Write_Facet_Tuples;

  procedure Compute_Pretropisms
               ( A : in VecMat; v : out Lists_of_Integer64_Vectors.List ) is

  -- DESCRIPTION :
  --   Computes facet and edge pretropisms.
  --   Facet pretropisms are returned in v.

    f,g : Array_of_Facet_3d_Lists(A'range);

  begin
    put_line("The supports of the polynomial system : "); put(A,3);
    for i in A'range loop
      f(i) := Convex_Hull_3D(A(i).all);
      put("checking Euler characteristic for support "); put(i,1);
      put_line(" :");
      Check_Euler_Characteristic(A(i)'last(2),f(i));
    end loop;
    v := Facet_Pretropisms(f);
    put_line("The list of facet pretropisms : "); put(v);
    g := Pretropical_Facets(f,v);
    put_line("The list of pretropical facets : ");
    Write_Facet_Tuples(g);
    Edge_Pretropisms(A,f,g,v);
  end Compute_Pretropisms;

  procedure Compute_Pretropisms ( p : in Poly_Sys ) is

    s : constant Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);
    A : constant VecMat(p'range) := Lists2VecMat(s);
    v : Lists_of_Integer64_Vectors.List;

  begin
    Compute_Pretropisms(A,v);
    put_line("The list of facet pretropisms : "); put(v);
  end Compute_Pretropisms;

  procedure Compute_Pretropisms ( p : in Laur_Sys ) is

    s : constant Array_of_Lists(p'range)
      := Supports_of_Polynomial_Systems.Create(p);
    A : constant VecMat(p'range) := Lists2VecMat(s);
    v : Lists_of_Integer64_Vectors.List;

  begin
    Compute_Pretropisms(A,v);
    put_line("The list of facet pretropisms : "); put(v);
  end Compute_Pretropisms;

  procedure Test_Pretropisms is

  -- DESCRIPTION :
  --   Reads a polynomial system and calls the pretropism calculator.

    lp : Link_to_Poly_Sys;

  begin
    put_line("Computing pretropisms of Newton polytopes");
    new_line;
    put_line("Reading a polynomial system ..."); get(lp);
    Compute_Pretropisms(lp.all);
  end Test_Pretropisms;

  procedure Verify_Edge_Tropisms
               ( A,B : in Matrix;
                 fpt,p,q : in Lists_of_Integer64_Vectors.List ) is

  -- DESCRIPTION :
  --   Verifies whether all vectors of p are valid tropisms and
  --   also occur in the list q.  The list fpt contains the normals
  --   to the pretropical facets.

    bug : boolean := false;
    tmp : Lists_of_Integer64_Vectors.List := p;
    ltv : Standard_Integer64_Vectors.Link_to_Vector;
    sup : Standard_Integer_VecVecs.VecVec(1..2);

  begin
    while not Lists_of_Integer64_Vectors.Is_Null(tmp) loop
      ltv := Lists_of_Integer64_Vectors.Head_Of(tmp);
      put(ltv.all); put(" supports ");
      sup := Standard_Lattice_Supports.Support(A,B,ltv.all);
      put("("); put(sup(1)); put(","); put(sup(2)); put(")");
      if Is_Tropism(A,B,ltv.all)
       then put(" is tropism");
       else put(" is NOT a tropism"); bug := true;
      end if;
     -- exit when bug;
      if Lists_of_Integer64_Vectors.Is_In(q,ltv.all) then
        put_line(" occurs in the other list, OK");
      elsif On_Pretropical_Facet(A,B,fpt,ltv.all) then
        put_line(" lies on pretropical facet ...");
      else
        put_line(" does NOT occur in the other list, bug!"); bug := true;
      end if;
     -- exit when bug;
      Standard_Integer_VecVecs.Clear(sup);
      tmp := Lists_of_Integer64_Vectors.Tail_Of(tmp);
    end loop;
    if bug
     then put_line("A bug occurred!");
     else put_line("No bugs to report.");
    end if;
  end Verify_Edge_Tropisms;

  procedure Show_Edge_Tropisms
               ( p,q : in Poly; t : out Lists_of_Integer64_Vectors.List ) is 

  -- DESCRIPTION :
  --   For two Laurent polynomials p and q, returns in t the tropisms
  --   to the edges of the Newton polytopes of p and q.

    sp : Lists_of_Integer_Vectors.List
       := Supports_of_Polynomial_Systems.Create(p);
    sq : Lists_of_Integer_Vectors.List
       := Supports_of_Polynomial_Systems.Create(q);
    A : constant Matrix := List2Matrix(sp);
    B : constant Matrix := List2Matrix(sq);
    ans : character;
    fpt,w : Lists_of_Integer64_Vectors.List;

  begin
    put_line("first support matrix : "); put(A);
    put_line("second support matrix : "); put(B);
    new_line;
    put("Assume supports in general position ? (y/n) ");
    Ask_Yes_or_No(ans);
    t := Edge_Tropisms_by_Sum(A,B);
    if ans = 'y'
     then w := Edge_Tropisms_by_Wrap(A,B);
     else Edge_Tropisms_by_Wrap(A,B,fpt,w);
    end if;
    Lists_of_Integer_Vectors.Clear(sp);
    Lists_of_Integer_Vectors.Clear(sq);
    put_line("The edge tropisms via Minkowski sum : "); put(t);
    if not Lists_of_Integer64_Vectors.Is_Null(fpt) 
     then put_line("The facet pretropisms : "); put(fpt);
    end if;
    if not Lists_of_Integer64_Vectors.Is_Null(w) 
     then put_line("The edge tropisms via giftwrapping : "); put(w);
    end if;
    put_line("verifying Minkowski sum results with giftwrapping ...");
    Verify_Edge_Tropisms(A,B,fpt,t,w);
    put_line("verifying giftwrapping results with Minkowski sum ...");
    Verify_Edge_Tropisms(A,B,fpt,w,t);
  end Show_Edge_Tropisms;

  procedure Show_Facets
               ( p : in Poly; g : out Lists_of_Integer64_Vectors.List ) is

  -- DESCRIPTION :
  --   Shows the facet normals of the convex hull spanned by the support
  --   of the polynomial p.  If p is the common hypersurface, then this
  --   display is needed to check if all pretropical facets have been found.
  --   Returns in g the list of inner facet normals.

    s : Lists_of_Integer_Vectors.List
      := Supports_of_Polynomial_Systems.Create(p);
    A : constant Matrix := List2Matrix(s);
    f : constant Facet_3d_List := Convex_Hull_3D(A);

  begin
    Write_Facets(A,f);
    Check_Euler_Characteristic(A'last(2),f);
    Lists_of_Integer_Vectors.Clear(s);
    g := Facet_Normals(f);
  end Show_Facets;

  function Generate_Random_Laurent_Polynomial return Poly is

  -- DESCRIPTION :
  --   Interactive generation of a random Laurent polynomial:
  --   prompts the user for the number of monomials, lower and upper
  --   bounds for the exponents.

    res : Poly;
    n : constant natural32 := 3;
    m : natural32 := 0;
    lower,upper : integer32 := 0;

  begin
    loop
      put("Give #monomials for a random Laurent polynomial : "); get(m);
      exit when (m > 0);
      put("zero monomials is not acceptable, enter at least 1 ...");
    end loop;
    put("Give a lower bound for exponents : "); get(lower);
    put("Give an upper bound for exponents : "); get(upper);
    res := Random_Laurent_Polynomial(n,m,lower,upper);
    put_line("-> a random Laurent polynomial : "); put_line(res); new_line;
    return res;
  end Generate_Random_Laurent_Polynomial;

  procedure Default_Random_System 
              ( h : out Poly; c : out Laur_Sys; f : out Laur_Sys ) is

  -- DESCRIPTION :
  --   Generates a random polynomial system with a common hypersurface,
  --   a common curve and random cofactors, with default parameters:
  --   4 random monomials with exponents ranging between -9 and +9.

  -- ON RETURN :
  --   h        a random hypersurface with 4 monomials in 3 variables,
  --            with exponents ranging between -9 and +9;
  --   c        two random polynomials as h defining a curve;
  --   f        three random polynomials, all have h as factor,
  --            the first two have c(1) and c(2) as factor and
  --            the last one has c(1)*c(2) as factor.

    n : constant natural32 := 3;
    m : constant natural32 := 4;
    lower : constant integer32 := -9;
    upper : constant integer32 := +9;

  begin
    h := Random_Laurent_Polynomial(n,m,lower,upper);
    c(1) := Random_Laurent_Polynomial(n,m,lower,upper);
    c(2) := Random_Laurent_Polynomial(n,m,lower,upper);
    f(1) := Random_Laurent_Polynomial(n,m,lower,upper);
    Mul(f(1),h); Mul(f(1),c(1));
    f(2) := Random_Laurent_Polynomial(n,m,lower,upper);
    Mul(f(2),h); Mul(f(2),c(2));
    f(3) := Random_Laurent_Polynomial(n,m,lower,upper);
    Mul(f(3),h); Mul(f(3),c(1)); Mul(f(3),c(2));
  end Default_Random_System;

  procedure Interactive_Random_System
              ( h : out Poly; c : out Laur_Sys; f : out Laur_Sys ) is

  -- DESCRIPTION :
  --   Generates a random polynomial system with a common hypersurface,
  --   a common curve and random cofactors, prompting the user each
  --   time for the number of monomials and range of exponents.

  -- ON RETURN :
  --   h        a random hypersurface in 3 variables,
  --   c        two random polynomials as h defining a curve;
  --   f        three random polynomials, all have h as factor,
  --            the first two have c(1) and c(2) as factor and
  --            the last one has c(1)*c(2) as factor.

  begin
    put_line("generating a random hypersurface ...");
    h := Generate_Random_Laurent_Polynomial;
    put_line("generating a random curve ...");
    for i in c'range loop
      c(i) := Generate_Random_Laurent_Polynomial;
    end loop;
    put_line("generating random cofactors ...");
    for i in f'range loop
      f(i) := Generate_Random_Laurent_Polynomial;
      Mul(f(i),h);
      if i <= c'last
       then Mul(f(i),c(i));
       else Mul(f(i),c(1)); Mul(f(i),c(2));
      end if;
    end loop;
  end Interactive_Random_System;

  procedure Generate_Random_System is

    f : Laur_Sys(1..3);
    h : Poly;
    c : Laur_Sys(1..2);
    ans : character;
    hfs,ctp : Lists_of_Integer64_Vectors.List;

  begin
    Symbol_Table.Init(3);
    Symbol_Table.Add_String("x");
    Symbol_Table.Add_String("y");
    Symbol_Table.Add_String("z");
    put("Random system with default parameters ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Default_Random_System(h,c,f);
     else Interactive_Random_System(h,c,f);
    end if;
    put("Do you want to see the random system ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("The random system : "); put_line(f);
    end if;
    put("See the facets of the common hypersurface factor ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Show_Facets(h,hfs);
    end if;
    put("See the edge tropisms of the common curve ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Show_Edge_Tropisms(c(1),c(2),ctp);
    end if;
    put("Continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Compute_Pretropisms(f);
    end if;
    put_line("Facet normals to the common hypersurface : "); put(hfs);
  end Generate_Random_System;

  function Degree ( t : Standard_Integer64_Vectors.Vector ) return integer64 is

  -- DESCRIPTION :
  --   The degree of a tropism determines the degree of the initial form.

     min : integer64 := t(t'first);
     max : integer64 := t(t'first);

  begin
    for i in t'first+1..t'last loop
      if t(i) < min then
        min := t(i);
      elsif t(i) > max then
        max := t(i);
      end if;
    end loop;
    if min >= 0 then
      return max;
    elsif max <= 0 then
      return -min;
    else
      return max - min;
    end if;
  end Degree;

  function Positive ( t : Standard_Integer64_Vectors.Vector ) return boolean is

  -- DESCRIPTION :
  --   A tropism is positive if its first nonzero entry is positive.

  begin
    for i in t'range loop
      if t(i) > 0 then
        return true;
      elsif t(i) < 0 then
        return false;
      end if;
    end loop;
    return false;
  end Positive;

  function Degree ( t : Lists_of_Integer64_Vectors.List ) return integer64 is

  -- DESCRIPTION :
  --   Considers all tropisms for which the first nonzero entry is positive
  --   and returns the sum of their degrees.

    res : integer64 := 0;
    tmp : Lists_of_Integer64_Vectors.List := t;
    ltv : Standard_Integer64_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer64_Vectors.Is_Null(tmp) loop
      ltv := Lists_of_Integer64_Vectors.Head_Of(tmp);
      if Positive(ltv.all)
       then res := res + Degree(ltv.all);
      end if;
      tmp := Lists_of_Integer64_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Degree;

  procedure Pretropisms_for_Space_Curves is

    ans : character;
    p,q : Poly;
    t : Lists_of_Integer64_Vectors.List;

  begin
    put_line("Pretropisms for two polynomials in three variables.");
    new_line;
    put("Generate random polynomials for the space curve ? (y/n) ");
    Ask_Yes_or_No(ans);
    Symbol_Table.Init(3);
    Symbol_Table.Add_String("x");
    Symbol_Table.Add_String("y");
    Symbol_Table.Add_String("z");
    if ans = 'y' then
      p := Generate_Random_Laurent_Polynomial;
      q := Generate_Random_Laurent_Polynomial;
    else
      put_line("Reading polynomial p in x, y, z, terminate with semicolon ...");
      put("Give p : "); get(p);
      put_line("Reading polynomial q in x, y, z, terminate with semicolon ...");
      put("Give q : "); get(q);
    end if;
    Show_Edge_Tropisms(p,q,t);
    put("Computed ");
    put(Lists_of_Integer64_Vectors.Length_Of(t),1);
    put_line(" tropisms");
    put("The degree of the space curve : "); put(Degree(t),1); new_line;
  end Pretropisms_for_Space_Curves;

  procedure Main is

    ans : character;
 
  begin
    new_line;
    put_line("MENU to test the computation of pretropisms :");
    put_line("  1. compute pretropisms for a space curve;");
    put_line("  2. provide a polynomial system for pretropisms;");
    put_line("  3. generate a random system to test pretropisms.");
    put("Type 1, 2, or 3 to make a choice : ");
    Ask_Alternative(ans,"123");
    new_line;
    case ans is
       when '1' => Pretropisms_for_Space_Curves;
       when '2' => Test_Pretropisms;
       when '3' => Generate_Random_System;
       when others => null;
    end case;
  end Main;

begin
  Main;
end ts_pretrop;
