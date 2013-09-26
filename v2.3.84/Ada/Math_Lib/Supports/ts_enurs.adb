with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Lists_of_Floating_Vectors;
with Lists_of_Floating_Vectors_io;       use Lists_of_Floating_Vectors_io;
with Lists_of_Integer_Vectors;
with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;
with Facet_Vertex_Enumeration;           use Facet_Vertex_Enumeration;

procedure ts_enurs is

-- DESCRIPTION :
--   Test on the enumeration of vertices and facets by reverse search.

  tol : constant double_float := 10.0**(-8);

  procedure Test_Vertex_Enumeration ( n,m : in integer32 ) is

    cff : Matrix(1..n,1..m);
    rhs : Vector(1..m);
    pts : Lists_of_Floating_Vectors.List;
    lbl : Lists_of_Integer_Vectors.List;
    fail,infty : boolean;
    ans : character;

  begin
    new_line;
    put("Give a "); put(n,1); put("*"); put(m,1);
    put_line("-matrix with coefficients of constraints in columns:");
    get(cff);
    put("Give a "); put(m,1);
    put_line("-vector with right-hand sides :"); get(rhs);
    Enumerate_Vertices(cff,rhs,tol,pts,lbl,fail,infty);
    if fail
     then put_line("The problem is unfeasible.");
    end if;
    if infty
     then put_line("The problem is unbounded.");
    end if;
   -- pts := Enumerate_Vertex_Points(cff,rhs,tol);
    put("Found "); put(Lists_of_Floating_Vectors.Length_Of(pts),1);
    put_line(" vertex points : "); put(pts);
   -- lbl := Enumerate_Vertex_Labels(cff,rhs,tol);
    put("Found "); put(Lists_of_Integer_Vectors.Length_Of(lbl),1);
    put_line(" vertex labels : "); put(lbl);
    put("Do you wish to continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        lenpts : constant integer32 
               := integer32(Lists_of_Floating_Vectors.Length_Of(pts));
        ptsmat : Matrix(1..n,1..lenpts) := List_to_Matrix(n,pts);
        fcs : Lists_of_Floating_Vectors.List;
      begin
        put_line("The matrix of vertex points : "); put(ptsmat,3);
        fcs := Enumerate_Facet_Inequalities(ptsmat,tol);
        put("Found "); put(Lists_of_Floating_Vectors.Length_Of(fcs),1);
        put_line(" facet inequalities : "); put(fcs);
      end;
    end if;
  end Test_Vertex_Enumeration;

  procedure Vertex_Enumeration is

    n,m : integer32 := 0;

  begin
    new_line;
    put("Give dimension, length of points : "); get(n);
    put("Give number of facet inequalities : "); get(m);
    Test_Vertex_Enumeration(n,m);
  end Vertex_Enumeration;

  procedure Test_Facet_Enumeration ( n,m : in integer32 ) is

    cff : Matrix(1..n,1..m);
    fcs : Lists_of_Floating_Vectors.List;
    lbl : Lists_of_Integer_Vectors.List;
    fail,infty : boolean;
    ans : character;

  begin
    new_line;
    put("Give a "); put(n,1); put("*"); put(m,1);
    put_line("-matrix with points in its columns :");
    get(cff);
    Enumerate_Facets(cff,tol,fcs,lbl,fail,infty);
    if fail
     then put_line("The problem is unfeasible.");
    end if;
    if infty 
     then put_line("The problem is unbounded.");
    end if;
   -- fcs := Enumerate_Facet_Inequalities(cff,tol);
    put("Found "); put(Lists_of_Floating_Vectors.Length_Of(fcs),1);
    put_line(" facet inequalities :"); put(fcs);
   -- lbl := Enumerate_Facet_Labels(cff,tol);
    put("Found "); put(Lists_of_Integer_Vectors.Length_Of(lbl),1);
    put_line(" point labels that span the facets :"); put(lbl);
    put("Do you wish to continue ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y' then
      declare
        lenfcs : constant integer32
               := integer32(Lists_of_Floating_Vectors.Length_Of(fcs));
        fcsmat : Matrix(1..n,1..lenfcs) := List_to_Matrix(n,fcs);
        fcsrhs : Vector(1..lenfcs) := List_to_Vector(n+1,fcs);
        pts : Lists_of_Floating_Vectors.List;
      begin
        put_line("The matrix of facet inequalities : "); put(fcsmat,3);
        put_line("The right-hand side vector : "); put(fcsrhs,3); new_line;
        pts := Enumerate_Vertex_Points(fcsmat,fcsrhs,tol);
        put("Found "); put(Lists_of_Floating_Vectors.Length_Of(pts),1);
        put_line(" vertex points : "); put(pts);
      end;
    end if;
  end Test_Facet_Enumeration;

  procedure Facet_Enumeration is

    n,m : integer32 := 0;

  begin
    new_line;
    put("Give dimension, length of points : "); get(n);
    put("Give the number of vertex points : "); get(m);
    Test_Facet_Enumeration(n,m);
  end Facet_Enumeration;

  procedure Test_Lower_Facet_Enumeration ( n,m : in integer32 ) is

    cff : Matrix(1..n,1..m);
    fcs : Lists_of_Floating_Vectors.List;
    lbl : Lists_of_Integer_Vectors.List;
    fail,infty : boolean;

  begin
    new_line;
    put("Give a "); put(n,1); put("*"); put(m,1);
    put_line("-matrix with points in its columns :");
    get(cff);
    Enumerate_Lower_Facets(cff,tol,fcs,lbl,fail,infty);
    if fail
     then put_line("The problem is unfeasible.");
    end if;
    if infty
     then put_line("The problem is unbounded.");
    end if;
   -- fcs := Enumerate_Lower_Facet_Inequalities(cff,tol);
    put("Found "); put(Lists_of_Floating_Vectors.Length_Of(fcs),1);
    put_line(" lower facet inequalities :"); put(fcs);
   -- lbl := Enumerate_Lower_Facet_Labels(cff,tol);
    put("Found "); put(Lists_of_Integer_Vectors.Length_Of(lbl),1);
    put_line(" point labels that span the lower facets :"); put(lbl);
  end Test_Lower_Facet_Enumeration;

  procedure Lower_Facet_Enumeration is

    n,m : integer32 := 0;

  begin
    new_line;
    put("Give dimension, length of points : "); get(n);
    put("Give the number of vertex points : "); get(m);
    Test_Lower_Facet_Enumeration(n,m);
  end Lower_Facet_Enumeration;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Facet and Vertex Enumeration by Reverse Search.");
    loop
      new_line;
      put_line("Choose one of the following options :");
      put_line("  0. exit this program");
      put_line("  1. enumerate vertices from given facet inequalities"); 
      put_line("  2. enumerate facet inequalities from given vertices"); 
      put_line("  3. enumerate lower facet inequalities from given vertices"); 
      put("Type 0,1,2 or 3 to select : "); get(ans);
      exit when (ans = '0');
      case ans is 
        when '1' => Vertex_Enumeration;
        when '2' => Facet_Enumeration;
        when '3' => Lower_Facet_Enumeration;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_enurs;
