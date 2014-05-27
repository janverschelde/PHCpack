with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers; 
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io; 
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Natural_VecVecs_io;        use Standard_Natural_VecVecs_io;
with Standard_Natural_Matrices;          use Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with Standard_Random_Matrices;
with Standard_Complex_Linear_Solvers;    use Standard_Complex_Linear_Solvers;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Checker_Boards,Checker_Boards_io;   use Checker_Boards,Checker_Boards_io;
with Checker_Moves;                      use Checker_Moves;
with Checker_Localization_Patterns;      use Checker_Localization_Patterns;
with Checker_Posets,Checker_Posets_io;   use Checker_Posets,Checker_Posets_io;
with Checker_Homotopies;
with Intersection_Posets;                use Intersection_Posets;
with Intersection_Posets_io;             use Intersection_Posets_io;
with Moving_Flag_Homotopies;             use Moving_Flag_Homotopies;

procedure ts_induce is

  procedure Generalizing_Moving_Flags ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Shows all boards, flags and transformations in a generalizing
  --   moving flag homotopy towards a general triangular matrix.

    fg : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Random_Flag(n);
    m : constant integer32 := integer32(Number_of_Moves(natural32(n)));
    all_moves : constant Standard_Natural_VecVecs.VecVec(1..m)
              := Checker_Posets.Generalizing_Moves(n);
    ans : character;
    p,p1 : Standard_Natural_Vectors.Vector(1..n);
    mf,t : Standard_Natural_Matrices.Matrix(1..n,1..n);
    nt : Standard_Complex_Matrices.Matrix(1..n,1..n);
    acc : Standard_Complex_Matrices.Matrix(1..n,1..n) := Identity(n);
    f,a : integer32;
    gamma : Complex_Number;

    use Standard_Complex_Matrices;

  begin
    put(all_moves);
    put("Do you wish to see all the boards ? (y/n) ");
    Ask_Yes_or_No(ans);
    put_line("A random flag :"); put(fg,3);
    if ans = 'y' then
      for i in all_moves'first..all_moves'last-1 loop
        p := all_moves(i).all; mf := Moving_Flag(p);
        f := Falling_Checker(p); a := Ascending_Checker(p,f);
        p1 := all_moves(i+1).all;
        t := Transformation(n,integer32(p1(f)));
        Write_Permutation(p,mf,t);
        put("After Move #"); put(i,1); 
        put(" (a,f) = ("); put(a,1); put(","); put(f,1);
        put(") : ");
        gamma := fg(n+1-a,f);
        nt := Numeric_Transformation(t,gamma);
        put_line("The numeric form of the transformation : "); put(nt,3);
        acc := acc*nt;
        put_line("The transformed moving flag : "); put(acc,3);
      end loop;
      put_line("The final configuration : ");
      p := all_moves(all_moves'last).all; mf := Moving_Flag(p);
      Write_Permutation(p,mf);
      put_line("The final flag : "); put(fg,3);
    end if;
  end Generalizing_Moving_Flags;

  procedure Test_Coordinates ( n,k : in integer32 ) is
 
  -- DESCRIPTION :
  --   Reads in a checker board and sets up the coordinates
  --   associated with the positions of white and black checkers.

    p : Vector(1..n);
    r,c : Vector(1..k);
    b : Board(1..n,1..n);
    f : Matrix(1..n,1..n);
    cm : Matrix(1..n,1..k);
    ht,ctr : integer32;

  begin
    put_line("Reading the permutation for the black checkers...");
    Read_Permutation(p);
    f := Moving_Flag(p);
    Write_Permutation(p,f);
    put_line("Reading rows and columns for the white checkers...");
    Read_Permutation(r);
    Read_Permutation(c);
    b := Configuration(p);
    Place_White(b,r,c);
    put_line("The board with black and white checkers : ");
    f := Moving_Flag(p);
    Write(b,f,p,r,c);
    cm := Column_Pattern(n,k,p,r,c);
    put_line("The pattern of the coordinate matrix : ");
    put(cm);
    Checker_Homotopies.Define_Specializing_Homotopy(n,p,r,c,ht,ctr);
  end Test_Coordinates;

  procedure Specializing_Homotopies ( n,k : in integer32 ) is

  -- DESCRIPTION :
  --   Walks through the poset starting at the root, printing information
  --   about the deformations involving a k-plane in n-space.

    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    f : Matrix(1..n,1..n) := Moving_Flag(p);
    b : Board(1..n,1..n) := Configuration(p);
    rows,cols : Vector(1..k);
    ps : Poset;
    ans : character;
    ptr : Link_to_Node;
    ht,ctr : integer32;
 
  begin
    put_line("The board with black checkers :");
    Write_Permutation(p,f);
    put_line("Reading rows and columns for white checkers...");
    Read_Permutation(rows);
    Read_Permutation(cols);
    Place_White(b,rows,cols);
    put_line("The board with black and white checkers :");
    Write(b,f,p,rows,cols);
    ps := Create(n,rows,cols);
    put_line("The complete poset : "); Write(ps);
    put_line("Running through the poset ...");
    for i in ps.black'first..ps.black'last-1 loop
      put(i,2); put(" : "); put(ps.black(i)); put(" : ");
      ptr := ps.white(i);
      while ptr /= null loop
        put("define homotopy for "); Write(ptr.rows,ptr.cols);
        new_line;
        b := Configuration(ps.black(i).all,ptr.rows,ptr.cols);
        f := Moving_Flag(ps.black(i).all);
        Write(b,f,ps.black(i).all,ptr.rows,ptr.cols);
        Checker_Homotopies.Define_Specializing_Homotopy
          (n,ps.black(i).all,ptr.rows,ptr.cols,ht,ctr);
        if ptr.next_sibling /= null then
          put("Continue to next sibling ? (y/n) ");
          Ask_Yes_or_No(ans);
          exit when (ans /= 'y');
        end if;
        ptr := ptr.next_sibling;
      end loop;
      if i < ps.black'last-1 then
        put("Continue to next level ? (y/n) ");
        Ask_Yes_or_No(ans);
        exit when (ans /= 'y');
      end if;
    end loop;
  end Specializing_Homotopies;

  procedure Descend_Poset
              ( ps : in Poset; n,k,leaf : in integer32;
                nd : in Link_to_Node ) is

  -- DESCRIPTION :
  --   Interactive descent in the poset from a leaf to the root.

  -- ON ENTRY :
  --   ps       poset created to deform k-planes in n-space;
  --   n        dimension of the ambient space;
  --   k        dimension of the subspace;
  --   leaf     label of a leaf to start at;
  --   nd       node at the leaf, should not be null.

    lvl : integer32 := ps.black'last;
    ind : integer32 := leaf;
    lnd : Link_to_Node := nd;
    stay : boolean;
    ans : character;
    np : natural32;
    kp,ht,ct : integer32 := 0;
    parent : Link_to_Node;
    b : Board(1..n,1..n);
    f : Matrix(1..n,1..n);

  begin
    loop
      new_line;
      put("Node "); put(ind,1); put(" at level "); put(lvl,1); put(" : ");
      put(ps.black(lvl).all); put(" : "); Write_Node(lnd.all); new_line;
      Write_Parents(ps.black(lvl-1).all,lnd.all);
      np := Number_of_Parents(lnd.all);
      put("#parents : "); put(np,1); new_line;
      if np = 1 then
        parent := lnd.first_parent;
      else
        loop
          put("  tell which parent to use (number <= ");
          put(np,1); put(") : "); get(kp);
          parent := Retrieve_Parent(lnd.all,kp);
          exit when (parent /= null);
          put_line("  wrong index to parent, please try again ...");
        end loop;
      end if;
      put_line("The configuration of checkers at child : ");
      b := Configuration(ps.black(lvl).all,lnd.rows,lnd.cols);
      f := Moving_Flag(ps.black(lvl).all);
      Write(b,f,ps.black(lvl).all,lnd.rows,lnd.cols);
      put_line("Localization pattern at child : ");
      put(Row_Pattern(n,k,ps.black(lvl).all,lnd.rows,lnd.cols));
      stay := Is_Stay_Child(parent.all,lnd.all);
      if stay
       then put_line("the child is a stay child");
       else put_line("the child is a swap child");
      end if;
      lnd := parent;
      lvl := lvl - 1;
      put_line("The configuration of checkers at parent : ");
      b := Configuration(ps.black(lvl).all,lnd.rows,lnd.cols);
      f := Moving_Flag(ps.black(lvl).all);
      Write(b,f,ps.black(lvl).all,lnd.rows,lnd.cols);
      put_line("Localization pattern at parent : ");
      put(Row_Pattern(n,k,ps.black(lvl).all,lnd.rows,lnd.cols));
      Checker_Homotopies.Define_Generalizing_Homotopy
        (n,ps.black(lvl).all,lnd.rows,lnd.cols,stay,ht,ct);
      exit when (lvl = 1);
      put("Move down to the next level ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      ind := Position(ps.white(lvl).all,lnd.all);
    end loop;
  end Descend_Poset;

  procedure Start_Descent_in_Poset ( n,k : in integer32; ps : in Poset ) is

  -- DESCRIPTION :
  --   Prompts the user to give a leaf in the poset to start the descent.

    ind : integer32 := 0;
    lnd : Link_to_Node;

  begin
    put_line("The complete poset : "); Write(ps);
    loop
      put("Type a number of a leaf to start deformations : "); get(ind);
      lnd := Retrieve(ps,ps.black'last,ind);
      exit when (lnd /= null);
      put("There is no leaf "); put(ind,1);
      put_line(".  Please try again...");
    end loop;
    Descend_Poset(ps,n,k,ind,lnd);
  end Start_Descent_in_Poset;

  procedure Generalizing_Homotopies ( n,k : in integer32 ) is

  -- DESCRIPTION :
  --   Creates a poset to deform k-planes in n-space.
  --   After prompting the user for the number of a leaf, when moving to 
  --   the root, information is printed about the homotopies.

    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    b : Board(1..n,1..n) := Configuration(p);
    f : Matrix(1..n,1..n) := Moving_Flag(p);
    rows,cols : Vector(1..k);
    ps : Poset;
 
  begin
    put_line("The board with black checkers :");
    Write_Permutation(p,f);
    put_line("Reading rows and columns for white checkers...");
    Read_Permutation(rows);
    Read_Permutation(cols);
    Place_White(b,rows,cols);
    f := Moving_Flag(p);
    put_line("The board with black and white checkers :");
    Write(b,f,p,rows,cols);
    ps := Create(n,rows,cols);
    Start_Descent_in_Poset(n,k,ps);
  end Generalizing_Homotopies;

  procedure Descend_Intersection_Poset
              ( n,k : in integer32; ips : in Intersection_Poset;
                leaf : in Link_to_Poset_Node ) is

  -- DESCRIPTION :
  --   Interactive descent in the intersection poset to the root,
  --   starting at the given leaf.

    cols : Vector(1..k);
    lvl : integer32 := ips.level;
    kp : integer32 := 0;
    pnd : Link_to_Poset_Node := leaf;
    lnd : Link_to_Node;
    ind : integer32;
    np : natural32;
    ans : character;

  begin
    loop
      if lvl = ips.level then
        Start_Descent_in_Poset(n,k,pnd.ps);
      else
        cols := Root_Rows(pnd.ps);
        if kp = 0
         then pnd := pnd.first_parent;
         else pnd := Retrieve_Parent(ips.nodes(lvl),pnd.all,kp);
        end if;
        Retrieve_Leaf(pnd.ps,cols,ind,lnd);
        if lnd /= null
         then Descend_Poset(pnd.ps,n,k,ind,lnd);
         else put("Leaf node with columns "); put(cols);
              put_line(" not found.  Abort descent."); return;
        end if;
      end if;
      exit when (pnd.first_parent = null);
      put("Current level is "); put(lvl,1);
      put(".  Descend further ? (y/n) ");  
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      lvl := lvl - 1;
      new_line;
      np := Number_of_Parents(ips.nodes(lvl),pnd.all);
      put("#parents : "); put(np,1); new_line;
      Write_Parents(ips.nodes(lvl),pnd.all);
      if np = 1 then
        put_line("Descending using first parent...");
        kp := 0;
      else
        loop
          put("  tell which parent to use (number <= ");
          put(np,1); put(") : "); get(kp);
          exit when (natural32(kp) <= np);
          put("  your number "); put(kp,1); put(" > "); put(np,1);
          put_line(" = #parents.  Please try again...");
        end loop;
      end if;
    end loop;
  end Descend_Intersection_Poset;

  procedure Intersecting_Homotopies ( n,k,m : in integer32 ) is

    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    rows,cols : Vector(1..k);
    ps : Poset;
    ips : Intersection_Poset(m-1);
    ind : integer32 := 0;
    leaf : Link_to_Poset_Node;

  begin
    put_line("Reading the first two intersection conditions...");
    Read_Permutation(rows); Read_Permutation(cols);
    if not Happy_Checkers(p,rows,cols) then
      put_line("Your conditions form an unhappy configuration.");
      return;
    else
      ps := Create(n,rows,cols);
      ips := Create(m-1,ps);
      for k in 3..m loop
        put("Reading intersection condition "); put(k,1); put_line("...");
        Read_Permutation(cols);
        Intersect(ips,cols,false);
      end loop;
    end if;
    put_line("All formal equations in the intersection poset :");
    Write_Formal_Equations(ips);
    put_line("The intersection condition resolved :");
    Write_Expansion(ips);
    loop
      put("Number of leaves in the intersecting poset : ");
      put(Length_Of(ips.nodes(ips.level)),1); new_line;
      put("Give a leaf to start the descent at : "); get(ind);
      leaf := Retrieve(ips.nodes(ips.level),ind);
      exit when (leaf /= null);
      put("There is no leaf "); put(ind,1);
      put_line(".  Please try again...");
    end loop;
    Descend_Intersection_Poset(n,k,ips,leaf);
  end Intersecting_Homotopies;

  procedure Localization_along_Paths ( n,k : integer32; ps : in Poset ) is

  -- DESCRIPTION :
  --   Writes for every node along the paths in the poset
  --   the localization pattern.

    cnt : integer32 := 0;

    procedure Write_Path ( nds : in Array_of_Nodes; ct : out boolean ) is
    begin
      cnt := cnt + 1;
      Checker_Posets_io.Write_Path_in_Poset(n,k,ps,nds,cnt);
      ct := true;
    end Write_Path;

    procedure Write_Paths is new Enumerate_Paths_in_Poset(Write_Path);

  begin
    Write_Paths(ps);
  end Localization_along_Paths;

  procedure Run_along_Paths_in_Poset ( n,k : integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a configuration of white checkers,
  --   creates a poset and writes all paths through the poset.
  --   For each node the localization pattern is written.

    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    rows,cols : Vector(1..k);
    ps : Poset;

  begin
    new_line;
    loop
      put_line("Reading two intersection conditions...");
      Read_Permutation(rows); Read_Permutation(cols);
      exit when Happy_Checkers(p,rows,cols);
      put("Your conditions form an unhappy configuration.");
      put_line("  Please try again...");
    end loop;
    ps := Create(n,rows,cols);
    put_line("The poset : "); Checker_Posets_io.Write(ps);
    put_line("Enumerating all paths : ");
    Localization_along_Paths(n,k,ps);
  end Run_along_Paths_in_Poset;

  function Flag_Transformation_Coefficient_Matrix
             ( n : integer32;
               f1,f2,g1,g2 : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coefficient matrix to compute the matrix A
  --   and the upper triangular matrices T1 and T2 in the transformation
  --   defined by A*f1 = g1*T1, A*f2 = g2*T2.
  --   The square matrix on return has dimension 2*n^2.

    dim : constant integer32 := 2*n*n;
    res : Standard_Complex_Matrices.Matrix(1..dim,1..dim);
    offset : integer32 := 0;

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    for k in 0..(n-1) loop
      for i in 1..n loop
        for j in 1..n loop
          res(k*n+i,k*n+j) := f1(j,i);
          res(n*(n+k)+i,k*n+j) := f2(j,i);
        end loop;
      end loop;
    end loop;
    for i in 1..n-1 loop
      for j in 1..n-i loop
        for k in 0..(n-1) loop
          put("i = "); put(i,1);
          put("  j = "); put(j,1); put("  k = "); put(k,1);
          put("  row : "); put(k*n+i+j,1);
          put("  col : "); put(offset+n*n+j,1); new_line;
          res(k*n+i+j,offset+n*n+j) := -g1(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i;
    end loop;
    new_line;
    for i in 1..n loop
      for j in 1..n-i+1 loop
        for k in 0..(n-1) loop
          put("i = "); put(i,1);
          put("  j = "); put(j,1); put("  k = "); put(k,1);
          put("  row : "); put(n*(n+k)+i+j-1,1);
          put("  col : "); put(offset+n*n+j,1); new_line;
          res(n*(n+k)+i+j-1,offset+n*n+j) := -g2(k+1,i);
        end loop;
      end loop;
      offset := offset + n - i + 1;
    end loop;
    return res;
  end Flag_Transformation_Coefficient_Matrix;

  function Flag_Transformation_Right_Hand_Side
             ( n : integer32; g : Standard_Complex_Matrices.Matrix )
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the right hand side vector in the linear system to
  --   compute a transformation between two flags.

    dim : constant integer32 := 2*n*n;
    res : Standard_Complex_Vectors.Vector(1..dim);
    ind : integer32 := 0;

  begin
    for i in g'range(1) loop
      for j in g'range(2) loop
        ind := ind + 1;
        res(ind) := g(i,j);
      end loop;
    end loop;
    for i in ind+1..res'last loop
      res(i) := Create(0.0);
    end loop;
    return res;
  end Flag_Transformation_Right_Hand_Side;

  procedure real_small_put ( A : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Writes the real part of the matrix A with only one decimal place
  --   of precision, enough to see the structure of the matrix.

    rp : double_float;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        rp := REAL_PART(A(i,j));
        if rp < 0.0
         then put(" "); 
         else put("  ");
        end if;
        put(REAL_PART(A(i,j)),1);
      end loop;
      new_line;
    end loop;
  end real_small_put;

  procedure Extract_Matrices
              ( n : in integer32; sol : in Standard_Complex_Vectors.Vector;
                A,T1,T2 : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Given a vector sol of size 2*n*n, extracts three n-dimensional 
  --   matrices: A, T1, and T2, where T1 and T2 are upper triangular
  --   and T1 has ones on its diagonal.

    idx : integer32 := 0;

  begin
    for i in 1..n loop
      for j in 1..n loop
        idx := idx + 1;
        A(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T1(i,j) := Create(0.0);
      end loop;
      T1(i,i) := Create(1.0);
      for j in (i+1)..n loop
        idx := idx + 1;
        T1(i,j) := sol(idx);
      end loop;
    end loop;
    for i in 1..n loop
      for j in 1..(i-1) loop
        T2(i,j) := Create(0.0);
      end loop;
      for j in i..n loop
        idx := idx + 1;
        T2(i,j) := sol(idx);
      end loop;
    end loop;
  end Extract_Matrices;

  procedure Test_Flag_Transformation
              ( f1,f2,g1,g2,A,T1,T2 : in Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Tests whether A*f1 = g1*T1 and A*f2 = g2*T2.

    use Standard_Complex_Matrices;

  begin
    put_line("A*f1 : "); put(A*f1);
    put_line("g1*T1 : "); put(g1*T1);
    put_line("A*f2 : "); put(A*f2);
    put_line("g2*T2 : "); put(g2*T2);
  end Test_Flag_Transformation;

  procedure Test_Random_Flags2Flags ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates two pairs of random flags (F1, F2) and (G1, G2) in n-space
  --   and computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.

    f1 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
    f2 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
    g1 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
    g2 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
    dim : constant integer32 := 2*n*n;
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Flag_Transformation_Coefficient_Matrix(n,f1,f2,g1,g2);
    rhs : constant Standard_Complex_Vectors.Vector(1..dim)
        := Flag_Transformation_Right_Hand_Side(n,g1);
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    sol : Standard_Complex_Vectors.Vector(1..dim) := rhs;
    res : Standard_Complex_Vectors.Vector(1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    nrm : double_float;
    A,T1,T2 : Standard_Complex_Matrices.Matrix(1..n,1..n);

    use Standard_Complex_Vectors,Standard_Complex_Matrices;

  begin
    put_line("The flag f1 : "); put(f1,1); new_line;
    put_line("The flag f2 : "); put(f2,1); new_line;
    put_line("The flag g1 : "); put(g1,1); new_line;
    put_line("The flag g2 : "); put(g2,1); new_line;
    put_line("The coefficient matrix "); real_small_put(mat); new_line;
    lufac(wrk,dim,ipvt,info);
    lusolve(wrk,dim,ipvt,sol);
    res := rhs - mat*sol;
    nrm := Max_Norm(res);
    put("The norm of the residual : "); put(nrm); new_line;
    Extract_Matrices(n,sol,A,T1,T2);
    Test_Flag_Transformation(f1,f2,g1,g2,A,T1,T2);
  end Test_Random_Flags2Flags;

  function Moved_Flag
             ( n : integer32 ) return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the coordinates for the moved flag in n-space.

    res : Standard_Complex_Matrices.Matrix(1..n,1..n);

  begin
    for i in 1..n loop
      for j in 1..n-i+1 loop
        res(i,j) := Create(1.0);
      end loop;
      for j in (n-i+2)..n loop
        res(i,j) := Create(0.0);
      end loop;
    end loop;
    return res;
  end Moved_Flag;

  procedure Test_Specific_Flags2Flags ( n : in integer32 ) is

  -- DESCRIPTION :
  --   For two pairs of flags (M, I) and (I, F) in n-space
  --   where M is the moved flags, I the identity matrix,
  --   and F some random n-dimensional matrix,
  --   computes the matrix A, upper triangular matrices T1 and T2
  --   so that A*F1 = G1*T1 and A*F2 = G2*T2.

    f1 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Moved_Flag(n);
    f2 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Moving_Flag_Homotopies.Identity(n);
    g1 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Moving_Flag_Homotopies.Identity(n);
    g2 : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
       := Standard_Random_Matrices.Random_Matrix(natural32(n),natural32(n));
    dim : constant integer32 := 2*n*n;
    mat : constant Standard_Complex_Matrices.Matrix(1..dim,1..dim)
        := Flag_Transformation_Coefficient_Matrix(n,f1,f2,g1,g2);
    rhs : constant Standard_Complex_Vectors.Vector(1..dim)
        := Flag_Transformation_Right_Hand_Side(n,g1);
    wrk : Standard_Complex_Matrices.Matrix(1..dim,1..dim) := mat;
    sol : Standard_Complex_Vectors.Vector(1..dim) := rhs;
    res : Standard_Complex_Vectors.Vector(1..dim);
    ipvt : Standard_Integer_Vectors.Vector(1..dim);
    info : integer32;
    nrm : double_float;
    A,T1,T2 : Standard_Complex_Matrices.Matrix(1..n,1..n);

    use Standard_Complex_Vectors,Standard_Complex_Matrices;

  begin
    put_line("The flag f1 : "); put(f1,1); new_line;
    put_line("The flag f2 : "); put(f2,1); new_line;
    put_line("The flag g1 : "); put(g1,1); new_line;
    put_line("The flag g2 : "); put(g2,1); new_line;
    put_line("The coefficient matrix "); real_small_put(mat); new_line;
    lufac(wrk,dim,ipvt,info);
    lusolve(wrk,dim,ipvt,sol);
    res := rhs - mat*sol;
    nrm := Max_Norm(res);
    put("The norm of the residual : "); put(nrm); new_line;
    Extract_Matrices(n,sol,A,T1,T2);
    Test_Flag_Transformation(f1,f2,g1,g2,A,T1,T2);
  end Test_Specific_Flags2Flags;

  procedure Main is
 
    n,k,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Choose one of the following operations : ");
    put_line("  1. see generalizing moving flag in numeric form;");
    put_line("  2. test homotopy definition for specific coordinates;");
    put_line("  3. create a poset and run through specialization;");
    put_line("  4. create a poset and deform from leaves to the root;");
    put_line("  5. let intersection poset define the deformations;");
    put_line("  6. run along paths in a checker poset;");
    put_line("  7. test transformation of pair of flags.");
    put("Type 1, 2, 3, 4, 5, 6, or 7 to make your choice : ");
    Ask_Alternative(ans,"1234567");
    new_line;
    put("Give n (ambient dimension) : "); get(n);
    if ans /= '1' and ans /= '7'
     then put("Give k (subspace dimension): "); get(k);
    end if;
    case ans is
      when '1' => Generalizing_Moving_Flags(n);
      when '2' => Test_Coordinates(n,k);
      when '3' => Specializing_Homotopies(n,k);
      when '4' => Generalizing_Homotopies(n,k);
      when '5' => new_line;
                  put("Give the total number of intersection conditions : ");
                  get(m);
                  Intersecting_Homotopies(n,k,m);
      when '6' => Run_along_Paths_in_Poset(n,k);
      when '7' => Test_Random_Flags2Flags(n);
                  Test_Specific_Flags2Flags(n);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_induce;
