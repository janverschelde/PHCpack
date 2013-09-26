with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;          use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;          use Standard_Natural_VecVecs;
with Standard_Natural_VecVecs_io;       use Standard_Natural_VecVecs_io;
with Standard_Natural_Matrices;         use Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;      use Standard_Natural_Matrices_io;
with Brackets,Brackets_io;              use Brackets,Brackets_io;
with Bracket_Monomials;                 use Bracket_Monomials;
with Bracket_Monomials_io;              use Bracket_Monomials_io;
with Checker_Boards,Checker_Boards_io;  use Checker_Boards,Checker_Boards_io;
with Checker_Moves;                     use Checker_Moves;
with Checker_Posets,Checker_Posets_io;  use Checker_Posets,Checker_Posets_io;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;
with Intersection_Posets;               use Intersection_Posets;
with Intersection_Posets_io;            use Intersection_Posets_io;

procedure ts_checkers is

-- DESCRIPTION :
--   Simple test on the operations in Checker_Posets.

-- TARGET TEST ROUTINES :

  procedure Interactive_Test ( n : in integer32 ) is

    p : Vector(1..n);
    ans : character;
    mf : Matrix(1..n,1..n);

  begin
    loop
      Read_Permutation(p);
      mf := Moving_Flag(p);
      Write_Permutation(p,mf);
      put("index of descending checker : ");
      put(Descending_Checker(p),1); new_line;
      put("index of rising checker : ");
      put(Rising_Checker(p,Descending_Checker(p)),1); new_line;
      put("Want another permutation ? (y/n) ");
      get(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Test;

  procedure Specialization_Order ( n : in integer32) is

    p : Vector(1..n) := Identity_Permutation(natural32(n));
    mf : Matrix(1..n,1..n);
    ans : character;
    desc,rise : integer32 := 0;
    cnt : natural32 := 0;

  begin
    put("Number of moves : "); 
    put(Number_of_Moves(natural32(n)),1); new_line;
    loop
      cnt := cnt + 1;
      put("permutation no ");
      put(cnt,1); put_line(" : ");
      mf := Moving_Flag(p);
      Write_Permutation(p,mf);
      put("index of descending checker : ");
      desc := Descending_Checker(p);
      put(desc,1); new_line;
      exit when (desc = 0);
      rise := Rising_Checker(p,desc);
      put("index of rising checker : ");
      put(rise,1); new_line;
      exit when (rise = 0);
      put("Do we make the move ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Move(p,desc,rise);
    end loop;
  end Specialization_Order;

  procedure Show_Specializing_Moves ( n : in integer32 ) is

    m : constant natural32 := Number_of_Moves(natural32(n));
    all_moves : constant VecVec(1..integer32(m)) := Specializing_Moves(n);
    ans : character;
    p : Vector(1..n);
    mf,t : Matrix(1..n,1..n);
    r,d : integer32;

  begin
    put(all_moves);
    put("Do you wish to see all the boards ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Initial configuration :");
      p := all_moves(all_moves'first).all;
      mf := Moving_Flag(p);
      d := Descending_Checker(p);
      r := Rising_Checker(p,d);
      Write_Permutation(p,mf);
      for i in all_moves'first+1..all_moves'last loop
        put("After Move #"); put(i-1,1); 
        put(" (r,d) = ("); put(r,1); put(","); put(d,1);
        put_line(") :");
        t := Transformation(n,integer32(p(d)));
        p := all_moves(i).all;
        mf := Moving_Flag(p);
        d := Descending_Checker(p);
        Write_Permutation(p,mf,t);
        exit when (d = 0);
        r := Rising_Checker(p,d);
      end loop;
    end if;
  end Show_Specializing_Moves;

  procedure Show_Generalizing_Moves ( n : in integer32 ) is

    m : constant natural32 := Number_of_Moves(natural32(n));
    all_moves : constant VecVec(1..integer32(m)) := Generalizing_Moves(n);
    ans : character;
    p,p1 : Vector(1..n);
    mf,t : Matrix(1..n,1..n);
    f,a : integer32;

  begin
    put(all_moves);
    put("Do you wish to see all the boards ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("Initial configuration :");
      p := all_moves(all_moves'first).all;
      mf := Moving_Flag(p);
      f := Falling_Checker(p);
      a := Ascending_Checker(p,f);
      p1 := all_moves(all_moves'first+1).all;
      t := Transformation(n,integer32(p1(f)));
      Write_Permutation(p,mf,t);
      for i in all_moves'first+1..all_moves'last loop
        put("After Move #"); put(i-1,1); 
        put(" (a,f) = ("); put(a,1); put(","); put(f,1);
        put_line(") :");
        p := all_moves(i).all;
        mf := Moving_Flag(p);
        f := Falling_Checker(p);
        if f = 0
         then Write_Permutation(p,mf);
        end if;
        exit when (f = 0);
        a := Ascending_Checker(p,f);
        p1 := all_moves(i+1).all;
        t := Transformation(n,integer32(p1(f)));
        Write_Permutation(p,mf,t);
      end loop;
    end if;
  end Show_Generalizing_Moves;

  procedure Show_All_Moves ( n : in natural32 ) is

    m : constant natural32 := Number_of_Moves(n);
    ans : character;

  begin
    put("The number of configurations : "); put(m,1); new_line;
    put("specializing or generalizing ? (s/g) ");
    Ask_Alternative(ans,"sg");
    if ans = 's'
     then Show_Specializing_Moves(integer32(n));
     else Show_Generalizing_Moves(integer32(n));
    end if;
  end Show_All_Moves;

  procedure Initialize_Checkerboard
               ( p,rows,cols : out Vector; b : out Board ) is

  -- DESCRIPTION :
  --   Initialization of the checkerboard for k-planes in n-space.
  --   Checks the happiness of the white checkers and forces the
  --   user to re-enter the rows and columns of the white checkers
  --   in case the white checkers are not happy.

  -- ON RETURN :
  --   p         permutation defines the n black checkers;
  --   rows      rows of the k white checkers;
  --   cols      columns of the k white checkers;
  --   b         a checkerboard with n black and k white checkers.

    happy : boolean;
    mf : Matrix(p'range,p'range);

  begin
    put_line("Reading coordinates for black checkers...");
    Read_Permutation(p);
    b := Configuration(p);
    put_line("The board with black checkers :");
    mf := Moving_Flag(p);
    Write_Permutation(p,mf);
    loop
      put_line("Reading rows and columns for white checkers...");
      Read_Permutation(rows);
      Read_Permutation(cols);
      Check_Happiness(p,rows,cols,happy);
      exit when happy;
      put_line("The white checkers are not happy.  Please try again...");
    end loop;
    Place_White(b,rows,cols);
    put_line("The board with black and white checkers :");
    Write(b,p,rows,cols);
  end Initialize_Checkerboard;

  procedure Test_White_Moves ( k,n : in integer32 ) is

    p : Vector(1..n);
    b : Board(1..n,1..n);
    rows,cols : Vector(1..k);
    r,d,cr,cd : integer32;
    ans : character;
    debug : constant boolean := true;
    we_swap : boolean;

  begin
    Initialize_Checkerboard(p,rows,cols,b);
    loop
      d := Descending_Checker(p);
      exit when (d = 0);
      we_swap := false;
      if debug then
        put("Descending black checker is in row ");
        put(p(d),1); put(" and column "); put(p'last-d+1,1);
        put_line(".");
      end if;
      r := Rising_Checker(p,d);
      exit when (r = 0);
      if debug then
        put("Rising black checker is in row ");
        put(p(r),1); put(" and column "); put(p'last-r+1,1);
        put_line(".");
        put("The critical row is "); put(p(d),1);
        put_line(".");
      end if;
      cr := Critical_Row(integer32(p(d)),p'last-d+1,rows,cols);
      if debug then
        if cr = 0 then
          put_line("There is no white checker in the critical row.");
          put_line("The white checkers stay at their current position.");
        else
          put("White checker at row "); put(rows(cr),1);
          put(" and column "); put(cols(cols'last-cr+1),1);
          put_line(" is in the critical row.");
          put("Start of critical diagonal : ");
          put(p(r),1); put(" "); put(n-r+1,1); new_line;
        end if;
      end if;
      if cr /= 0 then
        cd := Top_White_Checker(integer32(p(r)),n-r+1,n,rows,cols);
        if debug then
          if cd = 0 then
            put_line("No white checker is in the critical diagonal.");
            cr := 0;  -- forces a stay
          else
            put("Top white checker is in row ");
            put(rows(cd),1); put(" and column ");
            put(cols(cols'last-cd+1),1); new_line;
          end if;
        end if;
      end if;
      if cr > 0 and cd > 0 then
        case Central_Choice(p,d,r,rows,cols,cr,cd,true) is
          when 1 => we_swap := true;
                    put_line("There is no choice, we swap.");
          when 2 => put("Type 0 to stay, or 1 to swap : ");
                    Ask_Alternative(ans,"01");
                    we_swap := (ans = '1');
          when others => put_line("There is no choice, we stay.");
        end case;
      end if;
      if we_swap
       then Swap(rows,cr,cd);
      end if;
      Move(p,d,r);
      if debug then
        b := Configuration(p,rows,cols);
        put_line("The board after moves before making happy :");
        Write(b,p,rows,cols);
      end if;
      Make_Happy(p,rows,cols,true);
      b := Configuration(p,rows,cols);
      put_line("The board after all moves of black and white :");
      Write(b,p,rows,cols);
      put("Continue with next move ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
    end loop;
  end Test_White_Moves;

  procedure Group_Black_Checkers ( k,n : in integer32 ) is

    p : Vector(1..n);
    b : Board(1..n,1..n);
    rows,cols : Vector(1..k);
    r,d,cd,cr : integer32;
    debug : constant boolean := true;
    degenerating : constant boolean := false;

  begin
    Initialize_Checkerboard(p,rows,cols,b);
    if degenerating then
      d := Descending_Checker(p);
      if d /= 0 then
        if debug then
          put("Descending black checker is in row ");
          put(p(d),1); put(" and column "); put(p'last-d+1,1);
          put_line(".");
        end if;
        r := Rising_Checker(p,d);
        if r /= 0 then
          if debug then
            put("Rising black checker is in row ");
            put(p(r),1); put(" and column "); put(p'last-r+1,1);
            put_line(".");
            put("The critical row is "); put(p(d),1);
            put_line(".");
          end if;
          cr := Critical_Row(integer32(p(d)),p'last-d+1,rows,cols);
          if debug then
            if cr = 0 then
              put_line("There is no white checker in the critical row.");
              put_line("The white checkers stay at their current position.");
            else
              put("White checker at row "); put(rows(cr),1);
              put(" and column "); put(cols(cols'last-cr+1),1);
              put_line(" is in the critical row.");
              put("Start of critical diagonal : ");
              put(p(r),1); put(" "); put(n-r+1,1); new_line;
            end if;
          end if;
          if cr /= 0 then
            cd := Top_White_Checker(integer32(p(r)),n-r+1,n,rows,cols);
            if debug then
              if cd = 0 then
                put_line("No white checker is in the critical diagonal.");
                cr := 0;  -- forces a stay
              else
                put("Top white checker is in row ");
                put(rows(cd),1); put(" and column ");
                put(cols(cols'last-cd+1),1); new_line;
              end if;
            end if;
          end if;
          if cr > 0 and cd > 0 then
            case Central_Choice(p,d,r,rows,cols,cr,cd,true) is
              when 1 => put_line("There is no choice: must swap.");
              when 2 => put_line("There is choice between swap and stay.");
              when others => put_line("There is no choice: must stay.");
            end case;
          end if;
        end if;
      end if;
    else  -- generalizing move
      d := Falling_Checker(p);
      if d /= 0 then
        r := Ascending_Checker(p,d);
        if debug then
          put("Descending black checker is in row ");
          put(p(d)-1,1); put(" and column "); put(p'last-r+1,1);
          put_line(".");
        end if;
        if r /= 0 then
          if debug then
            put("Rising black checker is in row ");
            put(p(r)+1,1); put(" and column "); put(p'last-d+1,1);
            put_line(".");
            put("The critical row is "); put(p(d)-1,1);
            put_line(".");
          end if;
        end if;
        cr := integer32(p(d)) - 1;
        put_line("Checkers above critical row :");
        for i in p'range loop
          if integer32(p(i)) < cr then
            put("  checker "); put(p(i),1);
            if p'last+1-i < p'last-r+1
             then put_line(" is at left of descending checker");
            end if;
            if p'last+1-i > p'last-d+1
             then put_line(" is at right of rising checker");
            end if;
          elsif integer32(p(i)) > cr+1 then
            put("  checker "); put(p(i),1);
            put_line(" is below moving checkers");
          end if;
        end loop;
      end if;
    end if;
  end Group_Black_Checkers;

  procedure Write_Children ( p : in Vector; nd : in Node ) is

  -- DESCRIPTION :
  --   Writes the children of the node lnd, also displaying the
  --   configuration of black checkers encoded in the permutation p.

  begin
    if nd.stay_child /= null then
      put("  the stay child  : "); put(p); put(" : ");
      Write_Node(nd.stay_child.all); new_line;
    end if;
    if nd.swap_child /= null then
      put("  the swap child  : "); put(p); put(" : ");
      Write_Node(nd.swap_child.all); new_line;
    end if;
  end Write_Children;

  procedure Pick_a_Node ( n,k : in integer32; ps : in Poset ) is

  -- DESCRIPTION :
  --   Prompts the user for coordinates (level and position) of a node
  --   in the poset.  All relevant information of that node is displayed.

     i,j : integer32 := 0;
     lnd : Link_to_Node;
     ans : character := 'y';
     b : Board(1..n,1..n);
     f : Matrix(1..n,1..n);

  begin
    put_line("Retrieving nodes from the poset ...");
    loop
      put("  give the level of the node : "); get(i);
      put("  give the position of the node : "); get(j);
      lnd := Retrieve(ps,i,j);
      if lnd = null then
        put("There is no node at position "); put(j,1);
        put(" at level "); put(i,1); put_line(".  Please try again.");
      else
        put("Node "); put(j,1); put(" at level "); put(i,1); put(" : ");
        put(ps.black(i).all); put(" : ");  Write_Node(lnd.all); new_line;
        if i > ps.black'first then
          put("  number of parents : "); 
          put(Number_of_Parents(lnd.all),1); new_line;
          Write_Parents(ps.black(i-1).all,lnd.all);
        end if;
        if i < ps.black'last
         then Write_Children(ps.black(i+1).all,lnd.all);
        end if;
        put_line("The configuration of checkers at the node : ");
        b := Configuration(ps.black(i).all,lnd.rows,lnd.cols);
        f := Moving_Flag(ps.black(i).all);
        Write(b,f,ps.black(i).all,lnd.rows,lnd.cols);
        put_line("The localization pattern at the node :");
        put(Row_Pattern(n,k,ps.black(i).all,lnd.rows,lnd.cols));
        put("Do you wish to retrieve more nodes ? (y/n) ");
        Ask_Yes_or_No(ans);
      end if;
      exit when (ans /= 'y');
    end loop;
  end Pick_a_Node;

  procedure Create_Poset ( k,n : in integer32 ) is

  -- DESCRIPTION :
  --   Creates a poset for k-planes in n-space.

    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    b : Board(1..n,1..n) := Configuration(p);
    mf : constant Matrix(1..n,1..n) := Moving_Flag(p);
    rows,cols : Vector(1..k);
    ps : Poset;
    ans : character;
    happy : boolean;

  begin
    put_line("The board with black checkers :");
    Write_Permutation(p,mf);
    loop
      put_line("Reading coordinates for white checkers...");
      Read_Permutation(rows); Read_Permutation(cols);
      Check_Happiness(p,rows,cols,happy);
      exit when happy;
      put_line("The white checkers are not happy.  Please try again...");
    end loop;
    Place_White(b,rows,cols);
    put_line("The board with black and white checkers :");
    Write(b,mf,p,rows,cols);
    ps := Create(n,rows,cols);
    Write_Formal_Sums(ps);
    put_line("The complete poset : "); Write(ps);
    put("Do you wish to retrieve selected nodes ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Pick_a_Node(n,k,ps);
    end if;
    put("Do you wish to see all localization patterns ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      put_line("The localization patterns in the complete poset : ");
      Write_Patterns(ps);
    end if;
  end Create_Poset;

  procedure Create_Intersection_Poset
              ( n,m : in integer32; p,rows,cols : in Vector ) is

  -- DESCRIPTION :
  --   Creates an intersection poset for at most m intersection
  --   conditions on k-planes in n-space.  The initial configuration
  --   of black checkers is in the permutation p and the positions of
  --   the white checkers is defined by rows and cols.

    ps : Poset;
    ips : Intersection_Poset(m-1);
    d : natural32;
    ans : character;
    w : Vector(rows'range);

  begin
    ps := Create(n,rows,cols);              -- create the first poset
    Write_Formal_Equation(ps);
    d := Degree_of_Freedom(p,rows,cols);
    put("The degree of freedom : "); put(d,1); new_line;
    ips.level := 0;
    while (d > 0) and (ips.level < m-1) loop
      put("Do you want to add more conditions ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when ans /= 'y';
      if ans = 'y' then   
        if ips.level = 0
         then ips := Create(m-1,ps);        -- first poset becomes root
        end if;
        put_line("Reading an intersection condition...");
        Read_Permutation(w);
        Intersect(ips,w);
        put_line("The new formal equations : ");
        Write_Formal_Equations(ips,ips.level);
        d := Degree_of_Freedom(ips);
        put("The degree of freedom at the leaf posets : ");
        put(d,1); new_line;
      end if;
    end loop;
    put_line("All formal equations in the intersection poset :");
    Write_Formal_Equations(ips);
    put_line("The intersection condition resolved :");
    Write_Expansion(ips);
  end Create_Intersection_Poset;

  procedure Count_Planes ( k,n : in integer32 ) is

  -- DESCRIPTION :
  --   Counts the number of k-planes in n-space which satisfy
  --   certain user given Schubert conditions.

    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    rows,cols : Vector(1..k);
    m : integer32 := 0;

  begin
    loop
      put("Give the total number of intersection conditions : "); get(m);
      exit when (m >= 2);
      put("Your number "); put(m,1); put(" must be larger than 1.");
      put_line("  Please try again...");
    end loop;
    loop
      put_line("Reading two intersection conditions...");
      Read_Permutation(rows); Read_Permutation(cols);
      exit when Happy_Checkers(p,rows,cols);
      put("Your conditions form an unhappy configuration.");
      put_line("  Please try again...");
    end loop;
    Create_Intersection_Poset(n,m,p,rows,cols);
  end Count_Planes;

  procedure Create_Intersection_Poset
              ( n : in integer32; bm : in Bracket_Monomial ) is

    nb : constant integer32 := integer32(Number_of_Brackets(bm));
    cd : constant Array_of_Brackets(1..nb) := Create(bm);
    k : constant integer32 := cd(1)'last;
    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    r,c,w : Vector(1..k);
    ps : Poset;
    ips : Intersection_Poset(nb-1);

  begin
    put("  the dimension of the planes : "); put(k,1); new_line;
    put("  the number of conditions : "); put(nb,1); new_line;
    if nb >= 2 then
      r := Standard_Natural_Vectors.Vector(cd(1).all);
      c := Standard_Natural_Vectors.Vector(cd(2).all);
      put(cd(1).all); put(" and "); put(cd(2).all);
      if not Happy_Checkers(p,r,c) then
        put_line(" are not happy together.  Please try again.");
      else
        put_line(" form a happy configuration.");
        ps := Create(n,r,c);              -- create the first poset
        Write_Formal_Equation(ps);
        ips := Create(nb-1,ps);           -- first poset becomes root
        for i in 3..nb loop
          w := Standard_Natural_Vectors.Vector(cd(i).all);
          Intersect(ips,w);
          put_line("The new formal equations : ");
          Write_Formal_Equations(ips,ips.level);
        end loop;
        put_line("All formal equations in the intersection poset :");
        Write_Formal_Equations(ips);
        put_line("The intersection condition resolved :");
        Write_Expansion(ips);
      end if;
    end if;
  end Create_Intersection_Poset;

  procedure Resolve_Intersection_Condition ( n : in integer32 ) is

    bm : Bracket_Monomial;
    nn : integer32 := n;
    ans : character;

  begin
    loop
      new_line;
      put_line("A product of brackets is for example ");
      put_line("  [2 4 6]*[2 4 6]*[2 4 6]; or [2 4 6]^3;");
      put_line("Give a product of brackets, terminate by a semicolon ;");
      get(bm);
      put("Your product : "); put(bm); new_line;
      Create_Intersection_Poset(nn,bm);
      new_line;
      put("Resolve another condition ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      put("The current ambient dimension is "); put(nn,1); put_line(".");
      put("-> change ambient dimension ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then put("Enter new ambient dimension : "); get(nn);
      end if;
    end loop;
  end Resolve_Intersection_Condition;

  procedure Write_Paths_in_Poset ( ps : in Poset ) is

    cnt : natural32 := 0;

    procedure Write_Path ( nds : in Array_of_Nodes; ct : out boolean ) is
    begin
      cnt := cnt + 1;
      put("Path "); put(cnt,1); put(" starts at ");
      Checker_Posets_io.Write_Node(nds(nds'first).all); new_line;
      for i in nds'first+1..nds'last loop
        put(" -> "); put(i,1); put(" : ");
        Checker_Posets_io.Write_Node(nds(i).all); new_line;
      end loop;
      ct := true;
    end Write_Path;

    procedure Write_Paths is new Enumerate_Paths_in_Poset(Write_Path);

  begin
    Write_Paths(ps);
  end Write_Paths_in_Poset;

  procedure Walk_to_First_Parent ( ps : in Poset ) is

  -- DESCRIPTION :
  --   Walks from leaves to the root of the poset,
  --   using only the first parent.

    leaf : Link_to_Node := ps.white(ps.white'last);
    walk : Link_to_Node;
    np : natural32;
    i : integer32 := 1;

  begin
    loop
      put("Leaf "); put(i,1); put(" :");
      Checker_Posets_io.Write_Node(leaf.all);
      np := Number_of_Parents(leaf.all);
      put(" has "); put(np,1); put_line(" parents.");
      walk := leaf.first_parent;
      while walk /= null loop
        put("Node"); 
        Checker_Posets_io.Write_Node(walk.all);
        if walk.stay_child /= null then 
          put(" has stay child");
          Checker_Posets_io.Write_Node(walk.stay_child.all);
        end if;
        if walk.swap_child /= null then 
          put(" has swap child");
          Checker_Posets_io.Write_Node(walk.swap_child.all);
        end if;
        np := Number_of_Parents(walk.all);
        put("  #parents : "); put(np,1); new_line;
        walk := walk.first_parent;
      end loop;
      leaf := leaf.next_sibling;
      exit when (leaf = null);
      i := i + 1;
    end loop;
  end Walk_to_First_Parent;

  procedure Enumerate_Paths_in_Poset ( k,n : in integer32 ) is

  -- DESCRIPTION :
  --   Prompts the user for a configuration of white checkers,
  --   creates a poset and writes all paths through the poset.

    p : constant Vector(1..n) := Identity_Permutation(natural32(n));
    rows,cols : Vector(1..k);
    ps : Poset;
    ans : character;

  begin
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
    Write_Paths_in_Poset(ps);
    put("Walk to first parent ? (y/n) "); Ask_Yes_or_No(ans);
    if ans = 'y'
     then Walk_to_First_Parent(ps);
    end if;
  end Enumerate_Paths_in_Poset;

  procedure Main is

    n,k : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Littlewood-Richardson rule by checker games ...");
    new_line;
    put_line("Choose one of the following : ");
    put_line("  1. test calculation of rising and descending checker");
    put_line("  2. interactive generation of the specialization order");
    put_line("  3. show all permutations of moving black checkers");
    put_line("  4. test calculation of moves of white checkers");
    put_line("  5. create the poset of all black and white moves");
    put_line("  6. count k-planes satisfying intersection conditions");
    put_line("  7. resolve intersection condition given by brackets");
    put_line("  8. enumerate all paths from leaves to root of poset");
    put_line("  9. group black checkers into four zones");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to make a choice : ");
    Ask_Alternative(ans,"123456789");
    new_line;
    put("Give n (dimension of ambient space) : "); get(n);
    case ans is
      when '1' => Interactive_Test(n);
      when '2' => Specialization_Order(n);
      when '3' => Show_All_Moves(natural32(n));
      when '4' | '5' | '6' | '8' | '9' =>
        put("Give k (dimension of subspace) : "); get(k);
        case ans is
          when '4' => Test_White_Moves(k,n);
          when '5' => Create_Poset(k,n);
          when '6' => Count_Planes(k,n);
          when '8' => Enumerate_Paths_in_Poset(k,n);
          when '9' => Group_Black_Checkers(k,n);
          when others => null;
        end case;
      when '7' => Resolve_Intersection_Condition(n);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_checkers;
