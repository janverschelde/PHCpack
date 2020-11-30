with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Multprec_Natural_Numbers_io;       use Multprec_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;          use Standard_Natural_VecVecs;
with Standard_Natural_VecVecs_io;       use Standard_Natural_VecVecs_io;
with Standard_Natural_Matrices;         use Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;      use Standard_Natural_Matrices_io;
with Brackets_io;                       use Brackets_io;
with Bracket_Monomials_io;              use Bracket_Monomials_io;
with Symbolic_Schubert_Conditions;
with Checker_Boards_io;                 use Checker_Boards_io;
with Checker_Moves;                     use Checker_Moves;
with Checker_Posets_io;                 use Checker_Posets_io;
with Checker_Localization_Patterns;     use Checker_Localization_Patterns;
with Intersection_Posets;               use Intersection_Posets;
with Intersection_Posets_io;            use Intersection_Posets_io;
with Main_Schubert_Induction;

package body Test_Checker_Games is

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
        Intersect(ips,w,false); -- be not silent
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
    fs : Natural_Number;
    ans : character;

    use Main_Schubert_Induction;

  begin
    put("Do you want intermediate output ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Create_Intersection_Poset(n,nb,cd,false,fs);
     else Create_Intersection_Poset(n,nb,cd,true,fs);
    end if;
    put("The number of solutions : "); put(fs); new_line;
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

    procedure Write_Paths is
      new Checker_Posets.Enumerate_Paths_in_Poset(Write_Path);

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

  procedure Flip ( b : in out Bracket ) is

    tmp : natural32;

  begin
    for i in b'first..(b'last/2) loop
      tmp := b(i);
      b(i) := b(b'last - i + 1);
      b(b'last - i + 1) := tmp;
    end loop;
  end Flip;

  function Partition_to_Bracket
             ( k,n : integer32; p : Bracket ) return Bracket is

    res : Bracket(1..k);

  begin
    for j in 1..k loop 
      res(j) := natural32(n - k + j - integer32(p(j))); 
    end loop;
    return res;
  end Partition_to_Bracket;

  function Random_Partition ( k,n : integer32 ) return Bracket is

    res : Bracket(1..k);
    bnd : integer32;

  begin
    res(1) := natural32(Standard_Random_Numbers.Random(1,n-k-1));
    for i in 2..k loop
      bnd := integer32(res(i-1));
      res(i) := natural32(Standard_Random_Numbers.Random(0,bnd));
    end loop;
    return res;
  end Random_Partition;

  function Random_Complement ( k,n : integer32; p : Bracket ) return Bracket is

    res : Bracket(1..k);
    low,upp : integer32;

  begin
    upp := n - k - integer32(p(p'first)) - 1;
    res(1) := natural32(Standard_Random_Numbers.Random(0,upp));
    for i in 2..k loop
      low := integer32(res(i-1));
      upp := n - k - integer32(p(i)) - 1;
      res(i) := natural32(Standard_Random_Numbers.Random(low,upp));
    end loop;
    if res(k) = 0
     then res(k) := 1;
    end if;
    return res;
  end Random_Complement;

  function Difference ( k,n : integer32; p,q : Bracket ) return Bracket is

    res : Bracket(1..k);
    val : integer32;

  begin
    for i in 1..k loop
      val := n - k - integer32(p(i)) - integer32(q(i));
      res(i) := natural32(val);
    end loop;
    return res;
  end Difference;

  procedure Sort_into_Partition ( b : in out Bracket ) is

    max : natural32;

  begin
    for i in b'range loop
      max := b(i);
      for j in i+1..b'last loop
        if b(j) > max then
          b(i) := b(j);
          b(j) := max;
          max := b(i);
        end if;
      end loop;
    end loop;
  end Sort_into_Partition;

  function Random_Redistribute
             ( k,n : integer32; b : Bracket ) return Bracket is

    res : Bracket(b'range);
    s : integer32 := integer32(Sum(b));
    upp : integer32;

  begin
    for i in 1..(b'last-1) loop
      if s > n - k 
       then upp := n-k;
       else upp := s;
      end if;
      res(i) := natural32(Standard_Random_Numbers.Random(0,upp));
      s := s - integer32(res(i));
    end loop;
    res(b'last) := natural32(s);
    Sort_into_Partition(res);
    return res;
  end Random_Redistribute;

  procedure Random_Triplet ( k,n : in integer32; roco : out Natural_Number ) is

    p,q,r,bp,bq,br : Bracket(1..k);
    cd : Array_of_Brackets(1..3);

    use Main_Schubert_Induction;
 
  begin
    p := Random_Partition(k,n);
    q := Random_Complement(k,n,p);
    r := Difference(k,n,p,q);
    Flip(q);
    bp := Partition_to_Bracket(k,n,p);
    bq := Partition_to_Bracket(k,n,q);
    cd(1) := new Bracket'(bp);
    cd(2) := new Bracket'(bq);
    Sort_into_Partition(r);
    loop
      r := Random_Redistribute(k,n,r);
      put("The three partitions : ");
      put(p); put(q); put(r); new_line;
      br := Partition_to_Bracket(k,n,r);
      put("The three brackets : ");
      put(bp); put(bq); put(br); new_line;
      cd(3) := new Bracket'(br);
      Create_Intersection_Poset(n,3,cd,false,roco);
      exit when (roco > 0);
      Clear(cd(3));
    end loop;
  end Random_Triplet;

  procedure Generate_Intersection_Problem
              ( k,n,m : in integer32; roco : out Natural_Number ) is

    cd : Array_of_Brackets(1..m);
    choice : Standard_Natural_Vectors.Vector(1..k);
    invpart,cond : Bracket(1..k);
    upper : integer32 := n-k;
    restsum : integer32 := k*(n-k) - 1; -- degree of freedom minus one
    sign : integer32;
    cnt : integer32 := 0;

    use Main_Schubert_Induction;

  begin
    for i in 1..m-1 loop
      for j in 1..k loop
        if upper > restsum
         then upper := restsum;
        end if;
        exit when (upper = 0);
        if j = 1  -- prevent from generating zero condition
         then choice(j) := natural32(Standard_Random_Numbers.Random(1,upper));
         else choice(j) := natural32(Standard_Random_Numbers.Random(0,upper));
        end if;
        restsum := restsum - integer32(choice(j));
      end loop;
      Create(choice,invpart,sign);
      Flip(invpart);
      put("Partition : "); put(invpart);
      cond := Partition_to_Bracket(k,n,invpart);
      put(" -> condition : "); put(cond); new_line;
      cnt := cnt + 1;
      cd(cnt) := new Bracket'(cond);
      exit when (restsum = 0);
    end loop;
    restsum := restsum + 1;  -- at least one degree of freedom left
   -- put("restsum : "); put(restsum,1); new_line;
    invpart := (1..k => 0);
    for j in 1..k loop
      if restsum > n-k 
       then invpart(j) := natural32(n-k);
       else invpart(j) := natural32(restsum);
      end if;
      restsum := restsum - integer32(invpart(j));
      exit when (restsum = 0);
    end loop;
    put("Partition : "); put(invpart);
    cond := Partition_to_Bracket(k,n,invpart);
    put(" -> condition : "); put(cond); new_line;
    cnt := cnt + 1;
    cd(cnt) := new Bracket'(cond);
    Create_Intersection_Poset(n,cnt,cd(1..cnt),false,roco);
  end Generate_Intersection_Problem;

  procedure Random_Intersection_Problem ( k,n : integer32 ) is

    m : integer32 := 0;
    roco : Natural_Number;
    ans : character;

  begin
    loop
      loop
        put("Give the number of conditions (>= 3) : "); get(m);
        exit when (m >= 3);
       put_line("-> number must be at least 3.  Please try again ...");
      end loop;
      if m = 3
       then Random_Triplet(k,n,roco);
       else Generate_Intersection_Problem(k,n,m,roco);
      end if;
      new_line;
      put("The root count : "); put(roco); new_line;
      put("Generate another intersection problem ? (y/n) ");
      Ask_Yes_or_No(ans);
      exit when (ans /= 'y');
      Clear(roco);
    end loop;
  end Random_Intersection_Problem;

  procedure NotAbove ( k,n : in integer32 ) is

    use Symbolic_Schubert_Conditions;

    alpha,beta : Bracket(1..k);
    ans : character;
    cnt : natural32 := 0;

    procedure Write_and_Count ( b : in Bracket; c : out boolean ) is

    -- DESCRIPTION :
    --   Increases the global counter cnt, prints the counter
    --   and the new bracket in the instantiation of Process.

    begin
      cnt := cnt + 1;
      put(cnt,2); put(" : "); put(b); new_line;
      c := true;
    end Write_and_Count;
    procedure Enum is new Enumerate_NotAbove(Write_and_Count);

  begin
    put("Reading two ");
    put(k,1); put(" brackets of ");
    put(n,1); put_line(" elements ...");
    put("Give alpha : "); get(alpha);
    put("Give beta : "); get(beta);
    put(alpha); put(" <= "); put(beta);
    if alpha <= beta
     then put_line(" true");
     else put_line(" false");
    end if;
    new_line;
    put("Number of bracket not above : ");
    put(Number_of_NotAbove(natural32(n),alpha),1);
    put_line(".");
    put("Enumerate all brackets not above ");
    put(alpha); put(" ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Enum(natural32(n),alpha);
    end if;
  end NotAbove;

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
    put_line("  A. generate a random intersection problem");
    put_line("  B. check for brackets not above a given bracket");
    put("Type 1, 2, 3, 4, 5, 6, 7, 8, 9, A, or B to make a choice : ");
    Ask_Alternative(ans,"123456789AB");
    new_line;
    put("Give n (dimension of ambient space) : "); get(n);
    case ans is
      when '1' => Interactive_Test(n);
      when '2' => Specialization_Order(n);
      when '3' => Show_All_Moves(natural32(n));
      when '4' | '5' | '6' | '8' | '9' | 'A' | 'B' =>
        put("Give k (dimension of subspace) : "); get(k);
        case ans is
          when '4' => Test_White_Moves(k,n);
          when '5' => Create_Poset(k,n);
          when '6' => Count_Planes(k,n);
          when '8' => Enumerate_Paths_in_Poset(k,n);
          when '9' => Group_Black_Checkers(k,n);
          when 'A' => Random_Intersection_Problem(k,n);
          when 'B' => NotAbove(k,n);
          when others => null;
        end case;
      when '7' => Resolve_Intersection_Condition(n);
      when others => null;
    end case;
  end Main;

end Test_Checker_Games;
