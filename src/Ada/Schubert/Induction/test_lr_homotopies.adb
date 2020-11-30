with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Natural_Numbers;           use Multprec_Natural_Numbers;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_VecMats;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Solutions;         use Standard_Complex_Solutions;
with Brackets;                           use Brackets;
with Bracket_Monomials_io;               use Bracket_Monomials_io;
with Checker_Moves;                      use Checker_Moves;
with Checker_Posets,Checker_Posets_io;   use Checker_Posets;
with Intersection_Posets_io;             use Intersection_Posets_io;
with Wrapped_Path_Trackers;
with Setup_Flag_Homotopies;
with Standard_Solution_Posets;           use Standard_Solution_Posets;
with Resolve_Schubert_Problems;          use Resolve_Schubert_Problems;
with Main_Schubert_Induction;

package body Test_LR_Homotopies is

  procedure Read_Brackets ( bm : out Bracket_Monomial ) is
  begin
    new_line;
    put_line("A product of brackets is for example ");
    put_line("  [2 4 6]*[2 4 6]*[2 4 6]; or [2 4 6]^3;");
    put_line("Give a product of brackets, terminate by a semicolon ;");
    get(bm);
    put("Your product : "); put(bm); new_line;
  end Read_Brackets;

  procedure Walk_from_Root_to_Leaves
              ( n,k : in integer32; bm : in Bracket_Monomial ) is

    use Intersection_Posets;

    cnd : constant Array_of_Brackets := Create(bm);
    nbc : constant integer32 := cnd'last;
    ips : constant Intersection_Poset(nbc-1)
        := Main_Schubert_Induction.Process_Conditions(n,k,nbc,cnd);
    tmp : Poset_List;
    lpn : Link_to_Poset_Node;

    procedure Write_Parent ( node : in Link_to_Poset_Node ) is

    -- DESCRIPITION :
    --   Writes the root node of the parent given on input in node
    --   and also writes the leaf nodes of the parent.

    begin
      Checker_Posets_io.Write_Nodes_in_Poset(node.ps,node.ps.black'first);
      Checker_Posets_io.Write_Nodes_in_Poset(node.ps,node.ps.black'last);
    end Write_Parent;
    procedure Write_Parents is
      new Intersection_Posets.Enumerate_Parents(Write_Parent);

  begin
    put_line("Resolving the intersection conditions :");
    Write_Expansion(ips);
    for i in 1..ips.m loop
      put("The nodes at level "); put(i,1); put_line(" :");
      tmp := ips.nodes(i);
      for j in 1..Length_Of(ips.nodes(i)) loop
        lpn := Head_Of(tmp);
        put("-> poset node "); put(j,1); put_line(", root and leaves :");
       -- Checker_Posets_io.Write(lpn.ps);
        Checker_Posets_io.Write_Nodes_in_Poset(lpn.ps,lpn.ps.black'first);
        Checker_Posets_io.Write_Nodes_in_Poset(lpn.ps,lpn.ps.black'last);
        if i > 1 then
          put_line("-> its parents are listed by their root and leaves :");
          Write_Parents(ips.nodes(i-1),lpn.all);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
  end Walk_from_Root_to_Leaves;

  procedure Connect_Checker_Posets_to_Count
              ( pl : in Intersection_Posets.Poset_List;
                nd : in Intersection_Posets.Poset_Node ) is

    use Intersection_Posets;

    procedure Connect_Parent ( node : in Link_to_Poset_Node ) is

    -- DESCRIPTION :
    --   Connects the leaf of the parent given on input in node.

      child : constant Checker_Posets.Poset := nd.ps;
      childnode : constant Link_to_Node := child.white(child.white'first);
      childconds : constant Standard_Natural_Vectors.Vector
                 := childnode.rows;     -- root
      gamenode : Link_to_Node := node.ps.white(node.ps.white'last);
     -- parentconds : constant Standard_Natural_Vectors.Vector
     --             := node.ps.white(node.ps.white'last).cols; -- leaf

    begin
     -- Checker_Posets_io.Write_Nodes_in_Poset(node.ps,node.ps.black'first);
      put_line("** Before assigning coefficients at parent :");
      Checker_Posets_io.Write_Nodes_in_Poset(node.ps,node.ps.black'last);
      put("  conditions at child  : "); put(childconds); new_line;
      put("  conditions at parent : ");
      loop
        put(gamenode.cols); 
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds)
         then Add(gamenode.coeff,childnode.coeff);
        end if;
        exit when (gamenode.next_sibling = null);
        put(" |");
        gamenode := gamenode.next_sibling;
      end loop;
      new_line;
      put_line("** After assigning coefficients at parent :");
      Checker_Posets_io.Write_Nodes_in_Poset(node.ps,node.ps.black'last);
    end Connect_Parent;
    procedure Connect_Parents is
      new Intersection_Posets.Enumerate_Parents(Connect_Parent);

  begin
    Connect_Parents(pl,nd);
  end Connect_Checker_Posets_to_Count;

  procedure Walk_from_Leaves_to_Root
              ( n,k : in integer32; bm : in Bracket_Monomial ) is

    use Intersection_Posets;

    cnd : constant Array_of_Brackets := Create(bm);
    nbc : constant integer32 := cnd'last;
    ips : Intersection_Poset(nbc-1)
        := Main_Schubert_Induction.Process_Conditions(n,k,nbc,cnd);
    tmp : Poset_List;
    lpn : Link_to_Poset_Node;
    roco : Natural_Number;

  begin
   -- put_line("Resolving the intersection conditions :");
   -- Write_Expansion(ips);
    roco := Final_Sum(ips);
    put("The formal root count : "); put(roco); new_line;
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    for i in reverse 1..ips.m loop
      new_line;
      put("The nodes at level "); put(i,1); put_line(" :");
      tmp := ips.nodes(i);
      for j in 1..Length_Of(ips.nodes(i)) loop
        lpn := Head_Of(tmp);
       -- put_line("poset before adding from leaves to root : ");
       -- Checker_Posets_io.Write(lpn.ps);
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        put_line("poset after adding from leaves to root : ");
        Checker_Posets_io.Write(lpn.ps);
        put("-> poset node "); put(j,1); put_line(", root and leaves :");
        Checker_Posets_io.Write_Nodes_in_Poset(lpn.ps,lpn.ps.black'first);
        Checker_Posets_io.Write_Nodes_in_Poset(lpn.ps,lpn.ps.black'last);
        if i > 1 then
         -- put_line("-> its parents are listed by their root :");
          put_line("-> its parents are listed by their leaf :");
          Connect_Checker_Posets_to_Count(ips.nodes(i-1),lpn.all);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    put(" Top down root count : "); put(roco); new_line;
    lpn := Head_Of(ips.nodes(1));
    put("Bottom up root count : ");
    put(lpn.ps.white(lpn.ps.white'first).coeff); new_line;
  end Walk_from_Leaves_to_Root;

  procedure Stubbing_Littlewood_Richardson
              ( n,k : in integer32; bm : in Bracket_Monomial ) is

    use Intersection_Posets;

    cnd : constant Array_of_Brackets := Create(bm);
    nbc : constant integer32 := cnd'last;
    ips : Intersection_Poset(nbc-1)
        := Main_Schubert_Induction.Process_Conditions(n,k,nbc,cnd);
    top_roco,bottom_roco : Natural_Number;

  begin
    top_roco := Final_Sum(ips);
    put("The formal root count : "); put(top_roco); new_line;
    Count_Roots(standard_output,ips,bottom_roco);
    put(" Top down root count : "); put(top_roco); new_line;
    put("Bottom up root count : "); put(bottom_roco); new_line;
  end Stubbing_Littlewood_Richardson;

  procedure Make_Solution_Poset
              ( n,k : in integer32; bm : in Bracket_Monomial ) is

    use Intersection_Posets;

    cnd : constant Array_of_Brackets := Create(bm);
    nbc : constant integer32 := cnd'last;
    ips : constant Intersection_Poset(nbc-1)
        := Main_Schubert_Induction.Process_Conditions(n,k,nbc,cnd);
    sps : Solution_Poset(ips.m) := Create(ips);
    tmp : Solnode_List;
    snd : Link_to_Solution_Node;

  begin
    for i in sps.nodes'range loop
      put("Length of solution lists at level "); put(i,1); put(" :");
      tmp := sps.nodes(i);
      while not Is_Null(tmp) loop
        snd := Head_Of(tmp);
        put(" "); put(Length_Of(snd.sols),1);
        tmp := Tail_Of(tmp);
      end loop;
      new_line;
    end loop;
    put("Clearing the solution poset ...");
    Clear(sps);
    put_line(" done.");
  end Make_Solution_Poset;

  procedure Resolve_Schubert_Problem
              ( n,k : in integer32; bm : in Bracket_Monomial ) is

    use Intersection_Posets;

    file : file_type;
    cnd : constant Array_of_Brackets := Create(bm);
    nbc : constant integer32 := cnd'last;
    ips : Intersection_Poset(nbc-1)
        := Main_Schubert_Induction.Process_Conditions(n,k,nbc,cnd);
    sps : Solution_Poset(ips.m) := Create(ips);
    top_roco,bottom_roco : Natural_Number;
    q : constant Standard_Natural_Vectors.Vector
      := Identity_Permutation(natural32(n));
    rows : constant Standard_Natural_Vectors.Vector
         := Main_Schubert_Induction.Bracket_to_Vector(cnd(cnd'first).all);
    cols : constant Standard_Natural_Vectors.Vector
         := Main_Schubert_Induction.Bracket_to_Vector(cnd(cnd'first+1).all);
    conds : constant Standard_Natural_VecVecs.VecVec(1..nbc-2)
          := Main_Schubert_Induction.Remaining_Intersection_Conditions(cnd);
    link2conds : constant Standard_Natural_VecVecs.Link_to_VecVec
               := new Standard_Natural_VecVecs.VecVec'(conds);
    flags : Standard_Complex_VecMats.VecMat(1..nbc-2);
    fsys : Link_to_Poly_Sys;
    sols : Solution_List;
    tol : constant double_float := 1.0E-6;
    ans : character;
    monitor_games,report,verify,minrep,tosqr : boolean;
    timer : Timing_Widget;

  begin
    new_line;
    put_line("Reading a name for the output file ...");
    Read_Name_and_Create_File(file);
    new_line;
    put("Monitor Littlewood-Richardson homotopies"
      & " in each checker game ? (y/n) ");
    Ask_Yes_or_No(ans);
    new_line;
    monitor_games := (ans = 'y');
    new_line;
    put("Use an efficient problem formulation ? (y/n) ");
    Ask_Yes_or_No(ans);
    minrep := (ans = 'y');
    put("Square the overdetermined homotopies ? (y/n) ");
    Ask_Yes_or_No(ans);
    tosqr := (ans = 'y');
    new_line;
    put("Diagnostic verification needed ? (y/n) ");
    Ask_Yes_or_No(ans);
    verify := (ans = 'y');
    top_roco := Final_Sum(ips);
    put("The formal root count : "); put(top_roco); new_line;
    put_line("... running the root counting from the bottom up ...");
    Count_Roots(file,ips,bottom_roco);
    put(" Top down root count : "); put(top_roco); new_line;
    put("Bottom up root count : "); put(bottom_roco); new_line;
    for i in flags'range loop
      declare
        randflag : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
                 := Setup_Flag_Homotopies.Random_Flag(n);
      begin
        flags(i) := new Standard_Complex_Matrices.Matrix'(randflag);
      end;
    end loop;
   -- put_line("... resolving the Schubert problem ...");
   -- new_line;
   -- put_line("See the output file for results ...");
   -- new_line;
    Wrapped_Path_Trackers.Set_Parameters(file,report);
    tstart(timer);
    Resolve(file,monitor_games,report,n,k,0,tol,ips,sps,
            verify,minrep,tosqr,conds,flags,sols);
    tstop(timer);
    Main_Schubert_Induction.Write_Results
     (file,n,k,q,rows,cols,minrep,link2conds,flags,sols,fsys);
    new_line(file);
    print_times(file,timer,"resolving a Schubert problem");
  end Resolve_Schubert_Problem;

  procedure Main is

    ans : character;
    n,k : integer32 := 0;
    bm : Bracket_Monomial;

  begin
    new_line;
    put_line("MENU to solve Schubert problems with LR homotopies : ");
    put_line("  0. walk through intersection poset from root to leaves;");
    put_line("  1. walk through intersection poset from leaves to root;");
    put_line("  2. stubbing the run of Littlewood-Richardson homotopies;");
    put_line("  3. make an empty intersection solution poset.");
    put_line
      ("  4. resolve Schubert problem via Littlewood-Richardson homotopies;");
    put("Type 0, 1, 2, 3, or 4 to select : ");
    Ask_Alternative(ans,"01234");
    new_line;
    put("Give the ambient dimension : "); get(n);
    put("Give the dimension of the planes : "); get(k);
    Read_Brackets(bm);
    case ans is
      when '0' => Walk_from_Root_to_Leaves(n,k,bm);
      when '1' => Walk_from_Leaves_to_Root(n,k,bm);
      when '2' => Stubbing_Littlewood_Richardson(n,k,bm);
      when '3' => Make_Solution_Poset(n,k,bm);
      when '4' => Resolve_Schubert_Problem(n,k,bm);
      when others => null;
    end case;
  end Main;

end Test_LR_Homotopies;
