with text_io;                          use text_io;
with Standard_Natural_Numbers;         use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;         use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Communications_with_User;         use Communications_with_User;
with Brackets;                         use Brackets;
-- with Brackets_io;                      use Brackets_io;
with Localization_Posets;              use Localization_Posets;
with Localization_Posets_io;           use Localization_Posets_io;
with Make_Input_Planes;

procedure ts_posets is

-- DESCRIPTION :
--   Test on the construction of localization posets.

  --function Determine_Root ( m,p : natural32 ) return Node is

  -- DESCRIPTION :
  --   Proposes the trivial root to the user, allowing the user to
  --   modify this choice.

  --  root : Node(integer32(p)) := Trivial_Root(m,p);
  --  ans : character;

  --begin
  --  loop
  --    put("Top and bottom pivots of root are ");
  --    put(root.top); put(" and ");
  --    put(root.bottom); put_line(".");
  --    put("Level of the root : "); put(root.level,1); new_line;
  --    put("Do you want to use another root ? (y/n) "); get(ans);
  --    exit when (ans /= 'y');
  --    put("Give top pivots : "); get(root.top);
  --    put("Give bottom pivots : "); get(root.bottom);
  --    put("Give level of root : "); get(root.level);
  --  end loop;
  --  return root;
  --end Determine_Root;

  procedure Write_Poset
                ( file : in file_type;
                  lnkroot : in Link_to_Node; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates the posets and writes them onto the file.

    nq : constant integer32 := integer32(m*p + q*(m+p));
    level_poset : Array_of_Nodes(0..nq);
    index_poset : Array_of_Array_of_Nodes(0..nq);
    nbp : natural32;

  begin
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    index_poset := Create_Indexed_Poset(level_poset);
    put(file,index_poset);
    nbp := Root_Count_Sum(level_poset);
    put(file,"The number of paths : "); put(file,nbp,1); new_line(file);
  end Write_Poset;

  procedure Create_Top_Hypersurface_Poset ( m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing only top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Top_Create(lnkroot,m+p);
    put_line("The poset created from the top : ");
    Write_Poset(Standard_Output,lnkroot,m,p,0);
  end Create_Top_Hypersurface_Poset;

  procedure Create_Top_Hypersurface_Poset ( m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing only top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : Link_to_Node := new Node'(root);

  begin
    Q_Top_Create(lnkroot,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the top : ");
    Write_Poset(Standard_Output,lnkroot,m,p,q);
  end Create_Top_Hypersurface_Poset;

  procedure Create_Bottom_Hypersurface_Poset ( m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by decrementing only bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Bottom_Create(lnkroot);
    put_line("The poset created from the bottom : ");
    Write_Poset(Standard_Output,lnkroot,m,p,0);
  end Create_Bottom_Hypersurface_Poset;

  procedure Create_Bottom_Hypersurface_Poset ( m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by decrementing only bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Q_Bottom_Create(lnkroot,m+p);
    put_line("The poset created from the bottom : ");
    Write_Poset(Standard_Output,lnkroot,m,p,q);
  end Create_Bottom_Hypersurface_Poset;

  procedure Create_Mixed_Hypersurface_Poset ( m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing top and decrementing bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Top_Bottom_Create(lnkroot,m+p);
    put_line("The poset created in a mixed fashion : ");
    Write_Poset(Standard_Output,lnkroot,m,p,0);
  end Create_Mixed_Hypersurface_Poset;

  procedure Create_Mixed_Hypersurface_Poset ( m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Create the poset by incrementing top and decrementing bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);

  begin
    Q_Top_Bottom_Create(lnkroot,root.bottom(integer32(p)),m+p);
    put_line("The poset created in a mixed fashion : ");
    Write_Poset(Standard_Output,lnkroot,m,p,q);
  end Create_Mixed_Hypersurface_Poset;

  procedure Create_Top_General_Poset ( m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);

  begin
    Top_Create(lnkroot,codim,m+p);
    put_line("The poset created from the top : ");
    Write_Poset(Standard_Output,lnkroot,m,p,0);
  end Create_Top_General_Poset;

  procedure Create_Bottom_General_Poset ( m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);

  begin
    Bottom_Create(lnkroot,codim);
    put_line("The poset created from the bottom : ");
    Write_Poset(Standard_Output,lnkroot,m,p,0);
  end Create_Bottom_General_Poset;

  procedure Create_Mixed_General_Poset ( m,p : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by incrementing the top and decrementing the bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,0);

  begin
    Top_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created in a mixed fashion : ");
    Write_Poset(Standard_Output,lnkroot,m,p,0);
  end Create_Mixed_General_Poset;

  procedure Create_Top_General_Poset ( m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    Q_Top_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created from the top : ");
    Write_Poset(Standard_Output,lnkroot,m,p,q);
  end Create_Top_General_Poset;

  procedure Create_Bottom_General_Poset ( m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by consistently incrementing the top pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    Q_Bottom_Create(lnkroot,codim,m+p);
    put_line("The poset created from the bottom : ");
    Write_Poset(Standard_Output,lnkroot,m,p,q);
  end Create_Bottom_General_Poset;

  procedure Create_Mixed_General_Poset ( m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Creates a poset for counting general subspace intersections,
  --   by incrementing the top and decrementing the bottom pivots.

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    codim : constant Bracket := Make_Input_Planes.Read_Codimensions(m,p,q);

  begin
    Q_Top_Bottom_Create(lnkroot,codim,root.bottom(integer32(p)),m+p);
    put_line("The poset created in a mixed fashion : ");
    Write_Poset(Standard_Output,lnkroot,m,p,q);
  end Create_Mixed_General_Poset;

  procedure Test_Root_Counts
              ( file : in file_type;
                m,p,q : in natural32; codim : in Bracket;
                bug : out boolean ) is

  -- DESCRIPTION :
  --   Computes the root count in various ways for the given vector
  --   of co-dimensions.  Compares the results and reports bugs.

    mpq : constant integer32 := integer32(m*p + q*(m+p));
    top_root0,bot_root0,mix_root0 : Node(integer32(p));
    lnk_top_root0 : Link_to_Node := new Node'(top_root0);
    lnk_bot_root0 : Link_to_Node := new Node'(bot_root0);
    lnk_mix_root0 : Link_to_Node := new Node'(mix_root0);
    top_poset0,bot_poset0,mix_poset0 : Array_of_Nodes(0..mpq);
    top_rootq,bot_rootq,mix_rootq : Node(integer32(p));
    lnk_top_rootq : Link_to_Node := new Node'(top_rootq);
    lnk_bot_rootq : Link_to_Node := new Node'(bot_rootq);
    lnk_mix_rootq : Link_to_Node := new Node'(mix_rootq);
    top_posetq,bot_posetq,mix_posetq : Array_of_Nodes(0..mpq);

  begin
    bug := false;
    if q = 0 then
      top_root0 := Trivial_Root(m,p);
      bot_root0 := Trivial_Root(m,p);
      mix_root0 := Trivial_Root(m,p);
      lnk_top_root0 := new Node'(top_root0);
      lnk_bot_root0 := new Node'(bot_root0);
      lnk_mix_root0 := new Node'(mix_root0);
      Top_Create(lnk_top_root0,codim,m+p);
      Bottom_Create(lnk_bot_root0,codim);
      Top_Bottom_Create(lnk_mix_root0,codim,m+p);
      top_poset0 := Create_Leveled_Poset(lnk_top_root0);
      bot_poset0 := Create_Leveled_Poset(lnk_bot_root0);
      mix_poset0 := Create_Leveled_Poset(lnk_mix_root0);
      Count_Roots(top_poset0);
      Count_Roots(bot_poset0);
      Count_Roots(mix_poset0);
    end if;
    top_rootq := Trivial_Root(m,p,q);
    bot_rootq := Trivial_Root(m,p,q);
    mix_rootq := Trivial_Root(m,p,q);
    lnk_top_rootq := new Node'(top_rootq);
    lnk_bot_rootq := new Node'(bot_rootq);
    lnk_mix_rootq := new Node'(mix_rootq);
    Q_Top_Create(lnk_top_rootq,codim,top_rootq.bottom(integer32(p)),m+p);
    Q_Bottom_Create(lnk_bot_rootq,codim,m+p);
    Q_Top_Bottom_Create(lnk_mix_rootq,codim,mix_rootq.bottom(integer32(p)),m+p);
    top_posetq := Create_Leveled_Poset(lnk_top_rootq);
    bot_posetq := Create_Leveled_Poset(lnk_bot_rootq);
    mix_posetq := Create_Leveled_Poset(lnk_mix_rootq);
    Count_Roots(top_posetq); 
    Count_Roots(bot_posetq);
    Count_Roots(mix_posetq);
    if q = 0 then
      put(file,top_poset0(mpq).roco,1);
      if top_poset0(mpq).roco = bot_poset0(mpq).roco then
        put(file," = ");
        put(file,bot_poset0(mpq).roco,1); bug := false;
        if bot_poset0(mpq).roco = mix_poset0(mpq).roco then
          bug := false;
          put(file," = "); put(file,mix_poset0(mpq).roco,1);
        else
          bug := true;
          put(file," <> "); put(file,mix_poset0(mpq).roco,1);
          put_line(file,"  BUG !!!");
          put_line(file,"The poset created incrementing top pivots : ");
          Write_Poset(file,lnk_top_root0,m,p,q);
          put_line(file,"The poset created decrementing bottom pivots : ");
          Write_Poset(file,lnk_bot_root0,m,p,q);
          put_line(file,"The poset created in a mixed fashion : ");
          Write_Poset(file,lnk_mix_root0,m,p,q);
        end if;
      else
        bug := true;
        put(file," <> "); put(file,bot_poset0(mpq).roco,1);
        put_line(file,"  BUG !!!");
        put_line(file,"The poset created incrementing top pivots : ");
        Write_Poset(file,lnk_top_root0,m,p,q);
        put_line(file,"The poset created decrementing bottom pivots : ");
        Write_Poset(file,lnk_bot_root0,m,p,q);
      end if;
    end if;
    if q = 0 then
      if top_posetq(mpq).roco /= top_poset0(mpq).roco then
        bug := true;
        put(file," <> "); put(file,top_posetq(mpq).roco,1);
        put_line(file,"  BUG !!!");
        put_line(file,"The poset created without q = 0 : ");
        Write_Poset(file,lnk_top_root0,m,p,q);
        put_line(file,"The poset created with q = 0 : ");
        Write_Poset(file,lnk_bot_rootq,m,p,q);
      else
        put(file," = ");
      end if;
    end if;
    if not bug then
      put(file,top_posetq(mpq).roco,1);
      if top_posetq(mpq).roco = bot_posetq(mpq).roco then
        put(file," = ");
        put(file,bot_posetq(mpq).roco,1); bug := false;
        if bot_posetq(mpq).roco = mix_posetq(mpq).roco then
          bug := false;
          put(file," = ");
          put(file,mix_posetq(mpq).roco,1); new_line(file);
        else
          bug := true;
          put(file," <> "); put(file,mix_posetq(mpq).roco,1);
          put_line(file,"  BUG !!!");
          put_line(file,"The poset created incrementing top pivots : ");
          Write_Poset(file,lnk_top_rootq,m,p,q);
          put_line(file,"The poset created decrementing bottom pivots : ");
          Write_Poset(file,lnk_bot_rootq,m,p,q);
          put_line(file,"The poset created in a mixed fashion : ");
          Write_Poset(file,lnk_mix_rootq,m,p,q);
        end if;
      else
        bug := true;
        put(file," <> "); put(file,bot_posetq(mpq).roco,1);
        put_line(file,"  BUG !!!");
        put_line(file,"The poset created incrementing top pivots : ");
        Write_Poset(file,lnk_top_rootq,m,p,q);
        put_line(file,"The poset created decrementing bottom pivots : ");
        Write_Poset(file,lnk_bot_rootq,m,p,q);
      end if;
    end if;
    Clear(top_poset0); Clear(bot_poset0); Clear(mix_poset0);
    Clear(top_posetq); Clear(bot_posetq); Clear(mix_posetq);
  end Test_Root_Counts;

  procedure Enumerate_Partitions
               ( file : in file_type; m,p,q : in natural32 ) is

  -- DESCRIPTION :
  --   Test the root counts for all partitions of m*p + q*(m+p).
  --   The results are written on file.

    n : constant integer32 := integer32(m*p + q*(m+p));
    accu : Bracket(1..n);
    bug : boolean := false;

    procedure Enumerate ( k,nk : in integer32 ) is
    begin
      if nk = 0 then
        put(file,n,1); put(file," = ");
        for i in 1..k-2 loop
          put(file,accu(i),1); put(file," + ");
        end loop;
        put(file,accu(k-1),1); put(file," : ");
        Test_Root_Counts(file,m,p,q,accu(1..k-1),bug);
      else
        for i in 1..nk loop
          exit when (i > integer32(m));
          accu(k) := natural32(i);
          Enumerate(k+1,nk-i);
          exit when bug;
        end loop;
      end if;
    end Enumerate;

  begin
    Enumerate(1,n);
  end Enumerate_Partitions;

  procedure Main is

    m,p,q : natural32 := 0;
    ans : character;
    file : file_type;

  begin
    loop
      new_line;
      put_line("MENU for posets for counting p-planes in (m+p)-space : ");
      put_line("  0. exit this program.");
      put_line("-------- the case q = 0 ------------------------------------");
      put_line("  1. k_i == 1 consistently incrementing the top pivots.");
      put_line("  2.          consistently decrementing the bottom pivots.");
      put_line("  3.          mixed top-bottom sequence for poset creation.");
      put_line("  4. k_i >= 1 consistently incrementing the top pivots.");
      put_line("  5.          consistently decrementing the bottom pivots.");
      put_line("  6.          mixed top-bottom sequence for poset creation.");
      put_line("  7. Enumerate all partitions of m*p and test root counts.");
      put_line("-------- the case q >= 0 -----------------------------------");
      put_line("  8. k_i == 1 consistently incrementing top pivots.");
      put_line("  9.          consistently decrementing bottom pivots.");
      put_line("  A.          mixed top-bottom sequence for pivots.");
      put_line("  B. k_i >= 1 consistently incrementing top pivots.");
      put_line("  C.          consistently decrementing bottom pivots.");
      put_line("  D.          mixed top-bottom sequence for pivots.");
      put_line("  E. Test root counts for all partitions of m*p + q*(m+p).");
      put_line("------------------------------------------------------------");
      put("Type 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, or E to choose : ");
      Ask_Alternative(ans,"0123456789ABCDE");
      exit when ans = '0';
      if ans = '7' or ans = 'E' then
        new_line;
        put_line("Reading the name of the output file.");
        Read_Name_and_Create_File(file);
      end if;
      new_line;
      put("Give p, the number of entries in bracket : "); get(p);
      put("Give m, the complementary dimension : "); get(m);
      new_line;
      case ans is
        when '1' => Create_Top_Hypersurface_Poset(m,p);
        when '2' => Create_Bottom_Hypersurface_Poset(m,p);
        when '3' => Create_Mixed_Hypersurface_Poset(m,p);
        when '4' => Create_Top_General_Poset(m,p);
        when '5' => Create_Bottom_General_Poset(m,p);
        when '6' => Create_Mixed_General_Poset(m,p);
        when '7' => Enumerate_Partitions(file,m,p,0);
        when '8' => put("Give q, the degree of the maps : "); get(q);
                    Create_Top_Hypersurface_Poset(m,p,q);
        when '9' => put("Give q, the degree of the maps : "); get(q);
                    Create_Bottom_Hypersurface_Poset(m,p,q);
        when 'A' => put("Give q, the degree of the maps : "); get(q);
                    Create_Mixed_Hypersurface_Poset(m,p,q);
        when 'B' => put("Give q, the degree of the maps : "); get(q);
                    Create_Top_General_Poset(m,p,q);
        when 'C' => put("Give q, the degree of the maps : "); get(q);
                    Create_Bottom_General_Poset(m,p,q);
        when 'D' => put("Give q, the degree of the maps : "); get(q);
                    Create_Mixed_General_Poset(m,p,q);
        when 'E' => put("Give q, the degree of the maps : "); get(q);
                    Enumerate_Partitions(file,m,p,q);
        when others => put_line("Option not recognized.  Please try again.");
      end case;
    end loop;
  end Main;

begin
  new_line;
  put_line("Test on localization posets for linear subspace intersections.");
  Main;
end ts_posets;
