with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Localization_Posets;               use Localization_Posets;
with Localization_Posets_io;            use Localization_Posets_io;
with Localization_Poset_Strings;        use Localization_Poset_Strings;
with Pieri_Root_Count;

procedure ts_piroco is

-- DESCRIPTION :
--   This interactive program reads the (m,p,q) input parameters and
--   invokes the combinatorial root counting based on Pieri's rule.

  procedure Write_Nodes_Strings ( p : in Array_of_Nodes ) is

    np : natural32;

  begin
    if p'last >= 10
     then np := 2;
     else np := 1;
    end if;
    for n in p'range loop
      put("n = "); put(n,np); put(" : ");
      if p(n) = null
       then new_line;
       else put_line(Nodes_to_String(p(n).all));
      end if;
    end loop;
  end Write_Nodes_Strings;

  procedure Show_the_Poset ( m,p,q : in natural32 ) is

    root : constant Node(integer32(p)) := Trivial_Root(m,p,q);
    lnkroot : constant Link_to_Node := new Node'(root);
    nq : constant natural32 := m*p + q*(m+p);
    level_poset : Array_of_Nodes(0..integer32(nq));
   -- index_poset : Array_of_Array_of_Nodes(0..integer32(nq));

  begin
    Q_Top_Bottom_Create(lnkroot,root.bottom(integer32(p)),m+p);
    level_poset := Create_Leveled_Poset(lnkroot);
    Count_Roots(level_poset);
    new_line;
    put("root node = ");
    put_line(Node_to_String(root.top,root.bottom,natural32(root.roco)));
    new_line;
    put_line(Poset_to_String(level_poset));
    new_line;
    Write_Nodes_Strings(level_poset);
    new_line;
    put_line("the root count poset :"); put(level_poset);
   -- index_poset := Create_Indexed_Poset(level_poset);
   -- new_line;
   -- put_line("the poset with links to children:");
   -- put(index_poset);
  end Show_the_Poset;

  procedure Main is

    m,p,q,n : natural32 := 0;

  begin
    new_line;
    put_line("Combinatorial root count in the quantum Schubert calculus.");
    new_line;
    put("Give the dimension of the input planes : "); get(m);
    put("Give the dimension of the output planes : "); get(p);
    put("Give the degree of the maps : "); get(q);
    new_line;
    n := Pieri_Root_Count(m,p,q);
    put("The number of degree "); put(q,1);
    put(" maps producing "); put(p,1); put_line("-planes");
    put("that meet "); put(m*p,1); put(" general ");
    put(m,1); put("-planes is "); put(n,1); put_line(".");
    Show_the_Poset(m,p,q);
  end Main;

begin
  Main;
end ts_piroco;
