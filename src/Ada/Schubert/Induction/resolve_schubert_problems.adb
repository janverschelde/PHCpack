with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Standard_Natural_Vectors;
with Checker_Posets,Checker_Posets_io;   use Checker_Posets_io;

package body Resolve_Schubert_Problems is

  procedure Initialize_Leaves ( pl : in out Poset_List ) is

    tmp : Intersection_Posets.Poset_List := pl;
    lpn : Intersection_Posets.Link_to_Poset_Node;
    cps : Checker_Posets.Poset;

  begin
    while not Is_Null(tmp) loop
      lpn := Head_Of(tmp);
      cps := lpn.ps;
      Checker_Posets.Set_Coefficients_to_Zero(cps);
      Clear(cps.white(cps.white'last).coeff);
      cps.white(cps.white'last).coeff := Create(natural32(1));
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Leaves;

  procedure Initialize_Nodes ( pl : in out Poset_List ) is

    tmp : Intersection_Posets.Poset_List := pl;
    lpn : Intersection_Posets.Link_to_Poset_Node;
    cps : Checker_Posets.Poset;

  begin
    while not Is_Null(tmp) loop
      lpn := Head_Of(tmp);
      cps := lpn.ps;
      Checker_Posets.Set_Coefficients_to_Zero(cps);
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Nodes;

  procedure Connect_Checker_Posets_to_Count
              ( file : in file_type;
                pl : in Poset_List; nd : in Poset_Node ) is

    procedure Connect_Parent ( node : in Link_to_Poset_Node ) is

    -- DESCRIPTION :
    --   Connects the leaf of the parent given on input in node.

      child : constant Checker_Posets.Poset := nd.ps;
      childnode : constant Checker_Posets.Link_to_Node
                := child.white(child.white'first);
      childconds : constant Standard_Natural_Vectors.Vector
                 := childnode.rows;     -- root
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
      parentconds : constant Standard_Natural_Vectors.Vector
                  := node.ps.white(node.ps.white'last).cols; -- leaf

      use Checker_Posets;

    begin
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          put(file,"*** number of paths from child to the parent : ");
          put(file,childnode.coeff); put_line(file," ***");
        end if;
        exit when (gamenode.next_sibling = null);
        gamenode := gamenode.next_sibling;
      end loop;
      new_line(file);
      put_line(file,"** After assigning coefficients at parent :");
      Write_Nodes_in_Poset(file,node.ps,node.ps.black'last);
    end Connect_Parent;
    procedure Connect_Parents is
      new Intersection_Posets.Enumerate_Parents(Connect_Parent);

  begin
    Connect_Parents(pl,nd);
  end Connect_Checker_Posets_to_Count;

  procedure Connect_Checker_Posets_to_Track
              ( file : in file_type; level : in integer32;
                pl : in Poset_List; snd : in Link_to_Solution_Node;
                sps : in out Solution_Poset;
                flags : in out Standard_Complex_VecMats.VecMat ) is

    nd : constant Link_to_Poset_Node := snd.lpnd;

    procedure Connect_Parent ( node : in Link_to_Poset_Node ) is

    -- DESCRIPTION :
    --   Connects the leaf of the parent given on input in node.

      child : constant Checker_Posets.Poset := nd.ps;
      parent : constant Checker_Posets.Poset := node.ps;
      parent_snd : constant Link_to_Solution_Node
                 := Retrieve(sps.nodes(level),
                             Checker_Posets.Root_Rows(parent),
                             Checker_Posets.Root_Columns(parent));
      childnode : constant Checker_Posets.Link_to_Node
                := child.white(child.white'first);
      childconds : constant Standard_Natural_Vectors.Vector
                 := childnode.rows;     -- root
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
      parentconds : constant Standard_Natural_Vectors.Vector
                  := node.ps.white(node.ps.white'last).cols; -- leaf

      use Checker_Posets;

    begin
      if parent_snd /= null
       then put_line(file,"Found node with start solutions at parent.");
       else put_line(file,"No node with start solutions at parent found!");
      end if;
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          put(file,"*** number of paths from child to the parent : ");
          put(file,childnode.coeff); put_line(file," ***");
        end if;
        exit when (gamenode.next_sibling = null);
        gamenode := gamenode.next_sibling;
      end loop;
      new_line(file);
      put_line(file,"** After assigning coefficients at parent :");
      Write_Nodes_in_Poset(file,node.ps,node.ps.black'last);
    end Connect_Parent;
    procedure Connect_Parents is
      new Intersection_Posets.Enumerate_Parents(Connect_Parent);

  begin
    Connect_Parents(pl,nd.all);
  end Connect_Checker_Posets_to_Track;

  procedure Count_Roots
              ( file : in file_type; ips : in out Intersection_Poset;
                roco : out Natural_Number ) is
   
    tmp : Poset_List;
    lpn : Link_to_Poset_Node;

  begin
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    for i in reverse 1..ips.m loop
      new_line(file);
      put(file,"Solving at level "); put(file,i,1); put_line(file," :");
      tmp := ips.nodes(i);
      for j in 1..Length_Of(ips.nodes(i)) loop
        lpn := Head_Of(tmp);
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        put(file,"-> poset node "); put(file,j,1);
        put_line(file,", root and leaves :");
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'first);
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'last);
        put(file,"*** number of paths tracking in checker game : ");
        put(file,lpn.ps.white(lpn.ps.white'first).coeff);
        put_line(file," ***");
        if i > 1 then
          put_line(file,"-> solving at the leaves of its parents :");
          Connect_Checker_Posets_to_Count(file,ips.nodes(i-1),lpn.all);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    Copy(lpn.ps.white(lpn.ps.white'first).coeff,roco);
  end Count_Roots;

  procedure Resolve
              ( file : in file_type;
                ips : in out Intersection_Poset;
                sps : in out Solution_Poset;
                flags : in out Standard_Complex_VecMats.VecMat;
                sols : out Solution_List ) is

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;

  begin
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    for i in reverse 1..sps.m loop
      new_line(file);
      put(file,"Solving at level "); put(file,i,1); put_line(file," :");
      tmp := sps.nodes(i);
      for j in 1..Length_Of(sps.nodes(i)) loop
        snd := Head_Of(tmp);
        lpn := snd.lpnd;
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        put(file,"-> poset node "); put(file,j,1);
        put_line(file,", root and leaves :");
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'first);
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'last);
        put(file,"*** number of paths tracking in checker game : ");
        put(file,lpn.ps.white(lpn.ps.white'first).coeff);
        put_line(file," ***");
        if i > 1 then
          put_line(file,"-> solving at the leaves of its parents :");
          Connect_Checker_Posets_to_Track
            (file,i-1,ips.nodes(i-1),snd,sps,flags);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      sps.level := i; -- note that level of completion goes in reverse!
    end loop;
    put(file,"The formal root count : ");
    put(file,lpn.ps.white(lpn.ps.white'first).coeff);
    new_line(file);
    snd := Head_Of(sps.nodes(sps.level));
    sols := snd.sols;
  end Resolve;
 
end Resolve_Schubert_Problems;
