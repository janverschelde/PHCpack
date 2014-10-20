with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Checker_Posets,Checker_Posets_io;   use Checker_Posets_io;
with Checker_Localization_Patterns;
with Moving_Flag_Homotopies;             use Moving_Flag_Homotopies;
with Moving_Flag_Continuation;

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

  function Flip ( b : Standard_Natural_Vectors.Vector )
                return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector b with its entries in reverse order.

    res : Standard_Natural_Vectors.Vector(b'range);

  begin
    for i in b'range loop
      res(i) := b(b'last-i+1);
    end loop;
    return res;
  end Flip;

  procedure Start_Solution 
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                snd : in out Link_to_Solution_Node;
                fail : out boolean;
                res : out double_float ) is

    slnp : constant Checker_Posets.Poset := snd.lpnd.ps;
    rows : Standard_Natural_Vectors.Vector
         := Checker_Posets.Root_Rows(slnp);
    cols : Standard_Natural_Vectors.Vector(rows'range) := Flip(rows);
    q : constant Standard_Natural_Vectors.Vector
      := slnp.black(slnp.black'last).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    x : Standard_Complex_Vectors.Vector(1..integer32(dim));
    eqs : Link_to_Poly_Sys;

  begin
    put_line(file,"The localization map : "); put(file,locmap);
    put(file,"The number of variables : ");
    put(file,dim,1); put_line(file,".");
    put(file,"q = "); put(file,q);
    put(file,"  rows = "); put(file,rows);
    put(file,"  cols = "); put(file,cols);
    put(file,"  conds = "); put(file,conds(conds'first)); new_line(file);
    Flag_Conditions(n,k,q,rows,cols,conds,flags,eqs);
    First_Solution(eqs.all,fail,x,res);
    put(file,"Residual of the solution : "); put(file,res,3);
    if fail
     then put_line(file," : failure.");
     else put_line(file," : success.");
    end if;
    declare
      sol : Solution(x'last);
      ls : Link_to_Solution;
    begin
      sol.t := Create(0.0); sol.m := 1; sol.v := x;
      sol.err := 0.0; sol.res := res; sol.rco := 1.0;
      ls := new Solution'(sol);
      Construct(ls,snd.sols);
    end;
  end Start_Solution;

  procedure Initialize_Solution_Nodes
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                nodes : in out Solnode_List;
                res : out double_float ) is

    tmp : Solnode_List := nodes;
    snd : Link_to_Solution_Node;
    res_node : double_float;
    fail : boolean;
    cnt : natural32 := 0;

  begin
    new_line(file);
    put_line(file,"Computing the solutions at the leaves ...");
    res := 0.0;
    while not Is_Null(tmp) loop
      snd := Head_Of(tmp);
      Start_Solution(file,n,k,conds,flags,snd,fail,res_node);
      Set_Head(tmp,snd);
      cnt := cnt + 1;
      if fail then
        put(file,"Failed to compute start solution at node ");
        put(file,cnt,1); new_line(file);
      end if;
      res := res + res_node;
      tmp := Tail_Of(tmp);
    end loop;
    put(file,"The sum of all residuals : "); put(file,res,3); new_line(file);
  end Initialize_Solution_Nodes;

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
              ( file : in file_type; n,k,level : in integer32;
                pl : in Poset_List; snd : in Link_to_Solution_Node;
                sps : in out Solution_Poset;
                conds : in Standard_Natural_VecVecs.VecVec;
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
      parent_snd_last : Solution_List;
      childnode : constant Checker_Posets.Link_to_Node
                := child.white(child.white'first);
      childconds : constant Standard_Natural_Vectors.Vector
                 := childnode.rows;     -- root
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
      parentconds : constant Standard_Natural_Vectors.Vector
                  := node.ps.white(node.ps.white'last).cols; -- leaf

      use Checker_Posets,Moving_Flag_Continuation;

    begin
      put(file,"Number of start solutions at child : ");
      put(file,Length_Of(snd.sols),1); new_line(file);
      if parent_snd = null then
        put_line(file,"No solution node at parent found!"); return;
      else
        put(file,"Number of solutions at parent node : ");
        put(file,Length_Of(parent_snd.sols),1); put_line(file,".");
        parent_snd_last := parent_snd.sols;
      end if;
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          put(file,"*** number of paths from child to the parent : ");
          put(file,childnode.coeff); put_line(file," ***");
          declare
            sols : Solution_List;
            totalflags : constant natural32 := natural32(flags'length);
            nbflags : constant integer32 := sps.m - level;
          begin
            put(file,"Calling Track_All_Paths_in_Posets with ");
            put(file,totalflags,1); put(file," fixed flags; level = ");
            put(file,level,1); 
            put(file,"  sps.m = "); put(file,sps.m,1); put_line(file,".");
            put(file,"number of flags = ");
            put(file,nbflags,1); put_line(file,".");
            Track_All_Paths_in_Poset
              (file,n,k,node.ps,
               conds(conds'last-nbflags+1..conds'last),
               flags(flags'last-nbflags+1..flags'last),snd,sols);
            Concat(parent_snd.sols,parent_snd_last,sols);
          end;
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
              ( file : in file_type; n,k : in integer32;
                ips : in out Intersection_Poset;
                sps : in out Solution_Poset;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in out Standard_Complex_VecMats.VecMat;
                sols : out Solution_List ) is

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;
    residual : double_float;

  begin
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    Initialize_Solution_Nodes(file,n,k,
      conds(conds'last..conds'last),
      flags(flags'last..flags'last),sps.nodes(sps.m),residual);
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
            (file,n,k,i-1,ips.nodes(i-1),snd,sps,conds,flags);
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
