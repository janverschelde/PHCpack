with unchecked_deallocation;
with Checker_Posets;

package body QuadDobl_Solution_Posets is

-- CONSTRUCTORS :

  function Create ( pnd : Link_to_Poset_Node ) return Solution_Node is

    res : Solution_Node;

  begin
    res.lpnd := pnd;
    return res;
  end Create;

  function Create ( ips : Intersection_Poset ) return Solution_Poset is
  
    res : Solution_Poset(ips.m);
    tmp : Poset_List;
    pnd : Link_to_Poset_Node;

  begin
    res.level := 0;
    for i in 1..ips.level loop
      if not Is_Null(ips.nodes(i)) then
        pnd := Head_Of(ips.nodes(i));
        declare
          snd : constant Solution_Node := create(pnd);
          lnd : constant Link_to_Solution_Node := new Solution_Node'(snd);
        begin
          Construct(lnd,res.nodes(i));
          res.last(i) := res.nodes(i);
        end;
        tmp := Tail_Of(ips.nodes(i));
        while not Is_Null(tmp) loop
          pnd := Head_Of(tmp);
          declare
            snd : constant Solution_Node := create(pnd);
            lnd : constant Link_to_Solution_Node := new Solution_Node'(snd);
          begin
            Append(res.nodes(i),res.last(i),lnd);
          end;
          tmp := Tail_Of(tmp);
        end loop;
      end if;
    end loop;
    return res;
  end Create;

-- SELECTOR :

  function Retrieve ( snl : Solnode_List;
                      rows,cols : Standard_Natural_Vectors.Vector )
                    return Link_to_Solution_Node is

    res : Link_to_Solution_Node := null;
    tmp : Solnode_List := snl;
    pnd : Link_to_Solution_Node;

  begin
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      declare
        p : constant Checker_Posets.Poset := pnd.lpnd.ps;
        r : constant Standard_Natural_Vectors.Vector
          := Checker_Posets.Root_Rows(p);
        c : constant Standard_Natural_Vectors.Vector
          := Checker_Posets.Root_Columns(p);
      begin
        if (Standard_Natural_Vectors.Equal(r,rows) 
            and then Standard_Natural_Vectors.Equal(c,cols))
         then return pnd;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Retrieve;

-- DESTRUCTORS :

  procedure free is
    new unchecked_deallocation(Solution_Node,Link_to_Solution_Node);

  procedure Clear ( snd : in out Solution_Node ) is
  begin
    clear(snd.sols);
  end Clear;

  procedure Clear ( snd : in out Link_to_Solution_Node ) is
  begin
    if snd /= null then
      Clear(snd.all);
      free(snd);
    end if;
  end Clear;

  procedure Clear ( nodes : in out Solnode_List ) is

    tmp : Solnode_List := nodes;
    pnd : Link_to_Solution_Node;

  begin
    while not Is_Null(tmp) loop
      pnd := Head_Of(tmp);
      Clear(pnd);
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Solution_Nodes.Clear(Lists_of_Solution_Nodes.List(nodes));
  end Clear;

  procedure Clear ( nodes : in out Array_of_Solnode_Lists ) is
  begin
    for i in nodes'range loop
      Clear(nodes(i));
    end loop;
  end Clear;

  procedure Clear ( sps : in out Solution_Poset ) is
  begin
    for i in sps.nodes'range loop
      Clear(sps.nodes(i));
    end loop;
  end Clear;

end QuadDobl_Solution_Posets;
