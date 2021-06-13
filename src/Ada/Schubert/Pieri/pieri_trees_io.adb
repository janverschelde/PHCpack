with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Brackets_io;                        use Brackets_io;

package body Pieri_Trees_io is

  type Boolean_Array is array ( integer32 range <> ) of boolean;

  procedure put ( nd : in Pieri_Node ) is
  begin
    put(Standard_Output,nd);
  end put;

  procedure put ( file : in file_type; nd : in Pieri_Node ) is
  begin
    put(file,nd.node);
    put(file,"(c="); put(file,nd.c,1);
    put(file,",i="); put(file,nd.i,1);
    put(file,",h="); put(file,nd.h,1); put(file,")");
  end put;

  procedure put ( lnd : in Link_to_Pieri_Node ) is
  begin
    put(Standard_Output,lnd);
  end put;

  procedure put ( file : in file_type; lnd : in Link_to_Pieri_Node ) is
  begin
    put(file,lnd.all);
    if lnd.ancestor /= null then
      put(file," > ");
      put(file,lnd.ancestor);
    end if;
  end put;

  procedure put ( t : in Pieri_Tree ) is
  begin
    put(Standard_Output,t);
  end put;

  procedure put ( t : in Pieri_Tree; level : in natural32 ) is
  begin
    put(Standard_Output,t,level);
  end put;

  procedure put ( file : in file_type;
                  t : in Pieri_Tree; level : in natural32 ) is

    procedure Write_Node ( lnd : in Link_to_Pieri_Node;
                           continue : out boolean ) is
    begin
      put(file,lnd.all);
      continue := true;
    end Write_Node;
    procedure Write_Nodes is new Enumerate_Nodes(Write_Node);

  begin
    Write_Nodes(t,level);
  end put;

  procedure put ( file : in file_type; t : in Pieri_Tree ) is

    h : constant natural32 := Height(t);

  begin
    put(file,"Branching at "); put(file,t.branches); new_line(file);
    put(file,"Root node : "); put(file,t.root.node); new_line(file);
    for i in reverse 1..h loop
      put(file,"Nodes at level "); put(file,i,1); put(file," :");
      put(file,t,i); new_line(file);
    end loop;
  end put;

  function Last_Child ( nd : Pieri_Node; i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the ith child is the last child of the node.

  begin
    for j in (i+1)..nd.children'last loop
      if nd.children(j) /= null
       then return false;
      end if;
    end loop;
    return true;
  end Last_Child;

  procedure Write_Labels ( file : in file_type;
                           nd : in Pieri_Node; jump : in integer32;
                           last : in Boolean_Array ) is

  -- DESCRIPTION :
  --   Writes the contents of the Pieri node with the jump, taking into
  --   account which children appeared last.

  begin
    if nd.h /= 0
     then put(file,"   ");
    end if;
    for i in 1..(integer32(nd.h)-1) loop
      if last(i)
       then put(file,"     ");
       else put(file,"|    ");
      end if;
    end loop;
    if nd.h /= 0
     then put(file,"!-+"); put(file,jump,1);
    end if;
    put(file,nd);
    new_line(file);
  end Write_Labels;

  procedure Write_Nodes ( file : in file_type;
                          nd : in Pieri_Node; jump : in integer32;
                          last : in out Boolean_Array ) is

  -- DESCRIPTION :
  --   Writes the contents of the Pieri node preceded with the index of
  --   increase (jump), followed by the information on the children.

  begin
    Write_Labels(file,nd,jump,last);
    for i in nd.children'range loop
      if nd.children(i) /= null then
        last(integer32(nd.h)+1) := Last_Child(nd,i);
        Write_Nodes(file,nd.children(i).all,i,last);
      end if;
    end loop;
  end Write_Nodes;

  procedure Write_Tree ( t : in Pieri_Tree ) is
  begin
    Write_Tree(Standard_Output,t);
  end Write_Tree;

  procedure Write_Tree ( file : in file_type; t : in Pieri_Tree ) is

    last : Boolean_Array(1..integer32(Height(t)));

  begin
    Write_Nodes(file,t.root.all,0,last);
  end Write_Tree;

end Pieri_Trees_io;
