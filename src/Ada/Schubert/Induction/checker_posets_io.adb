with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
--with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Checker_Boards_io;                  use Checker_Boards_io;
--with Checker_Moves;
with Checker_Localization_Patterns;      use Checker_Localization_Patterns;

package body Checker_Posets_io is

  procedure Write ( rows,cols : in Vector ) is
  begin
    Write(standard_output,rows,cols);
  end Write;

  procedure Write ( file : in file_type; rows,cols : in Vector ) is
  begin
    put(file,"("); Write_Bracket(file,rows);
    put(file,","); Write_Bracket(file,cols); put(file,")");
  end Write;

  procedure Write_Node ( nd : in Node ) is
  begin
    Write_Node(standard_output,nd);
  end Write_Node;

  procedure Write_Node ( file : in file_type; nd : in Node ) is
  begin
    put(file," + "); put(file,nd.coeff);
    Write(file,nd.rows,nd.cols);
  end Write_Node;

  procedure Write_Node_Labels ( level,pos : in integer32 ) is

  -- DESCRIPTION :
  --   Writes the labels of one node in the poset between round brackets.

  begin
    put("("); put(level,1); put(","); put(pos,1); put(")");
  end Write_Node_Labels;

  procedure Write_Parents ( p : in Vector; nd : in Node ) is

    cnt : integer32 := 0;

    procedure Write_Parent ( pnd : in Link_to_Node ) is
    begin
      cnt := cnt + 1;
      put("  parent node "); put(cnt,1); put(" : ");
      put(p); put(" : ");
      Write_Node(pnd.all); new_line;
    end Write_Parent;
    procedure Enum_Parents is new Enumerate_Parents(Write_Parent);

  begin
    Enum_Parents(nd);
  end Write_Parents;

  procedure Write_Nodes_in_Poset ( ps : in Poset; i : in integer32 ) is

    ptr : Link_to_Node;

  begin
    put(i,2); put(" : ");
    put(ps.black(i)); put(" : ");
    ptr := ps.white(i);             -- write coordinates
    while ptr /= null loop
      Write_Node(ptr.all);
      ptr := ptr.next_sibling;
    end loop;
    new_line;
  end Write_Nodes_in_Poset;

  procedure Write_Nodes_in_Poset
              ( file : in file_type; ps : in Poset; i : in integer32 ) is

    ptr : Link_to_Node;

  begin
    put(file,i,2); put(file," : ");
    put(file,ps.black(i)); put(file," : ");
    ptr := ps.white(i);             -- write coordinates
    while ptr /= null loop
      Write_Node(file,ptr.all);
      ptr := ptr.next_sibling;
    end loop;
    new_line(file);
  end Write_Nodes_in_Poset;

  procedure Write ( ps : in Poset ) is

    ptr : Link_to_Node;
    ind,pos : integer32;

  begin
    for i in ps.black'range loop
      Write_Nodes_in_Poset(ps,i);
      if i < ps.black'last then       -- write children
        ptr := ps.white(i);
        ind := 0;                     -- ind is index to current node
        while ptr /= null loop
          ind := ind + 1;
          Write_Node_Labels(i,ind); put("-> {");
          if ptr.stay_child /= null then
            pos := Position(ps.white(i+1).all,ptr.stay_child.all);
            Write_Node_Labels(i+1,pos);
            if ptr.swap_child /= null
             then put(",");
            end if;
          end if;
          if ptr.swap_child /= null then
            pos := Position(ps.white(i+1).all,ptr.swap_child.all);
            Write_Node_Labels(i+1,pos);
          end if;
          put("} ");
          ptr := ptr.next_sibling;
        end loop;
        new_line;
      end if;
    end loop;
  end Write;

  procedure Write_Patterns ( ps : in Poset ) is

    n : constant integer32 := ps.black(ps.black'first)'last;
    k : constant integer32 := ps.white(ps.white'first).rows'last;

  begin
    for i in ps.black'range loop
      put(i,2); put(" : ");
      put(ps.black(i)); put_line(" : ");
      declare
        tmp : Link_to_Node := ps.white(i);
      begin
        while tmp /= null loop
          Write_Node(tmp.all);
          put_line(" with pattern :");
          put(Row_Pattern(n,k,ps.black(i).all,tmp.rows,tmp.cols));
          tmp := tmp.next_sibling;
        end loop;
      end;
    end loop;
  end Write_Patterns;

  procedure Write_Formal_Sum ( nd : in Link_to_Node ) is

    tmp : Link_to_Node := nd;

  begin
    while tmp /= null loop
      put("+"); put(tmp.coeff);
      Write_Bracket(tmp.cols);
      tmp := tmp.next_sibling;
    end loop;
  end Write_Formal_Sum;

  procedure Write_Formal_Sum
              ( file : in file_type; nd : in Link_to_Node ) is

    tmp : Link_to_Node := nd;

  begin
    while tmp /= null loop
      put(file,"+");
      put(file,tmp.coeff);
      Write_Bracket(file,tmp.cols);
      tmp := tmp.next_sibling;
    end loop;
  end Write_Formal_Sum;

  procedure Write_Formal_Sums ( ps : in Poset ) is
  begin
    for i in ps.white'range loop
      put("Formal sum at level "); put(i,1); put(" : ");
      Write_Formal_Sum(ps.white(i));
      new_line;
    end loop;
  end Write_Formal_Sums;

  procedure Write_Final_Sum ( ps : in Poset ) is
  begin
    Write_Formal_Sum(ps.white(ps.white'last));
  end Write_Final_Sum;

  procedure Write_Final_Sum
             ( file : in file_type; ps : in Poset ) is
  begin
    Write_Formal_Sum(file,ps.white(ps.white'last));
  end Write_Final_Sum;

  procedure Write_Formal_Product ( ps : in Poset ) is
  begin
    put("+"); put(ps.white(ps.white'first).coeff);
    Write_Bracket(ps.white(ps.white'first).rows); put("*");
    Write_Bracket(ps.white(ps.white'first).cols);
  end Write_Formal_Product;

  procedure Write_Formal_Product
              ( file : in file_type; ps : in Poset ) is
  begin
    put(file,"+");
    put(file,ps.white(ps.white'first).coeff);
    Write_Bracket(file,ps.white(ps.white'first).rows);
    put(file,"*");
    Write_Bracket(file,ps.white(ps.white'first).cols);
  end Write_Formal_Product;

  procedure Write_Formal_Equation ( ps : in Poset ) is
  begin
    Write_Formal_Product(ps);
    put(" = ");
    Write_Final_Sum(ps);
    new_line;
  end Write_Formal_Equation;

  procedure Write_Formal_Equation
              ( file : in file_type; ps : in Poset ) is
  begin
    Write_Formal_Product(file,ps);
    put(file," = ");
    Write_Final_Sum(file,ps);
    new_line(file);
  end Write_Formal_Equation;

  procedure Write_Node_in_Path
               ( n,k : in integer32; ps : in Poset;
                 path : in Array_of_Nodes; index : in integer32 ) is
  begin
    Write_Node_in_Path(standard_output,n,k,ps,path,index);
  end Write_Node_in_Path;

  procedure Write_Node_in_Path
               ( file : in file_type; n,k : in integer32; ps : in Poset;
                 path : in Array_of_Nodes; index : in integer32 ) is

    previous : constant Link_to_Node := path(index-1);
    current : constant Link_to_Node := path(index);
    ptr : constant integer32 := ps.black'last - index + 1;
    p : constant Standard_Natural_Vectors.Vector(1..n) := ps.black(ptr+1).all;
    q : constant Standard_Natural_Vectors.Vector(1..n) := ps.black(ptr).all;
   -- fc : constant natural := Checker_Moves.Falling_Checker(p);
   -- tq : constant Standard_Natural_Matrices.Matrix(1..n,1..n) 
   --    := Checker_Localization_Patterns.Transformation(n,q(fc));
   -- mf : constant Standard_Natural_Matrices.Matrix(1..n,1..n)
   --    := Moving_Flag(p);

  begin
    put(file," -> node "); put(file,index,1); put(file," : ");
    put(file,"Move "); put(file,index-1,1);
    put(file," from"); put(file,p); put(file," to");
    put(file,q); put_line(file," :");
   -- Write_Permutations
   --   (file,p,q,previous.rows,previous.cols,current.rows,current.cols);
   -- Write_Patterns(file,Moving_Flag(p),Moving_Flag(q),
   --                Column_Pattern(n,k,p,previous.rows,previous.cols),
   --                Column_Pattern(n,k,q,current.rows,current.cols));
    Write_Permutations_and_Patterns
      (file,p,q,previous.rows,previous.cols,current.rows,current.cols,
       Moving_Flag(p),Moving_Flag(q),
       Column_Pattern(n,k,p,previous.rows,previous.cols),
       Column_Pattern(n,k,q,current.rows,current.cols));
   -- Write_Permutation(file,p,previous.rows,previous.cols,mf,tq);
   -- put_line(file,"Localization pattern at current node :");
   -- put(file,Row_Pattern(n,k,q,current.rows,current.cols));
  end Write_Node_in_Path;

  procedure Write_Path_in_Poset
              ( n,k : in integer32; ps : in Poset;
                path : in Array_of_Nodes; count : in integer32 ) is

  -- DESCRIPTION :
  --   Shows checker game information for one path through the poset.

  -- ON ENTRY :
  --   n        dimension of the ambient space;
  --   k        dimension of the k-plane;
  --   ps       checker poset needed for black checkers;
  --   path     sequence of nodes in a path of white checkers;
  --   count    counter for path number.

    leaf : constant Link_to_Node := path(path'first);
    root : constant Link_to_Node := path(path'last);
    p : constant Standard_Natural_Vectors.Vector(1..n)
      := ps.black(ps.black'last).all;
    q : constant Standard_Natural_Vectors.Vector(1..n)
      := ps.black(ps.black'first).all;

  begin
    put("Path "); put(count,1); put(" starts at ");
    Checker_Posets_io.Write_Node(leaf.all); new_line;
    put_line("Localization pattern at leaf :");
    put(Row_Pattern(n,k,p,leaf.rows,leaf.cols));
    for i in path'first+1..path'last loop
      Write_Node_in_Path(n,k,ps,path,i);
    end loop;
    put_line("Final checker configuration at root :");
    Write_Permutation(q,root.rows,root.cols,Moving_Flag(q));
  end Write_Path_in_Poset;

end Checker_Posets_io;
