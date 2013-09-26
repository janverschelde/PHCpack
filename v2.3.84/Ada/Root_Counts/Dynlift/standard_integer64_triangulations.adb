with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Lists_of_Integer64_Vectors;         use Lists_of_Integer64_Vectors;

--with text_io,Simplices_io;              use text_io,Simplices_io;
--with Integer_Vectors_io;                use Integer_Vectors_io;
--with integer_io;                        use integer_io;

package body Standard_Integer64_Triangulations is

-- AUXILIAIRY HASH TABLE FOR Update_One :

  type Matrix_of_Triangulations is
    array ( integer32 range <>, integer32 range <> ) of Triangulation;

  type Hash_Table ( n : integer32 ) is record
    weight1,weight2 : Link_to_Vector;
    data : Matrix_of_Triangulations(0..n,0..n);
  end record;

  procedure Init ( ht : in out Hash_Table; first,last : in integer32 ) is

  -- Initializes the hash table with Null_Triangulation.
  -- Determines the weights for the hash function.
  -- The numbers first and last determine the length of the weight function.

  begin
    ht.weight1 := new Vector(first..last);
    ht.weight2 := new Vector(first..last);
    for i in first..last loop
      ht.weight1(i) := integer64(Random(-last,last));
      ht.weight2(i) := integer64(Random(-last,last));
    end loop;
    ht.weight1(last) := 1;  -- TO AVOID PROBLEMS WITH HIGH LIFTING
    ht.weight2(last) := 1;
    for i in ht.data'range(1) loop
      for j in ht.data'range(2) loop
        ht.data(i,j) := Null_Triangulation;
      end loop;
    end loop;
  end Init;

  procedure Hash_Function ( ht : in Hash_Table; s : in Simplex;
                            i,j : out integer32 ) is

  -- DESCRIPTION :
  --   Computes the entries in the hash table for the simplex s.

    r1 : constant integer64 := ht.weight1.all*Vertex(s,1);
    r2 : constant integer64 := ht.weight2.all*Vertex(s,2);

  begin
    i := integer32(r1) mod (ht.n + 1);
    j := integer32(r2) mod (ht.n + 1);
  end Hash_Function;

  procedure Update ( ht : in out Hash_Table; s : in Simplex ) is

   -- DESCRIPTION :
   --   Updates the hash table with the given simplex s.

    i,j : integer32;

  begin
    Hash_Function(ht,s,i,j);
    Construct(s,ht.data(i,j));
  end Update;

  function Is_In ( ht : Hash_Table; s : Simplex ) return boolean is

  -- DESCRIPTION :
  --   Checks whether the simplex belongs to the hash table.

    i,j : integer32;

  begin
    Hash_Function(ht,s,i,j);
    return Is_In(ht.data(i,j),s);
  end Is_In;

--  procedure Write_Distribution ( ht : in Hash_Table ) is
--
--   -- DESCRIPTION :
--   --   Writes the lengths of the lists in ht.
--
--  begin
--    for i in ht.data'range(1) loop
--      for j in ht.data'range(2) loop
--        put(Length_Of(ht.data(i,j)),1);  put(' ');
--      end loop;
--      new_line;
--    end loop;
--  end Write_Distribution;

  procedure Clear ( ht : in out Hash_Table ) is

  -- DESCRIPTION :
  --   Clears the allocated memory space for the hash table,
  --   only a shallow copy is performed.

  begin
    Clear(ht.weight1);
    Clear(ht.weight2);
    for i in ht.data'range(1) loop
      for j in ht.data'range(2) loop
        Lists_of_Simplices.Clear(Lists_of_Simplices.List(ht.data(i,j)));
      end loop;
    end loop;
  end Clear;

-- CREATORS :

  function Create ( s : Simplex ) return Triangulation is

    res : Triangulation;

  begin
    Construct(s,res);
    return res;
  end Create;

  function Is_Inner ( n : integer32; s : Simplex ) return boolean is
  begin
    for k in 1..n loop
      if Neighbor(s,k) = Null_Simplex
       then return false;
      end if;
    end loop;
    return true;
  end Is_Inner;

  procedure Update ( t : in out Triangulation; s : in out Simplex;
                     x : in Link_to_Vector ) is

    pos : constant vector(0..x'last) := Position(s,x.all);

  begin
    for k in x'range loop
      if Neighbor(s,k) = Null_Simplex then
        if pos(k-1)*pos(pos'last) > 0 then
          Update(s,x,k);
          Construct(Neighbor(s,k),t);
        end if;
      end if;
    end loop;
  end Update;

--  procedure Connect_and_Update 
--                   ( s : in out Simplex; t : in out Triangulation ) is
--
--   -- DESCRIPTION :
--   --   Connects the simplex with each other simplex in the triangulation
--   --   and adds it to the list t.
--
--    tmp : Triangulation := t;
--
--  begin
--    while not Is_Null(tmp) loop
--      declare
--        s2 : Simplex := Head_Of(tmp);
--      begin
--        Connect(s,s2);
--        --Set_Head(tmp,s2);
--      end;
--      tmp := Tail_Of(tmp);
--    end loop;
--    Construct(s,t);
--  end Connect_and_Update;

  procedure Connect_and_Update 
                   ( t1 : in Triangulation; t2 : in out Triangulation ) is

   -- DESCRIPTION :
   --   Connects all simplices in t1 properly with each other
   --   and adds them to the list t2.

    tmp1 : Triangulation := t1;
    tmp2 : Triangulation;
    s1,s2 : Simplex;

  begin
    while not Is_Null(tmp1) loop
      s1 := Head_Of(tmp1);
      tmp2 := Tail_Of(tmp1);
      while not Is_Null(tmp2) loop
        s2 := Head_Of(tmp2);
        Connect(s1,s2);
        tmp2 := Tail_Of(tmp2);
      end loop;
      tmp1 := Tail_Of(tmp1);
      Construct(s1,t2);
    end loop;
  end Connect_and_Update;

  procedure Update ( t : in out Triangulation; x : in Link_to_Vector;
                     newt : out Triangulation ) is

    tmp : Triangulation := t;
    s : Simplex;
    nwt : Triangulation; -- for the new simplices that will contain x

  begin
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      if not Is_Inner(x'last,s)
       then Update(nwt,s,x);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Connect_and_Update(nwt,t);
    newt := nwt;
  end Update;

  procedure Update ( t : in out Triangulation; x : in Link_to_Vector ) is

    nwt : Triangulation; -- for the new simplices that will contain x

  begin
    Update(t,x,nwt);
   -- SHALLOW CLEAR :
    Lists_of_Simplices.Clear(Lists_of_Simplices.List(nwt));
  end Update;

  procedure Collect_Vertices 
                ( s : in Simplex; L : in out List ) is

   -- DESCRIPTION :
   --   Collects the vertices in s and puts them in the list L.
   --   Auxiliary routine for Update_One.

    pts : constant VecVec := Vertices(s);

  begin
    for k in pts'range loop
      if not Is_In(L,pts(k).all) then
        declare
          newpt : constant Link_to_Vector := new Vector'(pts(k).all);
        begin
          Construct(newpt,L);
        end;
      end if;
    end loop;
  end Collect_Vertices;

  function Has_Vertices_in_List ( s : Simplex; L : List ) return boolean is

   -- DESCRIPTION :
   --   Returns true if the simplex s has vertices in the list L.

   pts : constant VecVec := Vertices(s);

  begin
    for k in pts'range loop
      if Is_In(L,pts(k).all)
       then return true;
      end if;
    end loop;
    return false;
  end Has_Vertices_in_List;

  procedure Update_One 
              ( t : in out Triangulation; x : in Link_to_Vector ) is

    s : Simplex := Head_Of(t);
    news,nei : Simplex;
    nwt,leaves : Triangulation;
    root : Hash_Table(2*x'last-1);
    pos : Vector(0..x'last);
    border_vertices : List;

    procedure Process_Tree_of_Updates is

     -- DESCRIPTION :
     --   Constructs a tree of adjencencies which consists of
     --   connected simplices on the edge of the triangulation.

      tmp,newleaves : Triangulation;
      si,neisi : Simplex;

    begin
     -- INITIALIZATION :
      --Construct(s,root);
      Init(root,x'first,x'last);
      Update(root,s);
      for k in x'range loop
        declare
          nei : constant Simplex := Neighbor(s,k);
        begin
          if (nei /= Null_Simplex) and then not Is_Vertex(nei,x.all)
           then Construct(nei,leaves);
          end if;
        end;
      end loop;
     -- PERFORM THE UPDATES :
      loop  
       -- INVARIANT CONDITION : 
       --   the simplices considered have vertices on the border close to x
        newleaves := Null_Triangulation;
        tmp := leaves;
        while not Is_Null(tmp) loop
          si := Head_Of(tmp);
         -- UPDATE root TO PREVENT TURNING BACK :
          Update(root,si);
         -- COMPUTE NEW SIMPLICES AND LEAVES :
          if not Is_Inner(x'last,si) then
            pos := Position(si,x.all);
            for k in x'range loop
           -- LOOK IN THE DIRECTION TOWARDS x :
              if pos(k-1)*pos(pos'last) > 0 then
                neisi := Neighbor(si,k);
                if neisi = Null_Simplex then
                  Update(si,x,k);         -- NEW SIMPLEX
                  neisi := Neighbor(si,k);
                  Construct(neisi,nwt);
                  Collect_Vertices(neisi,border_vertices);
                end if;
              end if;
            end loop;
          end if;
          for L in x'range loop   -- GO FURTHER 
            neisi := Neighbor(si,L);
            if (neisi /= Null_Simplex)
              and then not Is_Vertex(neisi,x.all)
              and then Has_Vertices_in_List(neisi,border_vertices)
              and then not Is_In(leaves,neisi)
              and then not Is_In(newleaves,neisi)
              and then not Is_In(root,neisi) 
             then Construct(neisi,newleaves);
            end if;
          end loop;
          tmp := Tail_Of(tmp);
        end loop;
        exit when (newleaves = Null_Triangulation);
        Lists_of_Simplices.Clear(Lists_of_Simplices.List(leaves));
        leaves := newleaves;
      end loop;
     -- SHALLOW CLEAR :
     -- Lists_of_Simplices.Clear(Lists_of_Simplices.List(root));
     -- put_line("Distribution of root : "); Write_Distribution(root);
      Clear(root); 
      Lists_of_Simplices.Clear(Lists_of_Simplices.List(leaves));
    end Process_Tree_of_Updates;

  begin
    nei := s;                 -- USED TO SEE WHETHER s WILL CHANGE
    pos := Position(s,x.all);
    Update_One(s,x,pos,news);
    Construct(news,nwt);
    Collect_Vertices(news,border_vertices);
    if s /= nei
     then pos := Position(s,x.all);  -- BECAUSE s IS CHANGED !!!
    end if;
    for k in x'range loop      -- COMPUTE OTHER NEIGHBORS 
      if pos(k-1)*pos(pos'last) > 0 then
        nei := Neighbor(s,k);
        if nei = Null_Simplex then
          Update(s,x,k);
          nei := Neighbor(s,k);
          Construct(nei,nwt);
          Collect_Vertices(nei,border_vertices);
        end if;
      end if;
    end loop;
    Process_Tree_of_Updates;
    Clear(border_vertices);
    Connect_and_Update(nwt,t);
    Lists_of_Simplices.Clear(Lists_of_Simplices.List(nwt));
  end Update_One;

  procedure Connect ( t : in out Triangulation ) is

    tmp1 : Triangulation := t;
    tmp2 : Triangulation;
    s1,s2 : Simplex;

  begin
    while not Is_Null(tmp1) loop
      s1 := Head_Of(tmp1);
      tmp2 := Tail_Of(tmp1);
      while not Is_Null(tmp2) loop
        s2 := Head_Of(tmp2);
        Connect(s1,s2);
        --Set_Head(tmp2,s2);
        tmp2 := Tail_Of(tmp2);
      end loop;
      --Set_Head(tmp1,s1);
      tmp1 := Tail_Of(tmp1);
    end loop;
  end Connect;

-- THE OPTIMAL DISCRETE CONSERVATIVE LIFTING FUNCTION :

--  procedure Write_Neighbors ( s : in Simplex ) is
--  begin
--    put("The normal of simplex s : "); put(Normal(s)); new_line;
--    put_line(" with vertices : "); put(s);
--    for k in Normal(s)'range loop
--      if Neighbor(s,k) /= Null_Simplex
--       then put("normal of a not null neighbors of s :"); 
--            put(Normal(Neighbor(s,k))); new_line;
--      end if;
--    end loop;
--  end Write_Neighbors;

--  procedure Write ( t : in Triangulation ) is
--
--    tmp : Triangulation := t;
--
--  begin
--    while not Is_Null(tmp) loop
--      Write_Neighbors(Head_Of(tmp));
--      tmp := Tail_Of(tmp);
--    end loop;
--  end Write;

  procedure Flatten ( t : in out Triangulation ) is

    tmp : Triangulation := t;
    s : Simplex;

  -- IMPLEMENTATION REQUIREMENT :
  --   Cells that are already flattened are grouped at the end of the list.

  begin
   -- put_line("THE TRIANGULATION BEFORE LIFTING : "); Write(t);
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      exit when Is_Flat(s);
     -- put_line("NEIGHBORS BEFORE FLATTENING : "); Write_Neighbors(s);
      Flatten(s);
     -- put_line("NEIGHBORS AFTER FLATTENING : "); Write_Neighbors(s);
      --Set_Head(tmp,s);
      tmp := Tail_Of(tmp);
    end loop;
   -- put_line("THE TRIANGULATION AFTER LIFTING : "); Write(t);
  end Flatten;

-- SELECTORS :

  function Is_Vertex ( t : Triangulation; x : Vector ) return boolean is
   
    tmp : Triangulation := t;

  begin
    while not Is_Null(tmp) loop
      if Is_Vertex(Head_Of(tmp),x)
       then return true;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return false;
  end Is_Vertex;

  function Vertices ( t : Triangulation ) return List is

   -- DESCRIPTION :
   --   Returns a list with all vertices in the simplices of t.

    res,res_last : List;
    tmp : Triangulation := t;

  begin
    res_last := res;
    while not Is_Null(tmp) loop
      declare
        s : constant Simplex := Head_Of(tmp);
        v : constant VecVec := Vertices(s);
      begin
        for i in v'range loop
          if not Is_In(res,v(i))
           then Append(res,res_last,v(i));
          end if;
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Vertices;

  function Vertices ( t : Triangulation; x : vector ) return List is

    res,res_last : List;
    tmp : Triangulation := t;

  begin
    res_last := res;
    while not Is_Null(tmp) loop
      declare
        s : constant Simplex := Head_Of(tmp);
      begin
        if Is_Vertex(s,x)
         then 
           declare
             v : constant VecVec := Vertices(s);
           begin
             for i in v'range loop
               if not Is_In(res,v(i))
                then Append(res,res_last,v(i));
               end if;
             end loop;
           end;
        end if;
        tmp := Tail_Of(tmp);
      end;
    end loop;
    return res;
  end Vertices;

  function Vertices ( t : Triangulation ) return VecVec is

    vertri : List := Vertices(t);
    res : VecVec(1..integer32(Length_Of(vertri)));
    tmp : List := vertri;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(vertri);
   -- Lists_of_Link_to_Integer_Vectors.Clear
   --                    (Lists_of_Link_to_Integer_Vectors.List(vertri));
    return res;
  end Vertices;

  function Vertices ( t : Triangulation; x : vector ) return VecVec is

    vertri : List := Vertices(t,x);
    res : VecVec(1..integer32(Length_Of(vertri)));
    tmp : List := vertri;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(vertri);
   -- Lists_of_Link_to_Integer_Vectors.Clear
   --                    (Lists_of_Link_to_Integer_Vectors.List(vertri));
    return res;
  end Vertices;

--  function Is_In1 ( t : Triangulation; x : Vector ) return boolean is
--
--    tmp : Triangulation := t;
--
--  begin
--    while not Is_Null(tmp) loop
--      if Is_In(Head_Of(tmp),x)
--       then return true;
--      end if;
--      tmp := Tail_Of(tmp);
--    end loop;
--    return false;
--  end Is_In1;

  function Is_In ( t : Triangulation; x : Vector ) return boolean is
  begin
    return Is_In_All(Head_Of(t),x);
  end Is_In;

  function Is_In ( t : Triangulation; x : Vector ) return Simplex is
  begin
    return Is_In_All(Head_Of(t),x);
  end Is_In;

  function Is_In ( t : Triangulation; s : Simplex ) return boolean is

    tmp : Triangulation := t;

  begin
   -- put("Length of triangulation : "); put(Length_Of(t),1); new_line;
    while not Is_Null(tmp) loop
      if Equal(s,Head_Of(tmp))
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Volume ( t : Triangulation ) return natural64 is

    tmp : Triangulation := t;
    vol : natural64 := 0;
 
  begin
    while not Is_Null(tmp) loop
      vol := vol + Volume(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return vol;
  end Volume;

-- DESTRUCTOR :

  procedure Clear ( t : in out Triangulation ) is

    tmp : Triangulation := t;

  begin
    while not Is_Null(tmp) loop
      declare
        s : Simplex := Head_Of(tmp);
      begin
        Clear(s);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Lists_of_Simplices.Clear(Lists_of_Simplices.List(t));
  end Clear;

end Standard_Integer64_Triangulations;
