with unchecked_deallocation;
with Standard_Natural_Numbers_io;      use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;      use Standard_Integer_Numbers_io;
with Multprec_Complex_Numbers;         use Multprec_Complex_Numbers;
with Monomial_Hashing;                 use Monomial_Hashing;
with Multprec_Deflation_Matrices;      use Multprec_Deflation_Matrices;
with Multprec_Evaluate_Deflation_io;   use Multprec_Evaluate_Deflation_io;
with Standard_Natural_Vectors_io;      use Standard_Natural_Vectors_io;

package body Multprec_Evaluate_Deflation is

  function Create ( k : natural32 ) return Link_to_Eval_Tree is

    use Standard_Natural_Vectors;

    res : Link_to_Eval_Tree := null;
    first,last : Node_List;
    key_cnt : integer32 := 0;

    d : constant Vector(0..integer32(k)) := (0..integer32(k) => 0);
 
    procedure Generate ( i,m : in integer32; d : in Vector;
                         label : out integer32;
                         evt : in out Link_to_Eval_Tree ) is

    -- DESCRIPTION :
    --   Unwinds the application of d to A(i) at depth m in the tree.
    --   Only generates nodes which do not yet occur.

      to_generate : boolean;
      index : integer32;

    begin
      if Is_Null(first) then
        to_generate := true;
        declare
          res_gen : Eval_Tree(integer32(k),i);
        begin
          res_gen.key := 0;
          evt := new Eval_Tree'(res_gen);
          Append(first,last,evt);
        end;
        evt.d := d;
        evt.key := 0;
        label := 0;
        res := evt;
      else
        index := Key_In(first,d,i,key_cnt);
        if index = -1 then
          to_generate := true;
          declare
            res_gen : Eval_Tree(integer32(k),i);
          begin
            res_gen.key := 0;
            evt := new Eval_Tree'(res_gen);
            Append(first,last,evt);
          end;
          evt.d := d;
          for j in evt.e'range loop
            evt.e(j) := -1;
          end loop;
          key_cnt := key_cnt + 1;
          evt.key := natural32(key_cnt);
          label := integer32(evt.key);
        else
          to_generate := false;
          label := index;
          evt := null;
        end if;
      end if;
      if to_generate and then i > 0 then
        if d(i) = 0
         then Generate(i-1,m+1,d,evt.e(0),evt.c(0));
         else evt.c(0) := null;
              evt.e(0) := -1;
        end if;
        declare
          nd : Vector(d'range) := d;
        begin
          nd(0) := d(0)+1;
          nd(i) := 0;
          Generate(i-1,m+1,nd,evt.e(1),evt.c(1));
          nd(0) := d(0);
          for j in 1..i-1 loop
            if d(j) < 1 then
              nd(j) := d(j)+1;
              Generate(i-1,m+1,nd,evt.e(j+1),evt.c(j+1));
              nd(j) := d(j);
            else
              evt.e(j+1) := -1;
            end if;
          end loop;
        end;
      end if;
    end Generate;

  begin
    declare
      w : integer32;
    begin
      Generate(integer32(k),0,d,w,res);
    end;
    return res;
  end Create;

  function Is_In ( t : Eval_Tree; d : Standard_Natural_Vectors.Vector;
                   m : natural32 ) return boolean is

    res : boolean;

  begin
    if t.m < integer32(m) then
      res := false;  -- all children of t have lower m-value
    elsif t.m = integer32(m) then
      res := Standard_Natural_Vectors.Equal(t.d,d);
    else 
      res := false;
      for i in t.c'range loop
        if t.c(i) /= null then
          res := Is_In(t.c(i).all,d,m);
        end if;
        exit when res;
      end loop;
    end if;
    return res;
  end Is_In;

  function Key_In ( t : Eval_Tree; d : Standard_Natural_Vectors.Vector;
                    m,max_key : integer32 ) return integer32 is

    res : integer32 := -1;

  begin
    if (t.m < m) or (integer32(t.key) > max_key) then
      null;  -- children of t have lower m-value and larger key
    elsif t.m = m then
      if Standard_Natural_Vectors.Equal(t.d,d)
       then res := integer32(t.key);
      end if;
    else 
      for i in t.c'range loop
        if t.c(i) /= null then
          res := Key_In(t.c(i).all,d,m,max_key);
        end if;
        exit when (res > -1);
      end loop;
    end if;
    return res;
  end Key_In;

  function Key_In ( l : Node_List; d : Standard_Natural_Vectors.Vector;
                    m,max_key : integer32 ) return integer32 is

    res : integer32 := -1;
    tmp : Node_List := l;
    cnt : integer32 := 0;

  begin
    while not Is_Null(tmp) loop
      declare
        t : constant Link_to_Eval_Tree := Head_Of(tmp);
      begin
        if t.m = m then
          if Standard_Natural_Vectors.Equal(t.d,d)
           then res := integer32(t.key);
          end if;
        end if;
      end;
      exit when (res > -1);
      tmp := Tail_Of(tmp);
      cnt := cnt + 1;
      exit when cnt >= max_key;
    end loop;
    return res;
  end Key_In;

  function Look_Up ( t : Link_to_Eval_Tree; key : integer32 )
                   return Link_to_Eval_Tree is

    res : Link_to_Eval_Tree := null;

  begin
    if integer32(t.key) = key then
      return t;
    else
      for i in t.c'range loop
        if t.c(i) /= null
         then res := Look_Up(t.c(i),key);
        end if;
        exit when (res /= null);
      end loop;
      return res;
    end if;
  end Look_Up;

  function Is_Leaf ( t : Eval_Tree ) return boolean is
  begin
    for i in t.c'range loop
      if t.c(i) /= null or t.e(i) /= -1
       then return false;
      end if;
    end loop;
    return true;
  end Is_Leaf;

  function Node_Count ( evt : Eval_Tree ) return natural32 is

    res : natural32 := 0;

    procedure Count_Nodes ( t : in Eval_Tree; m : in natural32 ) is
    begin
      res := res + 1;
      if t.m > 0 then
        for i in t.c'range loop
          if t.c(i) /= null
           then Count_Nodes(t.c(i).all,m+1);
          end if;
        end loop;
      end if;
    end Count_Nodes;

  begin
    Count_Nodes(evt,0);
    return res;
  end Node_Count;

  function Different_Node_Count ( evt : Eval_Tree ) return natural32 is

    res : natural32 := 0;

    procedure Count_Nodes ( t : in Eval_Tree; m : in natural32 ) is
    begin
      if Key_In(evt,t.d,t.m,integer32(t.key)-1) = -1 then
        res := res + 1;
        if t.m > 0 then
          for i in t.c'range loop
            if t.c(i) /= null
             then Count_Nodes(t.c(i).all,m+1);
            end if;
          end loop;
        end if;
      end if;
    end Count_Nodes;

  begin
    Count_Nodes(evt,0);
    return res;
  end Different_Node_Count;

-- ENUMERATORS :

  procedure Enumerate_Multiplier_Derivatives ( k : in natural32 ) is   

    use Standard_Natural_Vectors;

    procedure Generate ( i,m : in integer32; d : in Vector ) is

    -- DESCRIPTION :
    --   Unwinds the application of d to A(i), at level m in the tree.

      cont : boolean;

    begin
      Multiplier_Derivative(d,natural32(i),natural32(m),false,cont);
      if i > 0 and cont then
        if d(i) = 0
         then Generate(i-1,m+1,d);
         else Multiplier_Derivative(d,natural32(i),natural32(m)+1,true,cont);
        end if;
        declare
          nd : Vector(d'range) := d;
        begin
          nd(0) := d(0)+1;
          nd(i) := 0;
          Generate(i-1,m+1,nd);
          nd(0) := d(0);
          for j in 1..i-1 loop
            if d(j) < 1 then
              nd(j) := d(j)+1;
              Generate(i-1,m+1,nd);
              nd(j) := d(j);
            end if;
          end loop;
        end;
      end if;
    end Generate;

  begin
    declare
      accu : constant Vector(0..integer32(k)) := (0..integer32(k) => 0);
    begin
      Generate(integer32(k),0,accu);
    end;
  end Enumerate_Multiplier_Derivatives;

  procedure Update_Stack
              ( sd : in out Standard_Natural_VecVecs.VecVec;
                sk : in out Standard_Natural_Vectors.Vector;
                sz : in out natural32;
                d : in Standard_Natural_Vectors.Vector;
                k : in natural32; found : out boolean ) is

    use Standard_Natural_Vectors;

  begin
    found := false;
    for i in sd'first..integer32(sz) loop
      if sk(i) = k
       then found := Equal(sd(i).all,d);
      end if;
      exit when found;
    end loop;
    if not found then
      sz := sz + 1;
      sd(integer32(sz)) := new Vector'(d);
      sk(integer32(sz)) := k; 
    end if;
  end Update_Stack;

-- AUXILIARIES to the EVALUATORS :

  procedure Display_Dimensions
               ( file : in file_type; i,nbc : in natural32;
                 d,nq,nv : in Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Displays the dimensions of the matrix.

    moncnt : natural32;

  begin
    put(file,"  ");
    if ((d(d'first) = 0) or (i > 0)) then
      put(file,nq(integer32(i)),1); put("-by-"); put(file,nbc,1);
      put_line(file," matrix");
    else
      moncnt := Monomial_Count(d(d'first),nv(0));
      put(file,moncnt,1); put(file," ");
      put(file,nq(integer32(i)),1); put(file,"-by-");
      put(file,nv(integer32(i)),1); put_line(file," matrices.");
    end if;
  end Display_Dimensions;

  function Lookup ( evt : Link_to_Eval_Tree;
                    key : integer32 ) return Link_to_Matrix is

  -- DESCRIPTION :
  --   Writes the derivative operator in evt corresponding to the key.

    ft : constant Link_to_Eval_Tree := Look_Up(evt,key);
    null_ftv : constant Link_to_Matrix := null;

  begin
    if ft /= null
     then return ft.v;
     else return null_ftv;
    end if;
  end Lookup;

  function Lookup ( file : file_type; evt : Link_to_Eval_Tree;
                    key : integer32 ) return Link_to_Matrix is

  -- DESCRIPTION :
  --   Writes the derivative operator in evt corresponding to the key.

    ft : constant Link_to_Eval_Tree := Look_Up(evt,key);
    null_ftv : constant Link_to_Matrix := null;

  begin
    if ft /= null then
      Write_Derivative_Operator(file,ft.d,natural32(ft.m)); new_line(file);
      return ft.v;
    else
      put_line(file,"  lookup failed!");
      return null_ftv;
    end if;
  end Lookup;

  procedure Lookup ( evt : in Link_to_Eval_Tree;
                     key : in integer32; value : out Link_to_Matrix;
                     d : out Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Searches the tree for the node with the given key and returns in
  --   value the value of the matrix and in d the derivative operator.

    ft : constant Link_to_Eval_Tree := Look_Up(evt,key);

  begin
    if ft /= null then
      value := ft.v;
      d := new Standard_Natural_Vectors.Vector'(ft.d);
    else
      value := null;
      d := null;
    end if;
  end Lookup;

  procedure Lookup ( file : in file_type; evt : in Link_to_Eval_Tree;
                     key : in integer32; value : out Link_to_Matrix;
                     d : out Standard_Natural_Vectors.Link_to_Vector ) is

  -- DESCRIPTION :
  --   Searches the tree for the node with the given key and returns in
  --   value the value of the matrix and in d the derivative operator.

    ft : constant Link_to_Eval_Tree := Look_Up(evt,key);

  begin
    if ft /= null then
      Write_Derivative_Operator(file,ft.d,natural32(ft.m)); new_line(file);
      value := ft.v;
      d := new Standard_Natural_Vectors.Vector'(ft.d);
    else
      put_line(file,"  lookup failed!");
      value := null;
      d := null;
    end if;
  end Lookup;

  procedure Evaluate_Nodes
              ( root,t : in Link_to_Eval_Tree;
                jrt : in Jacobian_Remember_Table;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R0 : in Standard_Natural_Vectors.Vector;
                B : in VecMat; Bl,h,x : in Multprec_Complex_VecVecs.VecVec;
                depth : in natural32 );

  procedure Evaluate_Nodes
              ( file : in file_type; root,t : in Link_to_Eval_Tree;
                jrt : in Jacobian_Remember_Table;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R0 : in Standard_Natural_Vectors.Vector;
                B : in VecMat; Bl,h,x : in Multprec_Complex_VecVecs.VecVec;
                depth : in natural32 );

  -- DESCRIPTION :
  --   Traverses the directed acyclic graph and evaluates recursively.
  --   Except for R0, Bl, "root" and "depth", the same specifications as Eval.

  -- ON ENTRY :
  --   file     for writing diagnostics to file;
  --   root     directed acyclic graph with derivative operators;
  --   t        current node in the directed acyclic graph;
  --   nd       remember table for Jacobian matrices up to order k,
  --            created with the procedure Create_Remember_Derivatives;
  --   monkeys  monomial key codes for every Jacobian matrix;
  --   k        order of the deflation matrix;
  --   nv       nv(i) is the number of columns in the i-th deflation matrix,
  --            nv is short for the number of variables;
  --   nq       nq(i) is the number of row in i-th deflation matrix,
  --            nq is short for the number of equations;
  --   R0       R0(i) is number of multipliers at stage i in deflation,
  --            R0(0) is number of original variables;
  --   B        B(i) is the random matrix used in the i-th stage;
  --   h        h(i) is the random vector to scale the i-th multipliers;
  --   depth    depth in the tree structure.

  -- ON RETURN :
  --   root     updated with values of deflation matrices;
  --   t        updated with values of deflation matrices.

  procedure Recursive_Evaluate_Nodes
              ( root,t : in Link_to_Eval_Tree;
                jrt : in Jacobian_Remember_Table;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R0 : in Standard_Natural_Vectors.Vector;
                B : in VecMat; Bl,h,x : in Multprec_Complex_VecVecs.VecVec;
                depth : in natural32 ) is

  -- DESCRIPTION :
  --   Handles the case when there are children in the tree.

    vals : VecMat(t.c'range);      -- values of children
    order : natural32 := 0;          -- order of derivatives on children
    chtp : Standard_Natural_VecVecs.VecVec(t.c'range);
     -- derivative operators for the children

    use Standard_Natural_Vectors;

  begin
    for i in t.c'range loop
      if t.c(i) = null then        -- an empty child
        if t.e(i) > -1 then        -- either occurs somewhere else
          Lookup(root,t.e(i),vals(i),chtp(i));
        end if;                    -- or is the zero matrix
      else                         -- evaluate the nonempty children
        Evaluate_Nodes
          (root,t.c(i),jrt,monkeys,k,nv,nq,R0,B,Bl,h,x,depth+1);
        vals(i) := t.c(i).v;
        chtp(i) := new Standard_Natural_Vectors.Vector'(t.c(i).d);
      end if;
    end loop;
    if (vals(0) = null) and (t.m > 1) then
      Assign_Children(t.v,t.m,R0,vals,B(t.m));
    else
      for i in chtp'range loop
        if chtp(i) /= null then
          order := Standard_Natural_Vectors.Sum(chtp(i).all);
        end if;
      end loop;
      Assign_Children
        (t.v,nq(0),t.d(0),natural32(t.m),order,chtp,R0,vals,monkeys,
         jrt,B(t.m),Bl(t.m));
    end if;
    Standard_Natural_VecVecs.Clear(chtp);
  end Recursive_Evaluate_Nodes;

  procedure Recursive_Evaluate_Nodes
              ( file : in file_type; root,t : in Link_to_Eval_Tree;
                jrt : in Jacobian_Remember_Table;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R0 : in Standard_Natural_Vectors.Vector;
                B : in VecMat; Bl,h,x : in Multprec_Complex_VecVecs.VecVec;
                depth : in natural32 ) is

  -- DESCRIPTION :
  --   Handles the case when there are children in the tree.

    vals : VecMat(t.c'range);      -- values of children
    order : natural32 := 0;          -- order of derivatives on children
    chtp : Standard_Natural_VecVecs.VecVec(t.c'range);
     -- derivative operators for the children

    use Standard_Natural_Vectors;

  begin
    for i in t.c'range loop
      if t.c(i) = null then        -- an empty child
        if t.e(i) > -1 then        -- either occurs somewhere else
          Write_Spaces(file,natural32(t.m));
          put(file,"  occurs as node ");
          put(file,t.e(i),1); put(file," : ");
          Lookup(file,root,t.e(i),vals(i),chtp(i));
        elsif i = 0 then           -- or is the zero matrix
          Write_Zero(file,natural32(t.m)+1);
        end if;
      else                         -- evaluate the nonempty children
        Evaluate_Nodes
          (file,root,t.c(i),jrt,monkeys,k,nv,nq,R0,B,Bl,h,x,depth+1);
        vals(i) := t.c(i).v;
        chtp(i) := new Standard_Natural_Vectors.Vector'(t.c(i).d);
      end if;
    end loop;
    if (vals(0) = null) and (t.m > 1) then
      put(file,"Need to multiply with B(");
      put(file,t.m,1); put_line(file,"), calling Assign_Children 1.");
      Assign_Children(file,t.v,t.m,R0,vals,B(t.m));
      put_line(file,"Done with the call to Assign_Children 1.");
    else
      put(file,"Calling Assign_Children 2 with t.m = ");
      put(file,t.m,1); new_line(file);
      put_line(file,"Derivative types of children : ");
      for i in chtp'range loop
        if chtp(i) /= null then
          order := Standard_Natural_Vectors.Sum(chtp(i).all);
          put(file,"child "); put(file,i,1); put(file," : ");
          put(file,chtp(i).all); new_line(file);
        end if;
      end loop;
      put(file," order = "); put(file,order,1); new_line(file);
      Assign_Children(file,t.v,nq(0),t.d(0),natural32(t.m),order,chtp,R0,
                      vals,monkeys,jrt,B(t.m),Bl(t.m));
      put_line(file,"Done with the call to Assign_Children 2.");
    end if;
    Standard_Natural_VecVecs.Clear(chtp);
  end Recursive_Evaluate_Nodes;

  procedure Evaluate_Nodes
              ( root,t : in Link_to_Eval_Tree;
                jrt : in Jacobian_Remember_Table;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R0 : in Standard_Natural_Vectors.Vector;
                B : in VecMat; Bl,h,x : in Multprec_Complex_VecVecs.VecVec;
                depth : in natural32 ) is

    key,ind,col : integer32;
    nbc,nbrows,nbcols : natural32;

  begin
    key := Key_In(root.all,t.d,t.m,integer32(t.key)-1);
    if key > -1 then               -- occurs already and we assign its value
      t.v := Lookup(root,key);
    else              -- we have a new matrix and continue
      nbc := Number_of_Columns(t.d,nv,R0(1..k),natural32(t.m));
      if t.m = 0 then -- we look in the Jacobian remember table
        ind := integer32(t.d(t.d'first));
        if jrt(ind)'length = 1 then
          t.v := new Matrix'(Zero_Matrix(nq(t.m),nv(t.m)));
          col := t.v'first(2);
          Dimensions(jrt(ind),integer32(nbrows),integer32(nbcols));
          Assign_from_Jacobian_Matrices(t.v,col,jrt(ind),integer32(nbcols));
        else
          t.v := new Matrix(1..1,1..1);
          t.v(1,1) := Create(ind);
        end if;
      else
        t.v := new Matrix'(Zero_Matrix(nq(t.m),nbc));
        if Standard_Natural_Vectors.Sum(t.d) = 0
         then Assign_Scaling_Coefficients(t.v,h(t.m));
        end if;
        if not Is_Leaf(t.all) then     -- we have children
          Recursive_Evaluate_Nodes
            (root,t,jrt,monkeys,k,nv,nq,R0,B,Bl,h,x,depth);
        end if;
      end if;
    end if;
  end Evaluate_Nodes;

  procedure Evaluate_Nodes
              ( file : in file_type; root,t : in Link_to_Eval_Tree;
                jrt : in Jacobian_Remember_Table;
                monkeys : in Standard_Natural64_VecVecs.VecVec;
                k : in integer32;
                nv,nq,R0 : in Standard_Natural_Vectors.Vector;
                B : in VecMat; Bl,h,x : in Multprec_Complex_VecVecs.VecVec;
                depth : in natural32 ) is

    key,ind,col : integer32;
    nbc,nbrows,nbcols : natural32;

  begin
    Write_Derivative_Operator(file,t.d,natural32(t.m),depth);
    Write_Spaces(file,natural32(t.m));
    put(file,"  Node "); put(t.key,1);
    key := Key_In(root.all,t.d,t.m,integer32(t.key)-1);
    if key > -1 then               -- occurs already and we assign its value
      put(file," occurs as node "); put(file,key,1); put(file," : ");
      t.v := Lookup(file,root,key);
    else              -- we have a new matrix and continue
      nbc := Number_of_Columns(t.d,nv,R0(1..k),natural32(t.m));
      Display_Dimensions(file,natural32(t.m),nbc,t.d,nq,nv);
      if t.m = 0 then -- we look in the Jacobian remember table
        ind := integer32(t.d(t.d'first));
       -- put(file,jrt(ind)'length,1); put_line(file," Jacobian matrices");
        if jrt(ind)'length = 1 then
          t.v := new Matrix'(Zero_Matrix(nq(t.m),nv(t.m)));
          col := t.v'first(2);
          Dimensions(jrt(ind),integer32(nbrows),integer32(nbcols));
          Assign_from_Jacobian_Matrices(t.v,col,jrt(ind),integer32(nbcols));
        else
          t.v := new Matrix(1..1,1..1);
          t.v(1,1) := Create(ind);
        end if;
      else
        t.v := new Matrix'(Zero_Matrix(nq(t.m),nbc));
        if Standard_Natural_Vectors.Sum(t.d) = 0
         then Assign_Scaling_Coefficients(t.v,h(t.m));
        end if;
        if not Is_Leaf(t.all) then     -- we have children
          Recursive_Evaluate_Nodes
            (file,root,t,jrt,monkeys,k,nv,nq,R0,B,Bl,h,x,depth);
        end if;
      end if;
    end if;
  end Evaluate_Nodes;

-- EVALUATORS :

  function Eval ( f : Eval_Poly_Sys; jm : Eval_Jaco_Mat;
                  t : Link_to_Eval_Tree; nd : Link_to_Eval_Node;
                  monkeys : Standard_Natural64_VecVecs.VecVec;
                  k : integer32;
                  nv,nq,R1 : Standard_Natural_Vectors.Vector;
                  B : VecMat; h,x : Multprec_Complex_VecVecs.VecVec )
                return Multprec_Complex_Vectors.Vector is

    k1 : constant integer32 := k-1;
    y : Multprec_Complex_Vectors.Vector(1..integer32(nq(k)));
    y0,y1 : Multprec_Complex_Vectors.Vector(1..integer32(nq(k1)));
    A1 : Matrix(1..integer32(nq(k1)),1..integer32(nv(k1)));
    Bl : Multprec_Complex_Vectors.Vector(1..integer32(nv(k1)));
    acc : Complex_Number;

  begin
    if k = 1 then
      y0 := Eval(f,x(0).all);
      A1 := Eval(jm,x(0).all);
    else
      declare
        t1 : Link_to_Eval_Tree := Create(natural32(k1));
      begin
        y0 := Eval(f,jm,t1,nd,monkeys,k1,nv(0..k1),nq(0..k1),R1(1..k1),
                   B(1..k1),h(1..k1),x(0..k1));
        A1 := Eval(t1,nd,monkeys,k1,nv(0..k1),nq(0..k1),R1(1..k1),
                   B(1..k1),h(1..k1),x(0..k1));
        Clear(t1);
      end;
    end if;
    y(y0'range) := y0;
    Bl := B(k).all*x(k).all;
    y1 := A1*Bl;
    for i in y1'range loop
      y(y0'last+i) := y1(i);
    end loop;
    y(y'last) := Create(-integer(1));
    for i in h(k)'range loop
      acc := h(k)(i)*x(k)(i);
      Add(y(y'last),acc);
    end loop;
    Clear(A1); Multprec_Complex_Vectors.Clear(Bl);
    return y;
  end Eval;

  function Eval ( t : Link_to_Eval_Tree; nd : Link_to_Eval_Node;
                  monkeys : Standard_Natural64_VecVecs.VecVec;
                  k : integer32;
                  nv,nq,R1 : Standard_Natural_Vectors.Vector;
                  B : VecMat; h,x : Multprec_Complex_VecVecs.VecVec )
                return Matrix is

    wrk : constant Link_to_Eval_Tree := t;
    jrt : Jacobian_Remember_Table(0..k);
    Bl : Multprec_Complex_VecVecs.VecVec(1..k);
    R0 : Standard_Natural_Vectors.Vector(0..k);

  begin
    Bl := Multiply(B,x(1..k));
    jrt := Evaluate_Jacobian_Remember_Table(nd,k,integer32(nv(0)),x(0).all);
    R0(0) := nv(0);
    R0(1..k) := R1;
    Evaluate_Nodes(t,wrk,jrt,monkeys,k,nv,nq,R0,B,Bl,h,x,0);
    Multprec_Complex_VecVecs.Clear(Bl);
    Multprec_Jacobian_Trees.Clear(jrt);
    if wrk.v = null
     then return Zero_Matrix(nq(k),nv(k));
     else return wrk.v.all;
    end if;
  end Eval;

  function Eval ( file : file_type;
                  t : Link_to_Eval_Tree; nd : Link_to_Eval_Node;
                  monkeys : Standard_Natural64_VecVecs.VecVec;
                  k : integer32;
                  nv,nq,R1 : Standard_Natural_Vectors.Vector;
                  B : VecMat; h,x : Multprec_Complex_VecVecs.VecVec )
                return Matrix is

    wrk : constant Link_to_Eval_Tree := t;
    jrt : Jacobian_Remember_Table(0..k);
    Bl : Multprec_Complex_VecVecs.VecVec(1..k);
    R0 : Standard_Natural_Vectors.Vector(0..k);

  begin
    put_line(file,"multiplying the B-matrices with the multipliers...");
    Bl := Multiply(B,x(1..k));
    put_line(file,"evaluating the remember table of Jacobian matrices...");
    jrt := Evaluate_Jacobian_Remember_Table(nd,k,integer32(nv(0)),x(0).all);
    put(file,"#Jacobian matrices : ");
    for i in 0..k loop
      put(file," "); put(file,natural32(jrt(i)'length),1);
    end loop;
    new_line(file);
    R0(0) := nv(0);
    R0(1..k) := R1;
    Evaluate_Nodes(file,t,wrk,jrt,monkeys,k,nv,nq,R0,B,Bl,h,x,0);
    Multprec_Complex_VecVecs.Clear(Bl);
    Multprec_Jacobian_Trees.Clear(jrt);
    if wrk.v = null
     then return Zero_Matrix(nq(k),nv(k));
     else return wrk.v.all;
    end if;
  end Eval;

-- DESTRUCTORS :

  procedure Clear ( nl : in out Node_List ) is

    tmp : Node_List := nl;
    evt : Link_to_Eval_Tree;

  begin
    while not Is_Null(tmp) loop
      evt := Head_Of(tmp);
      Clear(evt);
      tmp := Tail_Of(tmp);
    end loop;
    List_of_Nodes.Clear(List_of_Nodes.List(nl));
  end Clear;

  procedure Clear ( evt : in out Eval_Tree ) is
  begin
    for i in evt.c'range loop
      Clear(evt.c(i));
    end loop;
    Multprec_Complex_Matrices.Clear(evt.v);
  end Clear;

  procedure Clear ( evt : in out Link_to_Eval_Tree ) is

    procedure free is new unchecked_deallocation(Eval_Tree,Link_to_Eval_Tree);

  begin
    if evt /= null
     then Clear(evt.all);
          free(evt);
    end if;
  end Clear;

end Multprec_Evaluate_Deflation;
