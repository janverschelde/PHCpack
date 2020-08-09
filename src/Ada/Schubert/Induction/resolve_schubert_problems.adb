with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Multprec_Natural_Numbers_io;        use Multprec_Natural_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Natural_Matrices;
with Standard_Natural_Matrices_io;       use Standard_Natural_Matrices_io;
with Standard_Complex_Matrices_io;       use Standard_Complex_Matrices_io;
with DoblDobl_Complex_Matrices_io;       use DoblDobl_Complex_Matrices_io;
with QuadDobl_Complex_Matrices_io;       use QuadDobl_Complex_Matrices_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Checker_Moves;
with Checker_Posets,Checker_Posets_io;   use Checker_Posets_io;
with Checker_Localization_Patterns;
with Checker_Homotopies;
with Flag_Transformations;
with Start_Flag_Homotopies;              use Start_Flag_Homotopies;
with Moving_Flag_Homotopies;             use Moving_Flag_Homotopies;
with Checker_Poset_Deformations;

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
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_float ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    slnp : constant Checker_Posets.Poset := snd.lpnd.ps;
    rows : constant Standard_Natural_Vectors.Vector
         := Checker_Posets.Root_Rows(slnp);
    cols : constant Standard_Natural_Vectors.Vector(rows'range) := Flip(rows);
    q : constant Standard_Natural_Vectors.Vector
      := slnp.black(slnp.black'last).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    x : Standard_Complex_Vectors.Vector(1..integer32(dim));
    eqs : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

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

  procedure Start_Solution 
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    slnp : constant Checker_Posets.Poset := snd.lpnd.ps;
    rows : constant Standard_Natural_Vectors.Vector
         := Checker_Posets.Root_Rows(slnp);
    cols : constant Standard_Natural_Vectors.Vector(rows'range) := Flip(rows);
    q : constant Standard_Natural_Vectors.Vector
      := slnp.black(slnp.black'last).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    x : DoblDobl_Complex_Vectors.Vector(1..integer32(dim));
    eqs : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

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
      sol.t := Create(integer(0)); sol.m := 1; sol.v := x;
      sol.err := create(0.0);
      sol.res := res;
      sol.rco := create(1.0);
      ls := new Solution'(sol);
      Construct(ls,snd.sols);
    end;
  end Start_Solution;

  procedure Start_Solution 
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out quad_double ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    slnp : constant Checker_Posets.Poset := snd.lpnd.ps;
    rows : constant Standard_Natural_Vectors.Vector
         := Checker_Posets.Root_Rows(slnp);
    cols : constant Standard_Natural_Vectors.Vector(rows'range) := Flip(rows);
    q : constant Standard_Natural_Vectors.Vector
      := slnp.black(slnp.black'last).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    x : QuadDobl_Complex_Vectors.Vector(1..integer32(dim));
    eqs : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

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
      sol.t := Create(integer(0)); sol.m := 1; sol.v := x;
      sol.err := create(0.0);
      sol.res := res;
      sol.rco := create(1.0);
      ls := new Solution'(sol);
      Construct(ls,snd.sols);
    end;
  end Start_Solution;

  procedure Start_Solution 
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_float ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Solutions;

    slnp : constant Checker_Posets.Poset := snd.lpnd.ps;
    rows : constant Standard_Natural_Vectors.Vector
         := Checker_Posets.Root_Rows(slnp);
    cols : constant Standard_Natural_Vectors.Vector(rows'range) := Flip(rows);
    q : constant Standard_Natural_Vectors.Vector
      := slnp.black(slnp.black'last).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    x : Standard_Complex_Vectors.Vector(1..integer32(dim));
    eqs : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Flag_Conditions(n,k,q,rows,cols,conds,flags,eqs);
    First_Solution(eqs.all,fail,x,res);
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

  procedure Start_Solution 
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Solutions;

    slnp : constant Checker_Posets.Poset := snd.lpnd.ps;
    rows : constant Standard_Natural_Vectors.Vector
         := Checker_Posets.Root_Rows(slnp);
    cols : constant Standard_Natural_Vectors.Vector(rows'range) := Flip(rows);
    q : constant Standard_Natural_Vectors.Vector
      := slnp.black(slnp.black'last).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    x : DoblDobl_Complex_Vectors.Vector(1..integer32(dim));
    eqs : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Flag_Conditions(n,k,q,rows,cols,conds,flags,eqs);
    First_Solution(eqs.all,fail,x,res);
    declare
      sol : Solution(x'last);
      ls : Link_to_Solution;
    begin
      sol.t := Create(integer(0)); sol.m := 1; sol.v := x;
      sol.err := create(0.0);
      sol.res := res;
      sol.rco := create(1.0);
      ls := new Solution'(sol);
      Construct(ls,snd.sols);
    end;
  end Start_Solution;

  procedure Start_Solution 
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                fail : out boolean; res : out quad_double ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Solutions;

    slnp : constant Checker_Posets.Poset := snd.lpnd.ps;
    rows : constant Standard_Natural_Vectors.Vector
         := Checker_Posets.Root_Rows(slnp);
    cols : constant Standard_Natural_Vectors.Vector(rows'range) := Flip(rows);
    q : constant Standard_Natural_Vectors.Vector
      := slnp.black(slnp.black'last).all;
    locmap : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
           := Checker_Localization_Patterns.Column_Pattern(n,k,q,rows,cols);
    dim : constant natural32
        := Checker_Localization_Patterns.Degree_of_Freedom(locmap);
    x : QuadDobl_Complex_Vectors.Vector(1..integer32(dim));
    eqs : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    Flag_Conditions(n,k,q,rows,cols,conds,flags,eqs);
    First_Solution(eqs.all,fail,x,res);
    declare
      sol : Solution(x'last);
      ls : Link_to_Solution;
    begin
      sol.t := Create(integer(0)); sol.m := 1; sol.v := x;
      sol.err := create(0.0);
      sol.res := res;
      sol.rco := create(1.0);
      ls := new Solution'(sol);
      Construct(ls,snd.sols);
    end;
  end Start_Solution;

  procedure Initialize_Solution_Nodes
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                nodes : in out Standard_Solution_Posets.Solnode_List;
                res : out double_float ) is

    use Standard_Solution_Posets;

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

  procedure Initialize_Solution_Nodes
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                nodes : in out DoblDobl_Solution_Posets.Solnode_List;
                res : out double_double ) is

    use DoblDobl_Solution_Posets;

    tmp : Solnode_List := nodes;
    snd : Link_to_Solution_Node;
    res_node : double_double;
    fail : boolean;
    cnt : natural32 := 0;

  begin
    new_line(file);
    put_line(file,"Computing the solutions at the leaves ...");
    res := create(0.0);
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

  procedure Initialize_Solution_Nodes
              ( file : in file_type; n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                nodes : in out QuadDobl_Solution_Posets.Solnode_List;
                res : out quad_double ) is

    use QuadDobl_Solution_Posets;

    tmp : Solnode_List := nodes;
    snd : Link_to_Solution_Node;
    res_node : quad_double;
    fail : boolean;
    cnt : natural32 := 0;

  begin
    new_line(file);
    put_line(file,"Computing the solutions at the leaves ...");
    res := create(0.0);
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

  procedure Initialize_Solution_Nodes
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                nodes : in out Standard_Solution_Posets.Solnode_List;
                res : out double_float ) is

    use Standard_Solution_Posets;

    tmp : Solnode_List := nodes;
    snd : Link_to_Solution_Node;
    res_node : double_float;
    fail : boolean;
    cnt : natural32 := 0;

  begin
    res := 0.0;
    while not Is_Null(tmp) loop
      snd := Head_Of(tmp);
      Start_Solution(n,k,conds,flags,snd,fail,res_node);
      Set_Head(tmp,snd);
      cnt := cnt + 1;
      res := res + res_node;
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Solution_Nodes;

  procedure Initialize_Solution_Nodes
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                nodes : in out DoblDobl_Solution_Posets.Solnode_List;
                res : out double_double ) is

    use DoblDobl_Solution_Posets;

    tmp : Solnode_List := nodes;
    snd : Link_to_Solution_Node;
    res_node : double_double;
    fail : boolean;
    cnt : natural32 := 0;

  begin
    res := create(0.0);
    while not Is_Null(tmp) loop
      snd := Head_Of(tmp);
      Start_Solution(n,k,conds,flags,snd,fail,res_node);
      Set_Head(tmp,snd);
      cnt := cnt + 1;
      res := res + res_node;
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Solution_Nodes;

  procedure Initialize_Solution_Nodes
              ( n,k : in integer32;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                nodes : in out QuadDobl_Solution_Posets.Solnode_List;
                res : out quad_double ) is

    use QuadDobl_Solution_Posets;

    tmp : Solnode_List := nodes;
    snd : Link_to_Solution_Node;
    res_node : quad_double;
    fail : boolean;
    cnt : natural32 := 0;

  begin
    res := create(0.0);
    while not Is_Null(tmp) loop
      snd := Head_Of(tmp);
      Start_Solution(n,k,conds,flags,snd,fail,res_node);
      Set_Head(tmp,snd);
      cnt := cnt + 1;
      res := res + res_node;
      tmp := Tail_Of(tmp);
    end loop;
  end Initialize_Solution_Nodes;

  procedure Transform_Start_Solutions
              ( file : in file_type; n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

     use Standard_Complex_Solutions;

     p : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Identity_Permutation(natural32(n));
     q : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Reverse_Permutation(natural32(n));
     src_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,p,r_src,c_src);
     tgt_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,q,r_tgt,c_tgt);
     tmp : Solution_List := sols;
     ls : Link_to_Solution;
     x,y : Standard_Complex_Matrices.Matrix(1..n,1..k);
     dim : constant natural32
         := Checker_Localization_Patterns.Degree_of_Freedom(tgt_pat);
    -- invtm : constant Standard_Complex_Matrices.Matrix(1..n,1..n)
    --       := Standard_Matrix_Inversion.Inverse(tm);

     use Standard_Complex_Matrices;

  begin
    put_line(file,"The transformation matrix :");
    put(file,tm,3);
    put(file,"At source : ");
    put(file," p = "); put(file,p);
    put(file," r = "); put(file,r_src);
    put(file," c = "); put(file,c_src); new_line(file);
    put(file,"At target : ");
    put(file," p = "); put(file,q);
    put(file," r = "); put(file,r_tgt);
    put(file," c = "); put(file,c_tgt); new_line(file);
    put_line(file,"The pattern at the source : "); put(file,src_pat);
    put_line(file,"The pattern at the target : "); put(file,tgt_pat);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      x := Checker_Localization_Patterns.Map(src_pat,ls.v);
      put_line(file,"The solution plane at the source : "); put(file,x,3);
      y := tm*x; --  y := invtm*x;
      put_line(file,"After the transformation : "); put(file,y,3);
      Checker_Homotopies.Normalize_and_Reduce_to_Fit(tgt_pat,y);
      put_line(file,"The transformed plane at target :"); put(file,y,3);
      declare
        z : constant Standard_Complex_Vectors.Vector(1..integer32(dim))
          := Checker_Localization_Patterns.Map(tgt_pat,y);
        s : Solution(integer32(dim));
      begin
        s.t := ls.t;
        s.m := ls.m;
        s.v := z;
        s.err := ls.err;
        s.rco := ls.rco;
        s.res := ls.res;
        Clear(ls);
        ls := new Solution'(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Transform_Start_Solutions;

  procedure Transform_Start_Solutions
              ( file : in file_type; n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

     use DoblDobl_Complex_Solutions;

     p : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Identity_Permutation(natural32(n));
     q : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Reverse_Permutation(natural32(n));
     src_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,p,r_src,c_src);
     tgt_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,q,r_tgt,c_tgt);
     tmp : Solution_List := sols;
     ls : Link_to_Solution;
     x,y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
     dim : constant natural32
         := Checker_Localization_Patterns.Degree_of_Freedom(tgt_pat);

     use DoblDobl_Complex_Matrices;

  begin
    put_line(file,"The transformation matrix :");
    put(file,tm,3);
    put(file,"At source : ");
    put(file," p = "); put(file,p);
    put(file," r = "); put(file,r_src);
    put(file," c = "); put(file,c_src); new_line(file);
    put(file,"At target : ");
    put(file," p = "); put(file,q);
    put(file," r = "); put(file,r_tgt);
    put(file," c = "); put(file,c_tgt); new_line(file);
    put_line(file,"The pattern at the source : "); put(file,src_pat);
    put_line(file,"The pattern at the target : "); put(file,tgt_pat);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      x := Checker_Localization_Patterns.Map(src_pat,ls.v);
      put_line(file,"The solution plane at the source : "); put(file,x,3);
      y := tm*x;
      put_line(file,"After the transformation : "); put(file,y,3);
      Checker_Homotopies.Normalize_and_Reduce_to_Fit(tgt_pat,y);
      put_line(file,"The transformed plane at target :"); put(file,y,3);
      declare
        z : constant DoblDobl_Complex_Vectors.Vector(1..integer32(dim))
          := Checker_Localization_Patterns.Map(tgt_pat,y);
        s : Solution(integer32(dim));
      begin
        s.t := ls.t;
        s.m := ls.m;
        s.v := z;
        s.err := ls.err;
        s.rco := ls.rco;
        s.res := ls.res;
        Clear(ls);
        ls := new Solution'(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Transform_Start_Solutions;

  procedure Transform_Start_Solutions
              ( file : in file_type; n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

     use QuadDobl_Complex_Solutions;

     p : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Identity_Permutation(natural32(n));
     q : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Reverse_Permutation(natural32(n));
     src_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,p,r_src,c_src);
     tgt_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,q,r_tgt,c_tgt);
     tmp : Solution_List := sols;
     ls : Link_to_Solution;
     x,y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
     dim : constant natural32
         := Checker_Localization_Patterns.Degree_of_Freedom(tgt_pat);

     use QuadDobl_Complex_Matrices;

  begin
    put_line(file,"The transformation matrix :");
    put(file,tm,3);
    put(file,"At source : ");
    put(file," p = "); put(file,p);
    put(file," r = "); put(file,r_src);
    put(file," c = "); put(file,c_src); new_line(file);
    put(file,"At target : ");
    put(file," p = "); put(file,q);
    put(file," r = "); put(file,r_tgt);
    put(file," c = "); put(file,c_tgt); new_line(file);
    put_line(file,"The pattern at the source : "); put(file,src_pat);
    put_line(file,"The pattern at the target : "); put(file,tgt_pat);
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      x := Checker_Localization_Patterns.Map(src_pat,ls.v);
      put_line(file,"The solution plane at the source : "); put(file,x,3);
      y := tm*x;
      put_line(file,"After the transformation : "); put(file,y,3);
      Checker_Homotopies.Normalize_and_Reduce_to_Fit(tgt_pat,y);
      put_line(file,"The transformed plane at target :"); put(file,y,3);
      declare
        z : constant QuadDobl_Complex_Vectors.Vector(1..integer32(dim))
          := Checker_Localization_Patterns.Map(tgt_pat,y);
        s : Solution(integer32(dim));
      begin
        s.t := ls.t;
        s.m := ls.m;
        s.v := z;
        s.err := ls.err;
        s.rco := ls.rco;
        s.res := ls.res;
        Clear(ls);
        ls := new Solution'(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Transform_Start_Solutions;

  procedure Transform_Start_Solutions
              ( n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in Standard_Complex_Matrices.Matrix;
                sols : in out Standard_Complex_Solutions.Solution_List ) is

     use Standard_Complex_Solutions;

     p : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Identity_Permutation(natural32(n));
     q : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Reverse_Permutation(natural32(n));
     src_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,p,r_src,c_src);
     tgt_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,q,r_tgt,c_tgt);
     tmp : Solution_List := sols;
     ls : Link_to_Solution;
     x,y : Standard_Complex_Matrices.Matrix(1..n,1..k);
     dim : constant natural32
         := Checker_Localization_Patterns.Degree_of_Freedom(tgt_pat);

     use Standard_Complex_Matrices;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      x := Checker_Localization_Patterns.Map(src_pat,ls.v);
      y := tm*x;
      Checker_Homotopies.Normalize_and_Reduce_to_Fit(tgt_pat,y);
      declare
        z : constant Standard_Complex_Vectors.Vector(1..integer32(dim))
          := Checker_Localization_Patterns.Map(tgt_pat,y);
        s : Solution(integer32(dim));
      begin
        s.t := ls.t;
        s.m := ls.m;
        s.v := z;
        s.err := ls.err;
        s.rco := ls.rco;
        s.res := ls.res;
        Clear(ls);
        ls := new Solution'(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Transform_Start_Solutions;

  procedure Transform_Start_Solutions
              ( n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in DoblDobl_Complex_Matrices.Matrix;
                sols : in out DoblDobl_Complex_Solutions.Solution_List ) is

     use DoblDobl_Complex_Solutions;

     p : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Identity_Permutation(natural32(n));
     q : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Reverse_Permutation(natural32(n));
     src_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,p,r_src,c_src);
     tgt_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,q,r_tgt,c_tgt);
     tmp : Solution_List := sols;
     ls : Link_to_Solution;
     x,y : DoblDobl_Complex_Matrices.Matrix(1..n,1..k);
     dim : constant natural32
         := Checker_Localization_Patterns.Degree_of_Freedom(tgt_pat);

     use DoblDobl_Complex_Matrices;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      x := Checker_Localization_Patterns.Map(src_pat,ls.v);
      y := tm*x;
      Checker_Homotopies.Normalize_and_Reduce_to_Fit(tgt_pat,y);
      declare
        z : constant DoblDobl_Complex_Vectors.Vector(1..integer32(dim))
          := Checker_Localization_Patterns.Map(tgt_pat,y);
        s : Solution(integer32(dim));
      begin
        s.t := ls.t;
        s.m := ls.m;
        s.v := z;
        s.err := ls.err;
        s.rco := ls.rco;
        s.res := ls.res;
        Clear(ls);
        ls := new Solution'(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Transform_Start_Solutions;

  procedure Transform_Start_Solutions
              ( n,k : in integer32;
                r_src,c_src,r_tgt,c_tgt : in Standard_Natural_Vectors.Vector;
                tm : in QuadDobl_Complex_Matrices.Matrix;
                sols : in out QuadDobl_Complex_Solutions.Solution_List ) is

     use QuadDobl_Complex_Solutions;

     p : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Identity_Permutation(natural32(n));
     q : constant Standard_Natural_Vectors.Vector
       := Checker_Moves.Reverse_Permutation(natural32(n));
     src_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,p,r_src,c_src);
     tgt_pat : constant Standard_Natural_Matrices.Matrix(1..n,1..k)
       := Checker_Localization_Patterns.Column_Pattern(n,k,q,r_tgt,c_tgt);
     tmp : Solution_List := sols;
     ls : Link_to_Solution;
     x,y : QuadDobl_Complex_Matrices.Matrix(1..n,1..k);
     dim : constant natural32
         := Checker_Localization_Patterns.Degree_of_Freedom(tgt_pat);

     use QuadDobl_Complex_Matrices;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      x := Checker_Localization_Patterns.Map(src_pat,ls.v);
      y := tm*x;
      Checker_Homotopies.Normalize_and_Reduce_to_Fit(tgt_pat,y);
      declare
        z : constant QuadDobl_Complex_Vectors.Vector(1..integer32(dim))
          := Checker_Localization_Patterns.Map(tgt_pat,y);
        s : Solution(integer32(dim));
      begin
        s.t := ls.t;
        s.m := ls.m;
        s.v := z;
        s.err := ls.err;
        s.rco := ls.rco;
        s.res := ls.res;
        Clear(ls);
        ls := new Solution'(s);
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Transform_Start_Solutions;

  procedure Connect_Checker_Posets_to_Count
              ( pl : in Poset_List; nd : in Poset_Node;
                vrblvl : in integer32 := 0 ) is

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
     -- parentconds : constant Standard_Natural_Vectors.Vector
     --             := node.ps.white(node.ps.white'last).cols; -- leaf

      use Checker_Posets;

    begin
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
        end if;
        exit when (gamenode.next_sibling = null);
        gamenode := gamenode.next_sibling;
      end loop;
    end Connect_Parent;
    procedure Connect_Parents is
      new Intersection_Posets.Enumerate_Parents(Connect_Parent);

  begin
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Count 1 ...");
    end if;
    Connect_Parents(pl,nd);
  end Connect_Checker_Posets_to_Count;

  procedure Connect_Checker_Posets_to_Count
              ( file : in file_type;
                pl : in Poset_List; nd : in Poset_Node;
                vrblvl : in integer32 := 0 ) is

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
     -- parentconds : constant Standard_Natural_Vectors.Vector
     --             := node.ps.white(node.ps.white'last).cols; -- leaf

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
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Count 2 ...");
    end if;
    Connect_Parents(pl,nd);
  end Connect_Checker_Posets_to_Count;

  procedure Connect_Checker_Posets_to_Track
              ( file : in file_type; n,k,level,nt : in integer32;
                tol : in double_float;
                pl : in Poset_List;
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                tmfo : in Standard_Complex_Matrices.Link_to_Matrix;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    nd : constant Link_to_Poset_Node := snd.lpnd;

    use Standard_Complex_Solutions;
    use Standard_Solution_Posets;

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
     -- childrows : constant Standard_Natural_Vectors.Vector
     --           := Flip(childconds);
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
      parentconds : constant Standard_Natural_Vectors.Vector
                  := node.ps.white(node.ps.white'last).cols; -- leaf
      parentrows : constant Standard_Natural_Vectors.Vector
                 := node.ps.white(node.ps.white'last).rows;
      startsols : Solution_List;

      use Standard_Complex_Matrices;
      use Checker_Posets,Checker_Poset_Deformations;

    begin
      put(file,"Number of start solutions at child : ");
      put(file,Length_Of(snd.sols),1); new_line(file);
     -- if parent_snd = null then
     --   put_line(file,"No solution node at parent found!"); return;
     -- else
     --   put(file,"Number of solutions at parent node : ");
     --   put(file,Length_Of(parent_snd.sols),1); put_line(file,".");
     -- end if;
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          put(file,"*** number of paths from child to the parent : ");
          put(file,childnode.coeff); put_line(file," ***");
          Copy(snd.sols,startsols);
          if tmfo /= null then
            put_line(file,"Transforming the start solutions ...");
            put(file,"cols at leaf of parent : ");
            put(file,parentconds); new_line(file);
            put(file,"rows at leaf of parent : ");
            put(file,parentrows); new_line(file);
            put(file,"cols at game node at root : ");
            Transform_Start_Solutions
              (file,n,k,childnode.rows,childnode.cols, -- childconds,
               gamenode.rows,gamenode.cols,
              -- Flip(gamenode.cols),gamenode.cols,
               tmfo.all,startsols);
          end if;
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
              (file,n,k,nt,node.ps,childconds,verify,minrep,tosqr,
               conds(conds'last-nbflags+1..conds'last),
               flags(flags'last-nbflags+1..flags'last),
               tol,startsols,sols,rpt,vrblvl-1);
            Push(sols,parent_snd.sols);
            put(file,"Before push : #sols returned = ");
            put(file,Length_Of(sols),1); new_line(file);
            put(file,"After push to parent_snd : #sols = ");
            put(file,Length_of(parent_snd.sols),1); new_line(file);
          end;
          Deep_Clear(startsols);
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
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Track 1 ...");
    end if;
    Connect_Parents(pl,nd.all);
  end Connect_Checker_Posets_to_Track;

  procedure Connect_Checker_Posets_to_Track
              ( file : in file_type; n,k,level,nt : in integer32;
                tol : in double_float;
                pl : in Poset_List;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in DoblDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 ) is

    nd : constant Link_to_Poset_Node := snd.lpnd;

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Posets;

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
     -- childrows : constant Standard_Natural_Vectors.Vector
     --           := Flip(childconds);
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
      parentconds : constant Standard_Natural_Vectors.Vector
                  := node.ps.white(node.ps.white'last).cols; -- leaf
      parentrows : constant Standard_Natural_Vectors.Vector
                 := node.ps.white(node.ps.white'last).rows;
      startsols : Solution_List;

      use DoblDobl_Complex_Matrices;
      use Checker_Posets,Checker_Poset_Deformations;

    begin
      put(file,"Number of start solutions at child : ");
      put(file,Length_Of(snd.sols),1); new_line(file);
     -- if parent_snd = null then
     --   put_line(file,"No solution node at parent found!"); return;
     -- else
     --   put(file,"Number of solutions at parent node : ");
     --   put(file,Length_Of(parent_snd.sols),1); put_line(file,".");
     -- end if;
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          put(file,"*** number of paths from child to the parent : ");
          put(file,childnode.coeff); put_line(file," ***");
          Copy(snd.sols,startsols);
          if tmfo /= null then
            put_line(file,"Transforming the start solutions ...");
            put(file,"cols at leaf of parent : ");
            put(file,parentconds); new_line(file);
            put(file,"rows at leaf of parent : ");
            put(file,parentrows); new_line(file);
            put(file,"cols at game node at root : ");
            Transform_Start_Solutions
              (file,n,k,childnode.rows,childnode.cols, -- childconds,
               gamenode.rows,gamenode.cols,
              -- Flip(gamenode.cols),gamenode.cols,
               tmfo.all,startsols);
          end if;
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
              (file,n,k,nt,node.ps,childconds,verify,minrep,tosqr,
               conds(conds'last-nbflags+1..conds'last),
               flags(flags'last-nbflags+1..flags'last),tol,startsols,sols,
               vrblvl-1);
            Push(sols,parent_snd.sols);
            put(file,"Before push : #sols returned = ");
            put(file,Length_Of(sols),1); new_line(file);
            put(file,"After push to parent_snd : #sols = ");
            put(file,Length_of(parent_snd.sols),1); new_line(file);
          end;
          Deep_Clear(startsols);
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
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Track 2 ...");
    end if;
    Connect_Parents(pl,nd.all);
  end Connect_Checker_Posets_to_Track;

  procedure Connect_Checker_Posets_to_Track
              ( file : in file_type; n,k,level,nt : in integer32;
                tol : in double_float;
                pl : in Poset_List;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in QuadDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 ) is

    nd : constant Link_to_Poset_Node := snd.lpnd;

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Posets;

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
     -- childrows : constant Standard_Natural_Vectors.Vector
     --           := Flip(childconds);
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
      parentconds : constant Standard_Natural_Vectors.Vector
                  := node.ps.white(node.ps.white'last).cols; -- leaf
      parentrows : constant Standard_Natural_Vectors.Vector
                 := node.ps.white(node.ps.white'last).rows;
      startsols : Solution_List;

      use QuadDobl_Complex_Matrices;
      use Checker_Posets,Checker_Poset_Deformations;

    begin
      put(file,"Number of start solutions at child : ");
      put(file,Length_Of(snd.sols),1); new_line(file);
     -- if parent_snd = null then
     --   put_line(file,"No solution node at parent found!"); return;
     -- else
     --   put(file,"Number of solutions at parent node : ");
     --   put(file,Length_Of(parent_snd.sols),1); put_line(file,".");
     -- end if;
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          put(file,"*** number of paths from child to the parent : ");
          put(file,childnode.coeff); put_line(file," ***");
          Copy(snd.sols,startsols);
          if tmfo /= null then
            put_line(file,"Transforming the start solutions ...");
            put(file,"cols at leaf of parent : ");
            put(file,parentconds); new_line(file);
            put(file,"rows at leaf of parent : ");
            put(file,parentrows); new_line(file);
            put(file,"cols at game node at root : ");
            Transform_Start_Solutions
              (file,n,k,childnode.rows,childnode.cols, -- childconds,
               gamenode.rows,gamenode.cols,
              -- Flip(gamenode.cols),gamenode.cols,
               tmfo.all,startsols);
          end if;
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
              (file,n,k,nt,node.ps,childconds,verify,minrep,tosqr,
               conds(conds'last-nbflags+1..conds'last),
               flags(flags'last-nbflags+1..flags'last),tol,startsols,sols,
               vrblvl-1);
            Push(sols,parent_snd.sols);
            put(file,"Before push : #sols returned = ");
            put(file,Length_Of(sols),1); new_line(file);
            put(file,"After push to parent_snd : #sols = ");
            put(file,Length_of(parent_snd.sols),1); new_line(file);
          end;
          Deep_Clear(startsols);
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
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Track 3 ...");
    end if;
    Connect_Parents(pl,nd.all);
  end Connect_Checker_Posets_to_Track;

  procedure Connect_Checker_Posets_to_Track
              ( n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in Standard_Solution_Posets.Link_to_Solution_Node;
                tmfo : in Standard_Complex_Matrices.Link_to_Matrix;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    nd : constant Link_to_Poset_Node := snd.lpnd;

    use Standard_Complex_Solutions;
    use Standard_Solution_Posets;

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
     -- childrows : constant Standard_Natural_Vectors.Vector
     --           := Flip(childconds);
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
     -- parentconds : constant Standard_Natural_Vectors.Vector
     --             := node.ps.white(node.ps.white'last).cols; -- leaf
     -- parentrows : constant Standard_Natural_Vectors.Vector
     --            := node.ps.white(node.ps.white'last).rows;
      startsols : Solution_List;

      use Standard_Complex_Matrices;
      use Checker_Posets,Checker_Poset_Deformations;

    begin
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          Copy(snd.sols,startsols);
          if tmfo /= null then
            Transform_Start_Solutions
              (n,k,childnode.rows,childnode.cols,
               gamenode.rows,gamenode.cols,tmfo.all,startsols);
          end if;
          declare
            sols : Solution_List;
           -- totalflags : constant natural32 := natural32(flags'length);
            nbflags : constant integer32 := sps.m - level;
          begin
            Track_All_Paths_in_Poset
              (n,k,nt,node.ps,childconds,minrep,tosqr,
               conds(conds'last-nbflags+1..conds'last),
               flags(flags'last-nbflags+1..flags'last),tol,startsols,sols,
               rpt,vrblvl-1);
            Push(sols,parent_snd.sols);
          end;
          Deep_Clear(startsols);
        end if;
        exit when (gamenode.next_sibling = null);
        gamenode := gamenode.next_sibling;
      end loop;
    end Connect_Parent;
    procedure Connect_Parents is
      new Intersection_Posets.Enumerate_Parents(Connect_Parent);

  begin
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Track 4 ...");
    end if;
    Connect_Parents(pl,nd.all);
  end Connect_Checker_Posets_to_Track;

  procedure Connect_Checker_Posets_to_Track
              ( n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in DoblDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in DoblDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 ) is

    nd : constant Link_to_Poset_Node := snd.lpnd;

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Solution_Posets;

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
     -- childrows : constant Standard_Natural_Vectors.Vector
     --           := Flip(childconds);
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
     -- parentconds : constant Standard_Natural_Vectors.Vector
     --             := node.ps.white(node.ps.white'last).cols; -- leaf
     -- parentrows : constant Standard_Natural_Vectors.Vector
     --            := node.ps.white(node.ps.white'last).rows;
      startsols : Solution_List;

      use DoblDobl_Complex_Matrices;
      use Checker_Posets,Checker_Poset_Deformations;

    begin
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          Copy(snd.sols,startsols);
          if tmfo /= null then
            Transform_Start_Solutions
              (n,k,childnode.rows,childnode.cols,
               gamenode.rows,gamenode.cols,tmfo.all,startsols);
          end if;
          declare
            sols : Solution_List;
           -- totalflags : constant natural32 := natural32(flags'length);
            nbflags : constant integer32 := sps.m - level;
          begin
            Track_All_Paths_in_Poset
              (n,k,nt,node.ps,childconds,minrep,tosqr,
               conds(conds'last-nbflags+1..conds'last),
               flags(flags'last-nbflags+1..flags'last),tol,startsols,sols,
               vrblvl-1);
            Push(sols,parent_snd.sols);
          end;
          Deep_Clear(startsols);
        end if;
        exit when (gamenode.next_sibling = null);
        gamenode := gamenode.next_sibling;
      end loop;
    end Connect_Parent;
    procedure Connect_Parents is
      new Intersection_Posets.Enumerate_Parents(Connect_Parent);

  begin
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Track 5 ...");
    end if;
    Connect_Parents(pl,nd.all);
  end Connect_Checker_Posets_to_Track;

  procedure Connect_Checker_Posets_to_Track
              ( n,k,level,nt : in integer32; tol : in double_float;
                pl : in Poset_List;
                snd : in QuadDobl_Solution_Posets.Link_to_Solution_Node;
                tmfo : in QuadDobl_Complex_Matrices.Link_to_Matrix;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                vrblvl : in integer32 := 0 ) is

    nd : constant Link_to_Poset_Node := snd.lpnd;

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Solution_Posets;

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
     -- childrows : constant Standard_Natural_Vectors.Vector
     --           := Flip(childconds);
      gamenode : Checker_Posets.Link_to_Node
               := node.ps.white(node.ps.white'last);
     -- parentconds : constant Standard_Natural_Vectors.Vector
     --             := node.ps.white(node.ps.white'last).cols; -- leaf
     -- parentrows : constant Standard_Natural_Vectors.Vector
     --            := node.ps.white(node.ps.white'last).rows;
      startsols : Solution_List;

      use QuadDobl_Complex_Matrices;
      use Checker_Posets,Checker_Poset_Deformations;

    begin
      loop
        if Standard_Natural_Vectors.Equal(gamenode.cols,childconds) then
          Add(gamenode.coeff,childnode.coeff);
          Copy(snd.sols,startsols);
          if tmfo /= null then
            Transform_Start_Solutions
              (n,k,childnode.rows,childnode.cols,
               gamenode.rows,gamenode.cols,tmfo.all,startsols);
          end if;
          declare
            sols : Solution_List;
           -- totalflags : constant natural32 := natural32(flags'length);
            nbflags : constant integer32 := sps.m - level;
          begin
            Track_All_Paths_in_Poset
              (n,k,nt,node.ps,childconds,minrep,tosqr,
               conds(conds'last-nbflags+1..conds'last),
               flags(flags'last-nbflags+1..flags'last),tol,startsols,sols,
               vrblvl-1);
            Push(sols,parent_snd.sols);
          end;
          Deep_Clear(startsols);
        end if;
        exit when (gamenode.next_sibling = null);
        gamenode := gamenode.next_sibling;
      end loop;
    end Connect_Parent;
    procedure Connect_Parents is
      new Intersection_Posets.Enumerate_Parents(Connect_Parent);

  begin
    if vrblvl > 0 then
      put("-> in resolve_schubert_problems.");
      put_line("Connect_Checker_Posets_to_Track 6 ...");
    end if;
    Connect_Parents(pl,nd.all);
  end Connect_Checker_Posets_to_Track;

  procedure Count_Roots
              ( ips : in out Intersection_Poset;
                roco : out Natural_Number; vrblvl : in integer32 := 0 ) is
   
    tmp : Poset_List;
    lpn : Link_to_Poset_Node;

  begin
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    for i in reverse 1..ips.m loop
      tmp := ips.nodes(i);
      for j in 1..Length_Of(ips.nodes(i)) loop
        lpn := Head_Of(tmp);
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        if i > 1 then
          Connect_Checker_Posets_to_Count(ips.nodes(i-1),lpn.all,vrblvl-1);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    Copy(lpn.ps.white(lpn.ps.white'first).coeff,roco);
  end Count_Roots;

  procedure Count_Roots
              ( file : in file_type; ips : in out Intersection_Poset;
                roco : out Natural_Number; vrblvl : in integer32 := 0 ) is
   
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
          Connect_Checker_Posets_to_Count
            (file,ips.nodes(i-1),lpn.all,vrblvl-1);
        end if;
        tmp := Tail_Of(tmp);
      end loop;
    end loop;
    Copy(lpn.ps.white(lpn.ps.white'first).coeff,roco);
  end Count_Roots;

  procedure Resolve
              ( file : in file_type; extopt,repcon : in boolean;
                n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                sols : out Standard_Complex_Solutions.Solution_List;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    use Standard_Solution_Posets;
    use Flag_Transformations;

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;
    residual : double_float;
    trans : Standard_Complex_Matrices.Link_to_Matrix := null;
    stack : standard_stack_of_flags(flags'first..flags'last-1);
    index : integer32 := stack'last;
    workf : Standard_Complex_VecMats.Link_to_VecMat;
    sqA,sqinvA,sqT : Standard_Complex_VecMats.VecMat(stack'range);

  begin
    if vrblvl > 0
     then put_line("-> in resolve_schubert_problems.Resolve 1 ...");
    end if;
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    if flags'last = 1 then
      Initialize_Solution_Nodes
        (file,n,k,conds(conds'last..conds'last),
         flags(flags'last..flags'last),sps.nodes(sps.m),residual);
    elsif flags'last > 1 then
      Flag_Transformations.Create(n,flags,stack,sqA,sqinvA,sqT);
      put_line(file,"Successively transformed a sequence of flags ...");
      put_line(file,"The stack of flags :");
      for i in stack'range loop
        put(file,"the stack at level "); put(file,i,1);
        put_line(file," :");
        for k in stack(i)'range loop
          put(file,"matrix "); put(file,k,1); put_line(file,":");
          put(file,stack(i)(k).all);
        end loop;
      end loop;
      workf := stack(index);
      Initialize_Solution_Nodes
        (file,n,k,conds(conds'last..conds'last),
         workf(workf'last..workf'last),sps.nodes(sps.m),residual);
    else -- flags'last < 1 ???
      return;
    end if;
    for i in reverse 1..sps.m loop
      new_line(file);
      put(file,"Solving at level "); put(file,i,1); put_line(file," :");
      tmp := sps.nodes(i);
      for j in 1..Length_Of(sps.nodes(i)) loop
        snd := Head_Of(tmp);
        lpn := snd.lpnd;
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        put(file,"-> Level "); put(file,i,1); 
        put(file," : poset node "); put(file,j,1);
        put_line(file,", root and leaves :");
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'first);
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'last);
        put(file,"*** number of paths tracking in checker game : ");
        put(file,lpn.ps.white(lpn.ps.white'first).coeff);
        put_line(file," ***");
        if i > 1 then
          put_line(file,"-> solving at the leaves of its parents :");
          if extopt or repcon then
            if i = 2 then -- use the original flags
              Connect_Checker_Posets_to_Track
                (file,n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 verify,minrep,tosqr,conds,flags,rpt,vrblvl-1);
            else
              Connect_Checker_Posets_to_Track
                (file,n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 verify,minrep,tosqr,conds,workf.all,rpt,vrblvl-1);
            end if;
          else
            if i = 2 then -- use the original flags
              Connect_Checker_Posets_to_Track
                (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 minrep,tosqr,conds,flags,rpt,vrblvl-1);
            else
              Connect_Checker_Posets_to_Track
                (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 minrep,tosqr,conds,workf.all,rpt,vrblvl-1);
            end if;
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      if flags'last > 1 then
        trans := sqT(index);
        if index > 1 then
          index := index - 1;
          workf := stack(index);
        end if;
      end if;
      sps.level := i; -- note that level of completion goes in reverse!
    end loop;
    put(file,"The formal root count : ");
    put(file,lpn.ps.white(lpn.ps.white'first).coeff);
    new_line(file);
    snd := Head_Of(sps.nodes(sps.level));
    sols := snd.sols;
  end Resolve;

  procedure Resolve
              ( file : in file_type; extopt,repcon : in boolean;
                n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Solution_Posets;
    use Flag_Transformations;

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;
    residual : double_double;
    trans : DoblDobl_Complex_Matrices.Link_to_Matrix := null;
    stack : dobldobl_stack_of_flags(flags'first..flags'last-1);
    index : integer32 := stack'last;
    workf : DoblDobl_Complex_VecMats.Link_to_VecMat;
    sqA,sqinvA,sqT : DoblDobl_Complex_VecMats.VecMat(stack'range);

  begin
    if vrblvl > 0
     then put_line("-> in resolve_schubert_problems.Resolve 2 ...");
    end if;
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    if flags'last = 1 then
      Initialize_Solution_Nodes
        (file,n,k,conds(conds'last..conds'last),
         flags(flags'last..flags'last),sps.nodes(sps.m),residual);
    elsif flags'last > 1 then
      Flag_Transformations.Create(n,flags,stack,sqA,sqinvA,sqT);
      workf := stack(index);
      Initialize_Solution_Nodes
        (file,n,k,conds(conds'last..conds'last),
         workf(workf'last..workf'last),sps.nodes(sps.m),residual);
    else -- flags'last < 1 ???
      return;
    end if;
    for i in reverse 1..sps.m loop
      new_line(file);
      put(file,"Solving at level "); put(file,i,1); put_line(file," :");
      tmp := sps.nodes(i);
      for j in 1..Length_Of(sps.nodes(i)) loop
        snd := Head_Of(tmp);
        lpn := snd.lpnd;
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        put(file,"-> Level "); put(file,i,1); 
        put(file," : poset node "); put(file,j,1);
        put_line(file,", root and leaves :");
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'first);
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'last);
        put(file,"*** number of paths tracking in checker game : ");
        put(file,lpn.ps.white(lpn.ps.white'first).coeff);
        put_line(file," ***");
        if i > 1 then
          put_line(file,"-> solving at the leaves of its parents :");
          if extopt or repcon then
            if i = 2 then -- use the original flags
              Connect_Checker_Posets_to_Track
                (file,n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 verify,minrep,tosqr,conds,flags,vrblvl-1);
            else
              Connect_Checker_Posets_to_Track
                (file,n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 verify,minrep,tosqr,conds,workf.all,vrblvl-1);
            end if;
          else
            if i = 2 then -- use the original flags
              Connect_Checker_Posets_to_Track
                (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 minrep,tosqr,conds,flags,vrblvl-1);
            else
              Connect_Checker_Posets_to_Track
                (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 minrep,tosqr,conds,workf.all,vrblvl-1);
            end if;
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      if flags'last > 1 then
        trans := sqT(index);
        if index > 1 then
          index := index - 1;
          workf := stack(index);
        end if;
      end if;
      sps.level := i; -- note that level of completion goes in reverse!
    end loop;
    put(file,"The formal root count : ");
    put(file,lpn.ps.white(lpn.ps.white'first).coeff);
    new_line(file);
    snd := Head_Of(sps.nodes(sps.level));
    sols := snd.sols;
  end Resolve;

  procedure Resolve
              ( file : in file_type; extopt,repcon : in boolean;
                n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                verify,minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Solution_Posets;
    use Flag_Transformations;

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;
    residual : quad_double;
    trans : QuadDobl_Complex_Matrices.Link_to_Matrix := null;
    stack : QuadDobl_stack_of_flags(flags'first..flags'last-1);
    index : integer32 := stack'last;
    workf : QuadDobl_Complex_VecMats.Link_to_VecMat;
    sqA,sqinvA,sqT : QuadDobl_Complex_VecMats.VecMat(stack'range);

  begin
    if vrblvl > 0
     then put_line("-> in resolve_schubert_problems.Resolve 3 ...");
    end if;
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    if flags'last = 1 then
      Initialize_Solution_Nodes
        (file,n,k,conds(conds'last..conds'last),
         flags(flags'last..flags'last),sps.nodes(sps.m),residual);
    elsif flags'last > 1 then
      Flag_Transformations.Create(n,flags,stack,sqA,sqinvA,sqT);
      workf := stack(index);
      Initialize_Solution_Nodes
        (file,n,k,conds(conds'last..conds'last),
         workf(workf'last..workf'last),sps.nodes(sps.m),residual);
    else -- flags'last < 1 ???
      return;
    end if;
    for i in reverse 1..sps.m loop
      new_line(file);
      put(file,"Solving at level "); put(file,i,1); put_line(file," :");
      tmp := sps.nodes(i);
      for j in 1..Length_Of(sps.nodes(i)) loop
        snd := Head_Of(tmp);
        lpn := snd.lpnd;
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        put(file,"-> Level "); put(file,i,1); 
        put(file," : poset node "); put(file,j,1);
        put_line(file,", root and leaves :");
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'first);
        Write_Nodes_in_Poset(file,lpn.ps,lpn.ps.black'last);
        put(file,"*** number of paths tracking in checker game : ");
        put(file,lpn.ps.white(lpn.ps.white'first).coeff);
        put_line(file," ***");
        if i > 1 then
          put_line(file,"-> solving at the leaves of its parents :");
          if extopt or repcon then
            if i = 2 then -- use the original flags
              Connect_Checker_Posets_to_Track
                (file,n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 verify,minrep,tosqr,conds,flags,vrblvl-1);
            else
              Connect_Checker_Posets_to_Track
                (file,n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 verify,minrep,tosqr,conds,workf.all,vrblvl-1);
            end if;
          else
            if i = 2 then -- use the original flags
              Connect_Checker_Posets_to_Track
                (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 minrep,tosqr,conds,flags,vrblvl-1);
            else
              Connect_Checker_Posets_to_Track
                (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
                 minrep,tosqr,conds,workf.all,vrblvl-1);
            end if;
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      if flags'last > 1 then
        trans := sqT(index);
        if index > 1 then
          index := index - 1;
          workf := stack(index);
        end if;
      end if;
      sps.level := i; -- note that level of completion goes in reverse!
    end loop;
    put(file,"The formal root count : ");
    put(file,lpn.ps.white(lpn.ps.white'first).coeff);
    new_line(file);
    snd := Head_Of(sps.nodes(sps.level));
    sols := snd.sols;
  end Resolve;

  procedure Resolve
              ( n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out Standard_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in Standard_Complex_VecMats.VecMat;
                sols : out Standard_Complex_Solutions.Solution_List;
                rpt : in boolean := true; vrblvl : in integer32 := 0 ) is

    use Standard_Solution_Posets;
    use Flag_Transformations;

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;
    residual : double_float;
    trans : Standard_Complex_Matrices.Link_to_Matrix := null;
    stack : standard_stack_of_flags(flags'first..flags'last-1);
    index : integer32 := stack'last;
    workf : Standard_Complex_VecMats.Link_to_VecMat;
    sqA,sqinvA,sqT : Standard_Complex_VecMats.VecMat(stack'range);

  begin
    if vrblvl > 0
     then put_line("-> in resolve_schubert_problems.Resolve 4 ...");
    end if;
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    if flags'last = 1 then
      Initialize_Solution_Nodes
        (n,k,conds(conds'last..conds'last),
         flags(flags'last..flags'last),sps.nodes(sps.m),residual);
    elsif flags'last > 1 then
      Flag_Transformations.Create(n,flags,stack,sqA,sqinvA,sqT);
      workf := stack(index);
      Initialize_Solution_Nodes
        (n,k,conds(conds'last..conds'last),
         workf(workf'last..workf'last),sps.nodes(sps.m),residual);
    else -- flags'last < 1 ???
      return;
    end if;
    for i in reverse 1..sps.m loop
      tmp := sps.nodes(i);
      for j in 1..Length_Of(sps.nodes(i)) loop
        snd := Head_Of(tmp);
        lpn := snd.lpnd;
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        if i > 1 then
          if i = 2 then -- use the original flags
            Connect_Checker_Posets_to_Track
              (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
               minrep,tosqr,conds,flags,rpt,vrblvl-1);
          else
            Connect_Checker_Posets_to_Track
              (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
               minrep,tosqr,conds,workf.all,rpt,vrblvl-1);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      if flags'last > 1 then
        trans := sqT(index);
        if index > 1 then
          index := index - 1;
          workf := stack(index);
        end if;
      end if;
      sps.level := i; -- note that level of completion goes in reverse!
    end loop;
    snd := Head_Of(sps.nodes(sps.level));
    sols := snd.sols;
  end Resolve;

  procedure Resolve
              ( n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out DoblDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in DoblDobl_Complex_VecMats.VecMat;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use DoblDobl_Solution_Posets;
    use Flag_Transformations;

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;
    residual : double_double;
    trans : DoblDobl_Complex_Matrices.Link_to_Matrix := null;
    stack : DoblDobl_Stack_of_Flags(flags'first..flags'last-1);
    index : integer32 := stack'last;
    workf : DoblDobl_Complex_VecMats.Link_to_VecMat;
    sqA,sqinvA,sqT : DoblDobl_Complex_VecMats.VecMat(stack'range);

  begin
    if vrblvl > 0
     then put_line("-> in resolve_schubert_problems.Resolve 5 ...");
    end if;
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    if flags'last = 1 then
      Initialize_Solution_Nodes
        (n,k,conds(conds'last..conds'last),
         flags(flags'last..flags'last),sps.nodes(sps.m),residual);
    elsif flags'last > 1 then
      Flag_Transformations.Create(n,flags,stack,sqA,sqinvA,sqT);
      workf := stack(index);
      Initialize_Solution_Nodes
        (n,k,conds(conds'last..conds'last),
         workf(workf'last..workf'last),sps.nodes(sps.m),residual);
    else -- flags'last < 1 ???
      return;
    end if;
    for i in reverse 1..sps.m loop
      tmp := sps.nodes(i);
      for j in 1..Length_Of(sps.nodes(i)) loop
        snd := Head_Of(tmp);
        lpn := snd.lpnd;
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        if i > 1 then
          if i = 2 then -- use the original flags
            Connect_Checker_Posets_to_Track
              (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
               minrep,tosqr,conds,flags,vrblvl-1);
          else
            Connect_Checker_Posets_to_Track
              (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
               minrep,tosqr,conds,workf.all,vrblvl-1);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      if flags'last > 1 then
        trans := sqT(index);
        if index > 1 then
          index := index - 1;
          workf := stack(index);
        end if;
      end if;
      sps.level := i; -- note that level of completion goes in reverse!
    end loop;
    snd := Head_Of(sps.nodes(sps.level));
    sols := snd.sols;
  end Resolve;

  procedure Resolve
              ( n,k,nt : in integer32; tol : in double_float;
                ips : in out Intersection_Poset;
                sps : in out QuadDobl_Solution_Posets.Solution_Poset;
                minrep,tosqr : in boolean;
                conds : in Standard_Natural_VecVecs.VecVec;
                flags : in QuadDobl_Complex_VecMats.VecMat;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                vrblvl : in integer32 := 0 ) is

    use QuadDobl_Solution_Posets;
    use Flag_Transformations;

    tmp : Solnode_List;
    snd : Link_to_Solution_Node;
    lpn : Link_to_Poset_Node;
    residual : quad_double;
    trans : QuadDobl_Complex_Matrices.Link_to_Matrix := null;
    stack : QuadDobl_Stack_of_Flags(flags'first..flags'last-1);
    index : integer32 := stack'last;
    workf : QuadDobl_Complex_VecMats.Link_to_VecMat;
    sqA,sqinvA,sqT : QuadDobl_Complex_VecMats.VecMat(stack'range);

  begin
    if vrblvl > 0
     then put_line("-> in resolve_schubert_problems.Resolve 6 ...");
    end if;
    Initialize_Leaves(ips.nodes(ips.m));
    for i in 1..ips.m-1 loop
      Initialize_Nodes(ips.nodes(i));
    end loop;
    if flags'last = 1 then
      Initialize_Solution_Nodes
        (n,k,conds(conds'last..conds'last),
         flags(flags'last..flags'last),sps.nodes(sps.m),residual);
    elsif flags'last > 1 then
      Flag_Transformations.Create(n,flags,stack,sqA,sqinvA,sqT);
      workf := stack(index);
      Initialize_Solution_Nodes
        (n,k,conds(conds'last..conds'last),
         workf(workf'last..workf'last),sps.nodes(sps.m),residual);
    else -- flags'last < 1 ???
      return;
    end if;
    for i in reverse 1..sps.m loop
      tmp := sps.nodes(i);
      for j in 1..Length_Of(sps.nodes(i)) loop
        snd := Head_Of(tmp);
        lpn := snd.lpnd;
        Checker_Posets.Add_from_Leaves_to_Root(lpn.ps);
        if i > 1 then
          if i = 2 then -- use the original flags
            Connect_Checker_Posets_to_Track
              (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
               minrep,tosqr,conds,flags,vrblvl-1);
          else
            Connect_Checker_Posets_to_Track
              (n,k,i-1,nt,tol,ips.nodes(i-1),snd,trans,sps,
               minrep,tosqr,conds,workf.all,vrblvl-1);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      if flags'last > 1 then
        trans := sqT(index);
        if index > 1 then
          index := index - 1;
          workf := stack(index);
        end if;
      end if;
      sps.level := i; -- note that level of completion goes in reverse!
    end loop;
    snd := Head_Of(sps.nodes(sps.level));
    sols := snd.sols;
  end Resolve;

end Resolve_Schubert_Problems;
