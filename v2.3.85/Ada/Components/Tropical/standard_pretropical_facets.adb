with text_io;                            use text_io;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Integer64_Vectors_io;      use Standard_Integer64_Vectors_io;
with Standard64_Common_Divisors;
with Standard_Integer_Norms;
with Lists_of_Integer64_Vectors_io;      use Lists_of_Integer64_Vectors_io;
with Standard_Lattice_Supports;
with Standard_Lattice_3d_Facets_io;      use Standard_Lattice_3d_Facets_io;

package body Standard_Pretropical_Facets is

-- converting supports into matrix formats :

  function List2Matrix ( A : Lists_of_Integer_Vectors.List ) return Matrix is

    tmp : Lists_of_Integer_Vectors.List := A;
    ls : Standard_Integer_Vectors.Link_to_Vector
       := Lists_of_Integer_Vectors.Head_Of(tmp);
    n : constant integer32 := ls'last;
    m : constant integer32 := integer32(Lists_of_Integer_Vectors.Length_Of(A));
    res : Matrix(1..n,1..m);

  begin
    for j in 1..m loop
      ls := Lists_of_Integer_Vectors.Head_Of(tmp);
      for i in 1..n loop
        res(i,j) := integer64(ls(i));
      end loop;
      tmp := Lists_of_Integer_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end List2Matrix;

  function Lists2VecMat ( A : Array_of_Lists ) return VecMat is

    res : VecMat(A'range);

  begin
    for i in A'range loop
      if not Lists_of_Integer_Vectors.Is_Null(A(i))
       then res(i) := new Matrix'(List2Matrix(A(i)));
      end if;
    end loop;
    return res;
  end Lists2VecMat;

-- computing pretropical facets :

  function Facet_Pretropisms
              ( f,g : Facet_3d_List ) return Lists_of_Integer64_Vectors.List is

    res,res_last : Lists_of_Integer64_Vectors.List;
    f_ptr : Facet_3d_List := f;
    g_ptr : Facet_3d_List;
    f_facet,g_facet : Link_to_3d_Facet;

  begin
    while not Is_Null(f_ptr) loop
      f_facet := Head_Of(f_ptr);
      g_ptr := g;
      while not Is_Null(g_ptr) loop
        g_facet := Head_Of(g_ptr);
        if Standard_Integer64_Vectors.Equal(f_facet.normal,g_facet.normal)
         then Lists_of_Integer64_Vectors.Append(res,res_last,f_facet.normal);
        end if;
        g_ptr := Tail_Of(g_ptr);
      end loop;
      f_ptr := Tail_Of(f_ptr);
    end loop;
    return res;
  end Facet_Pretropisms;

  function Facet_Pretropisms
              ( f : Array_of_Facet_3d_Lists )
              return Lists_of_Integer64_Vectors.List is

    res : Lists_of_Integer64_Vectors.List
        := Facet_Pretropisms(f(f'first),f(f'first+1));

  begin
    for i in f'first+2..f'last loop
      declare
        newres : constant Lists_of_Integer64_Vectors.List
               := Select_Facet_Normals(f(i),res);
      begin
        Lists_of_Integer64_Vectors.Clear(res);
        res := newres;
      end;
    end loop;
    return res;
  end Facet_Pretropisms;

  function Pretropical_Facets
               ( f : Array_of_Facet_3d_Lists;
                 v : Lists_of_Integer64_Vectors.List ) 
               return Array_of_Facet_3d_Lists is

    res,res_last : Array_of_Facet_3d_Lists(f'range);
    tmp : Lists_of_Integer64_Vectors.List := v;
    lv : Standard_Integer64_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer64_Vectors.Is_Null(tmp) loop
      lv := Lists_of_Integer64_Vectors.Head_Of(tmp);
      for i in res'range loop
        Append(res(i),res_last(i),Get_Facet(f(i),lv.all));
      end loop;
      tmp := Lists_of_Integer64_Vectors.Tail_Of(tmp);
    end loop;
    return res;
  end Pretropical_Facets;

  function Is_Subset
             ( p,s : Standard_Integer_Vectors.Vector ) return boolean is
  begin
    for i in s'range loop
      if Standard_Lattice_Supports.Member(p,s(i)) < p'first
       then return false;
      end if;
    end loop;
    return true;
  end Is_Subset;

  function On_Pretropical_Edge
                ( A,B : Matrix; e,v : Standard_Integer64_Vectors.Vector )
               return boolean is

    eA : constant Standard_Integer_Vectors.Vector
       := Standard_Lattice_Supports.Support(A,e);
    vA : constant Standard_Integer_Vectors.Vector
       := Standard_Lattice_Supports.Support(A,v);

  begin
    if not Is_Subset(eA,vA) then
      return false;
    else
      declare
        eB : constant Standard_Integer_Vectors.Vector
           := Standard_Lattice_Supports.Support(B,e);
        vB : constant Standard_Integer_Vectors.Vector
           := Standard_Lattice_Supports.Support(B,v);
      begin
        return Is_Subset(eB,vB);
      end;
    end if;
  end On_Pretropical_Edge;

  function On_Pretropical_Edge
                ( A,B : Matrix; e : Lists_of_Integer64_Vectors.List;
                  v : Standard_Integer64_Vectors.Vector )
               return boolean is

    tmp : Lists_of_Integer64_Vectors.List := e;
    lft : Standard_Integer64_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer64_Vectors.Is_Null(tmp) loop
      lft := Lists_of_Integer64_Vectors.Head_Of(tmp);
      if On_Pretropical_Edge(A,B,lft.all,v)
       then return true;     
       else tmp := Lists_of_Integer64_Vectors.Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end On_Pretropical_Edge;

  function Subset_Count
             ( p,s : Standard_Integer_Vectors.Vector ) return natural32 is

    res : natural32 := 0;

  begin
    for i in s'range loop
      if Standard_Lattice_Supports.Member(p,s(i)) > p'first-1
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Subset_Count;

  function On_Pretropical_Facet
             ( A,B : Matrix; f,v : Standard_Integer64_Vectors.Vector )
             return boolean is

    fA : constant Standard_Integer_Vectors.Vector
       := Standard_Lattice_Supports.Support(A,f);
    vA : constant Standard_Integer_Vectors.Vector
       := Standard_Lattice_Supports.Support(A,v);

  begin
    if Subset_Count(fA,vA) < 2 then
      return false;
    else
      declare
        fB : constant Standard_Integer_Vectors.Vector
           := Standard_Lattice_Supports.Support(B,f);
        vB : constant Standard_Integer_Vectors.Vector
           := Standard_Lattice_Supports.Support(B,v);
      begin
        if Subset_Count(fB,vB) < 2
         then return false;
         else return true;
        end if;
      end;
    end if;
  end On_Pretropical_Facet;

  function On_Pretropical_Facet
             ( A,B : Matrix; f : Lists_of_Integer64_Vectors.List;
               v : Standard_Integer64_Vectors.Vector )
             return boolean is

    tmp : Lists_of_Integer64_Vectors.List := f;
    lft : Standard_Integer64_Vectors.Link_to_Vector;

  begin
    while not Lists_of_Integer64_Vectors.Is_Null(tmp) loop
      lft := Lists_of_Integer64_Vectors.Head_Of(tmp);
      if On_Pretropical_Facet(A,B,lft.all,v)
       then return true;     
       else tmp := Lists_of_Integer64_Vectors.Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end On_Pretropical_Facet;

  function On_Pretropical_Facet
             ( A : VecMat; f : Array_of_3d_Facets;
               v : Standard_Integer64_Vectors.Vector )
             return boolean is

  begin
    for i in f'range loop
      declare
        s : constant Standard_Integer_Vectors.Vector
          := Standard_Lattice_Supports.Support(A(i).all,v);
      begin
        if s'length < 2 then
          return false;
        elsif Subset_Count(f(i).points,s) < 2 then
          return false;
        end if;
      end;
    end loop;
    return true;
  end On_Pretropical_Facet;

  function On_Pretropical_Facet
             ( A : VecMat; f : Array_of_Facet_3d_Lists;
               v : Standard_Integer64_Vectors.Vector )
             return boolean is

    tmp : Array_of_Facet_3d_Lists(f'range) := f;
    ftp : Array_of_3d_Facets(f'range);

  begin
    while not Is_Null(tmp(tmp'first)) loop
      for i in tmp'range loop
        ftp(i) := Head_Of(tmp(i));
      end loop;
      if On_Pretropical_Facet(A,ftp,v)
       then return true;
      end if;
      for i in tmp'range loop
        tmp(i) := Tail_Of(tmp(i));
      end loop;
    end loop;     
    return false;
  end On_Pretropical_Facet;

-- computing pretropical edges by Minkowski sum :

  function Minkowski_Sum ( A,B : Matrix ) return Matrix is

    res : Matrix(A'range(1),1..A'last(2)*B'last(2));
    ind : integer32 := 0;

  begin
    for i in A'range(2) loop
      for j in B'range(2) loop
        ind := ind + 1;
        for k in A'range(1) loop
          res(k,ind) := A(k,i) + B(k,j);
        end loop;
      end loop;
    end loop;
    return res;
  end Minkowski_Sum;

  function Is_Mixed_Facet_Normal
             ( A,B : Matrix; v : Standard_Integer64_Vectors.Vector )
             return boolean is

    sA : constant Standard_Integer_Vectors.Vector
       := Standard_Lattice_Supports.Support(A,v);
    sB : constant Standard_Integer_Vectors.Vector
       := Standard_Lattice_Supports.Support(B,v);

  begin
    return ((sA'last > 1) and (sB'last > 1));
  end Is_Mixed_Facet_Normal;

  function Mixed_Facet_Normals
              ( A,B : Matrix; f : Facet_3d_List )
              return Lists_of_Integer64_Vectors.List is

    res,res_last : Lists_of_Integer64_Vectors.List;
    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if Is_Mixed_Facet_Normal(A,B,lft.normal)
       then Lists_of_Integer64_Vectors.Append(res,res_last,lft.normal);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Mixed_Facet_Normals;

  function Edge_Tropisms_by_Sum
              ( A,B : Matrix ) return Lists_of_Integer64_Vectors.List is

    res : Lists_of_Integer64_Vectors.List;
    C : constant Matrix := Minkowski_Sum(A,B);
    f : Facet_3d_List := Convex_Hull_3D(C);

  begin
    res := Mixed_Facet_Normals(A,B,f);
    Clear(f);
    return res;
  end Edge_Tropisms_by_Sum;

-- pretropical edges of two supports via wrapping :

  procedure Check_Edge
              ( a,b,u,v : in Standard_Integer64_Vectors.Vector;
                uk,vk : out integer64 ) is

  -- DESCRIPTION :
  --   Computes coefficients uk and vk with u and v, so w = uk*u + vk*v 
  --   is perpendicular to the edge spanned by a and b.

    use Standard_Integer64_Vectors;

    amb : constant Vector := a - b;
    alpha : integer64 := u(1)*amb(1) + u(2)*amb(2) + u(3)*amb(3);
    beta : integer64 := v(1)*amb(1) + v(2)*amb(2) + v(3)*amb(3);
    g : constant integer64 := Standard64_Common_Divisors.gcd(alpha,beta);

  begin
    if g /= 0 and g /= 1
     then alpha := alpha/g; beta := beta/g;
    end if;
    if alpha = 0 then
      uk := 1; vk := 0;
    elsif beta = 0 then
      uk := 0; vk := 1;
    elsif beta > 0 then
      uk := beta; vk := -alpha;
    else
      uk := -beta; vk := alpha;
    end if;
  end Check_Edge;

  function Is_Tropism 
             ( A,B : Matrix;
               v : Standard_Integer64_Vectors.Vector ) return boolean is

    sA : constant Standard_Integer_Vectors.Vector
       := Standard_Lattice_Supports.Support(A,v);

  begin
    if sA'length < 2 then
      return false;
    else
      declare
        sB : constant Standard_Integer_Vectors.Vector
           := Standard_Lattice_Supports.Support(B,v);
      begin
        return (sB'length > 1);
      end;
    end if;
  end Is_Tropism;

  procedure Wrap_Edge_for_Tropism
              ( A,B : in Matrix; ip,iq : in integer32;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List;
                fail : out boolean ) is

    use Standard_Integer64_Vectors;
    p,q,t : Standard_Integer64_Vectors.Vector(B'range(1));
    uk,vk : integer64;

  begin
    put("testing edge ("); put(ip,1); put(","); put(iq,1); put(")");
    for i in p'range loop
      p(i) := B(i,ip);
      q(i) := B(i,iq);
    end loop;
    Check_Edge(p,q,u,v,uk,vk);
    put("  uk = "); put(uk,1); put("  vk = "); put(vk,1);
    if uk >= 0 and vk >= 0 then
      t := uk*u + vk*v; 
      Standard_Integer_Norms.Normalize(t);
      new_line;
      put("  the candidate tropism : "); put(t); 
      put(" supports (");
      put(Standard_Lattice_Supports.Support(A,t)); put(",");
      put(Standard_Lattice_Supports.Support(B,t)); put(")");
      if Is_Tropism(A,B,t) then
        fail := false; put_line("  is a tropism");
        if not Lists_of_Integer64_Vectors.Is_In(trp,t) then
          if not On_Pretropical_Edge(A,B,trp,t) 
           then Lists_of_Integer64_Vectors.Append(trp,trp_last,t);
          end if;
        end if;
      else
        fail := true; put_line("  is not a tropism");
      end if;
    else
      fail := true; put_line("  no tropism");
    end if;
  end Wrap_Edge_for_Tropism;

  procedure Wrap_Edge_for_Tropism
              ( A,B : in Matrix; ip,iq : in integer32;
                fpt : in Lists_of_Integer64_Vectors.List;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List;
                fail : out boolean ) is

    use Standard_Integer64_Vectors;
    p,q,t : Standard_Integer64_Vectors.Vector(B'range(1));
    uk,vk : integer64;

  begin
    put("testing edge ("); put(ip,1); put(","); put(iq,1); put(")");
    for i in p'range loop
      p(i) := B(i,ip);
      q(i) := B(i,iq);
    end loop;
    Check_Edge(p,q,u,v,uk,vk);
    put("  uk = "); put(uk,1); put("  vk = "); put(vk,1);
    if uk >= 0 and vk >= 0 then
      t := uk*u + vk*v; 
      Standard_Integer_Norms.Normalize(t);
      new_line;
      put("  the candidate tropism : "); put(t); 
      put(" supports (");
      put(Standard_Lattice_Supports.Support(A,t)); put(",");
      put(Standard_Lattice_Supports.Support(B,t)); put(")");
      if Is_Tropism(A,B,t) then
        fail := false;
        if  Lists_of_Integer64_Vectors.Is_In(trp,t) then
          put_line(" old tropism");
        else
          if On_Pretropical_Facet(A,B,fpt,t) then
            put_line(" on pretropical facet");
          elsif On_Pretropical_Edge(A,B,trp,t) then
            put_line(" on pretropical edge");
          else
            put_line(" is proper edge tropism"); 
            Lists_of_Integer64_Vectors.Append(trp,trp_last,t);
          end if;
        end if;
      else
        fail := true; put_line("  is not a tropism");
      end if;
    else
      fail := true; put_line("  no tropism");
    end if;
  end Wrap_Edge_for_Tropism;

  procedure Wrap_Edges
              ( A,B : in Matrix; f : in Facet_3d_List;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List ) is

    use Standard_Integer64_Vectors;
    sum : Standard_Integer64_Vectors.Vector(u'range) := u + v;
    sup : constant Standard_Integer_Vectors.Vector
        := Standard_Lattice_Supports.Support(B,sum);
    ind : constant integer32 := sup(sup'first);
    bom : Boolean_Matrix(B'range(2),B'range(2));

    procedure Test_Edge ( p,q : in integer32; continue : out boolean ) is

      notrop : boolean;

    begin
      Wrap_Edge_for_Tropism(A,B,p,q,u,v,trp,trp_last,notrop);
      continue := not notrop;
    end Test_Edge;
    procedure Test_Edges is new Crawl(Test_Edge);

  begin
    put("sum of u ="); put(u);
    put(" and v ="); put(v); put(" ="); put(sum);
    put(" supports"); put(sup);
    if sup'length > 1 then
      put_line(" found tropism!");
      Standard_Integer_Norms.Normalize(sum);
      if not Lists_of_Integer64_Vectors.Is_In(trp,sum) then
        if not On_Pretropical_Edge(A,B,trp,sum)
         then Lists_of_Integer64_Vectors.Append(trp,trp_last,sum);
        end if;
      end if;
    else
      put_line(" no tropism");
    end if;
    for i in bom'range(1) loop
      for j in bom'range(2) loop
        bom(i,j) := false;
      end loop;
    end loop;
    Test_Edges(f,ind,bom);
  end Wrap_Edges;

  procedure Wrap_Edges
              ( A,B : in Matrix; f : in Facet_3d_List;
                fpt : in Lists_of_Integer64_Vectors.List;
                u,v : in Standard_Integer64_Vectors.Vector;
                trp,trp_last : in out Lists_of_Integer64_Vectors.List ) is

    use Standard_Integer64_Vectors;
    sum : Standard_Integer64_Vectors.Vector(u'range) := u + v;
    sup : constant Standard_Integer_Vectors.Vector
        := Standard_Lattice_Supports.Support(B,sum);
    ind : constant integer32 := sup(sup'first);
    bom : Boolean_Matrix(B'range(2),B'range(2));

    procedure Test_Edge ( p,q : in integer32; continue : out boolean ) is

      notrop : boolean;

    begin
      Wrap_Edge_for_Tropism(A,B,p,q,fpt,u,v,trp,trp_last,notrop);
      continue := not notrop;
    end Test_Edge;
    procedure Test_Edges is new Crawl(Test_Edge);

  begin
    put("sum of u ="); put(u);
    put(" and v ="); put(v); put(" ="); put(sum);
    put(" supports"); put(sup);
    if sup'length > 1 then
      Standard_Integer_Norms.Normalize(sum);
      if Lists_of_Integer64_Vectors.Is_In(trp,sum) then
        put_line(" found old tropism");
      elsif On_Pretropical_Facet(A,B,fpt,sum) then
        put_line(" tropism on pretropical facet");
      elsif On_Pretropical_Edge(A,B,trp,sum) then
        put_line(" tropism on pretropical edge");
      else
        put_line(" tropism added to list");
        Lists_of_Integer64_Vectors.Append(trp,trp_last,sum);
      end if;
    else
      put_line(" no tropism");
    end if;
    for i in bom'range(1) loop
      for j in bom'range(2) loop
        bom(i,j) := false;
      end loop;
    end loop;
    Test_Edges(f,ind,bom);
  end Wrap_Edges;

  procedure Visit_Edges
              ( A,B : in Matrix; f : in Link_to_3d_Facet; g : in Facet_3d_List;
                w : in Lists_of_Integer64_Vectors.List;
                t,t_last : in out Lists_of_integer64_Vectors.List;
                cnt : in out natural32 ) is

    j : integer32;

  begin
    for i in f.points'range loop
      if i < f.points'last
       then j := i+1;
       else j := f.points'first;
      end if;
      put("edge "); put(f.points(i),1); put(" "); put(f.points(j),1);
      put(" is neighbor to facet with inner normal");
      put(f.neighbors(i).normal);
      if Lists_of_Integer64_Vectors.Is_In(w,f.neighbors(i).normal) then
        put_line(" visited");
      else
        put_line(" continue"); cnt := cnt + 1;
        Wrap_Edges(A,B,g,f.normal,f.neighbors(i).normal,t,t_last);
      end if;
    end loop;
  end Visit_Edges;

  procedure Visit_Edges
              ( A,B : in Matrix; f : in Link_to_3d_Facet; g : in Facet_3d_List;
                fpt,w : in Lists_of_Integer64_Vectors.List;
                t,t_last : in out Lists_of_integer64_Vectors.List;
                cnt : in out natural32 ) is

    j : integer32;

  begin
    for i in f.points'range loop
      if i < f.points'last
       then j := i+1;
       else j := f.points'first;
      end if;
      put("edge "); put(f.points(i),1); put(" "); put(f.points(j),1);
      put(" is neighbor to facet with inner normal");
      put(f.neighbors(i).normal);
      if Lists_of_Integer64_Vectors.Is_In(w,f.neighbors(i).normal) then
        put_line(" visited");
     -- elsif Lists_of_Integer64_Vectors.Is_In(fpt,f.neighbors(i).normal) then
     --   put_line(" neighbor is pretropical");
      else
        put_line(" continue"); cnt := cnt + 1;
        Wrap_Edges(A,B,g,fpt,f.normal,f.neighbors(i).normal,t,t_last);
      end if;
    end loop;
  end Visit_Edges;

  function Edge_Tropisms_by_Wrap
              ( A,B : Matrix; f,g : Facet_3d_List )
              return Lists_of_Integer64_Vectors.List is

    res,res_last,vis,vis_last : Lists_of_Integer64_Vectors.List;
    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      Visit_Edges(A,B,lft,g,vis,res,res_last,cnt);
      Lists_of_Integer64_Vectors.Append(vis,vis_last,lft.normal);
      tmp := Tail_Of(tmp);
    end loop;
    put("#good edges : "); put(cnt,1); new_line;
    return res;
  end Edge_Tropisms_by_Wrap;

  function Edge_Tropisms_by_Wrap
              ( A,B : Matrix; fpt : Lists_of_Integer64_Vectors.List;
                f,g : Facet_3d_List )
              return Lists_of_Integer64_Vectors.List is

    res,res_last,vis,vis_last : Lists_of_Integer64_Vectors.List;
    tmp : Facet_3d_List := f;
    lft : Link_to_3d_Facet;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if not Lists_of_Integer64_Vectors.Is_In(fpt,lft.normal) then
        Visit_Edges(A,B,lft,g,fpt,vis,res,res_last,cnt);
        Lists_of_Integer64_Vectors.Append(vis,vis_last,lft.normal);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put("#good edges : "); put(cnt,1); new_line;
    return res;
  end Edge_Tropisms_by_Wrap;

  function Edge_Tropisms_by_Wrap
              ( A,B : Matrix ) return Lists_of_Integer64_Vectors.List is

    res : Lists_of_Integer64_Vectors.List;
    Af : Facet_3d_List := Convex_Hull_3D(A);
    Bf : Facet_3d_List := Convex_Hull_3D(B);

  begin
    put_line("Facets of the first polytope : ");
    Write_Facets(A,Af);
    put_line("Facets of the second polytope : ");
    Write_Facets(B,Bf);
    res := Edge_Tropisms_by_Wrap(A,B,Af,Bf);
    Clear(Af); Clear(Bf);
    return res;
  end Edge_Tropisms_by_Wrap;

  procedure Edge_Tropisms_by_Wrap
              ( A,B : in Matrix;
                fpt,ept : out Lists_of_Integer64_Vectors.List ) is

    Af : Facet_3d_List := Convex_Hull_3D(A);
    Bf : Facet_3d_List := Convex_Hull_3D(B);

  begin
    put_line("Facets of the first polytope : ");
    Write_Facets(A,Af);
    put_line("Facets of the second polytope : ");
    Write_Facets(B,Bf);
    fpt := Facet_Pretropisms(Af,Bf);
    if Lists_of_Integer64_Vectors.Is_Null(fpt)
     then put_line("There are no pretropical facets.");
     else put_line("Inner normals to pretropical facets :"); put(fpt);
    end if;
    ept := Edge_Tropisms_by_Wrap(A,B,fpt,Af,Bf);
    Clear(Af); Clear(Bf);
  end Edge_Tropisms_by_Wrap;

-- computing pretropical edges :

  procedure Wrap ( A : in VecMat; g : in Array_of_Facet_3d_Lists;
                   u,v : in Standard_Integer64_Vectors.Vector ) is

    use Standard_Integer64_Vectors;
    sum : constant Standard_Integer64_Vectors.Vector := u + v;
    sup : constant Standard_Integer_Vectors.Vector
        := Standard_Lattice_Supports.Support(A(A'first+1).all,sum);
    done : boolean := false;

  begin
    put("normal"); put(u);
    if On_Pretropical_Facet(A,g,u) then
      put_line(" on pretropical facets...");
      put("normal "); put(v); 
      if On_Pretropical_Facet(A,g,v)
       then put_line(" also on pretropical facet!"); done := true;
       else put_line(" not on any pretropical facet.");
      end if;
    else
      put_line(" not on any pretropical facet.");
    end if;
    if not done then
      put("support of"); put(sum); put(" is "); put(sup); new_line;
      if sup'last > 1 then
        put("support with respect to last polytope : ");
        declare
          sp2 : constant Standard_Integer_Vectors.Vector
              := Standard_Lattice_Supports.Support(A(A'last).all,sum);
        begin
          put(sp2);
          if sp2'last > 1
           then put_line(" edge tropism !");
           else put_line(" no tropism");
          end if;
        end;
      end if;
    end if;
  end Wrap;

  procedure Report_Edges
               ( A : in VecMat; f : in Link_to_3d_Facet;
                 g : in Array_of_Facet_3d_Lists;
                 v,w : in Lists_of_Integer64_Vectors.List;
                 cnt : in out natural32 ) is

    j : integer32;

  begin
    for i in f.points'range loop
      if i < f.points'last
       then j := i+1;
       else j := f.points'first;
      end if;
      put("edge "); put(f.points(i),1); put(" "); put(f.points(j),1);
      put(" is neighbor to facet with inner normal");
      put(f.neighbors(i).normal);
      if Lists_of_Integer64_Vectors.Is_In(v,f.neighbors(i).normal) then
        put_line(" stop");
      elsif Lists_of_Integer64_Vectors.Is_In(w,f.neighbors(i).normal) then
        put_line(" visited");
      else
        put_line(" continue"); cnt := cnt + 1;
        Wrap(A,g,f.normal,f.neighbors(i).normal);
      end if;
    end loop;
  end Report_Edges;

  procedure Edge_Pretropisms
              ( A : in VecMat; f,g : in Array_of_Facet_3d_Lists;
                fpt : in Lists_of_Integer64_Vectors.List ) is

    tmp : Facet_3d_List := f(f'first);
    lft : Link_to_3d_Facet := null;
    vis,vis_last : Lists_of_Integer64_Vectors.List;
    cnt : natural32 := 0;

  begin
    put_line("Facets of first polytope not in facet pretropism list :");
    while not Is_Null(tmp) loop
      lft := Head_Of(tmp);
      if not Lists_of_Integer64_Vectors.Is_In(fpt,lft.normal) then
        Write_Facet(A(A'first).all,lft.all);
        Report_Edges(A,lft,g,fpt,vis,cnt);
        Lists_of_Integer64_Vectors.Append(vis,vis_last,lft.normal);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    put("#Good edges : "); put(cnt,1); new_line;
  end Edge_Pretropisms;

end Standard_Pretropical_Facets;
