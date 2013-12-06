with unchecked_deallocation;
with Standard_Integer_Norms;             use Standard_Integer_Norms;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Standard_Integer32_Transformations;
 use Standard_Integer32_Transformations;

package body Standard_Integer32_Simplices is

-- DATA STRUCTURES :

  type Point is record
    pt : Link_to_Vector;    -- the point
    si : Simplex;           -- by pivoting a new simplex is obtained
  end record;
  type Points is array ( integer32 range <> ) of Point;

  type Simplex_Rep ( n : integer32 ) is record
    nor : vector(1..n);     -- normal to the simplex
    tra : Transfo;          -- transformation to triangulate the cell
    pts : points(1..n);     -- vector of points
  end record;

-- AUXILIAIRIES :

  function Diagonalize ( s : simplex ) return matrix is

  -- DESCRIPTION :
  --   Places the vertices of the simplex, shifted w.r.t. the first point,
  --   in a matrix.  Returns the triangulated matrix.

    a : matrix(1..s.n-1,1..s.n-1);
    x,y : Vector(1..s.n-1);

  begin
    for k in 2..s.n loop
      x := s.pts(k).pt(x'range) - s.pts(1).pt(x'range);
      y := s.tra*x;
      for i in a'range(1) loop
        a(i,k-1) := y(i);
      end loop;
    end loop;
    return a;
  end Diagonalize;

  function Diagonalize ( s : simplex; pt : Vector ) return matrix is

  -- DESCRIPTION :
  --   Places the vertices of the simplex, shifted w.r.t. the first point,
  --   in a matrix.  Returns the triangulated matrix.

    a : matrix(1..s.n-1,1..s.n);
    x,y : Vector(1..s.n-1);

  begin
    for k in 2..s.n loop
      x := s.pts(k).pt(x'range) - s.pts(1).pt(x'range); y := s.tra*x;
      for i in a'range(1) loop
        a(i,k-1) := y(i);
      end loop;
    end loop;
    x := pt(x'range) - s.pts(1).pt(x'range); y := s.tra*x;
    for i in a'range(1) loop
      a(i,s.n) := y(i);
    end loop;
    return a;
  end Diagonalize;

  function Create ( pts : VecVec ) return Transfo is

    a,l : matrix(pts'first..pts'last-1,pts'first..pts'last-1);
    x : vector(pts(pts'first)'range);

  begin
    for k in pts'first+1..pts'last loop
      x := pts(k).all - pts(pts'first).all;
      for i in a'range(1) loop
        a(i,k-1) := x(i);
      end loop;
    end loop;
    Upper_Triangulate(l,a);
    return Create(l);
  end Create;

  function Create ( pts : VecVec ) return Vector is

    a : matrix(pts'first..pts'last - 1,pts'range);
    res : Vector(pts'range);

  begin
    for k in a'range(1) loop
      for i in a'range(2) loop
        a(k,i) := pts(k+1)(i) - pts(pts'first)(i);
      end loop;
    end loop;
    Upper_Triangulate(a);
    Scale(a);
    res := (res'range => 0);
    Solve0(a,res);
    Normalize(res);
    if res(res'last) < 0
     then return -res;
     else return res;
    end if;
  end Create;

-- CREATORS :

  function Create ( x : VecVec ) return Simplex is

    n : constant integer32 := x'last - x'first + 1;
    res : constant Simplex := new Simplex_Rep(n);
    cnt : integer32 := x'first;

  begin
    for k in res.pts'range loop
      res.pts(k).pt := x(cnt);
      cnt := cnt + 1;
      res.pts(k).si := Null_Simplex;
    end loop;
    res.tra := Create(x);
    res.nor := Create(x);
    return res;
  end Create;

  procedure Update ( s : in out Simplex; x : in Link_to_Vector;
                     k : in integer32 ) is

    pts : VecVec(1..s.n);
    nei : Simplex;

  begin
    if s.pts(k).si = Null_Simplex
     then for i in pts'range loop
            if i = k
             then pts(i) := x;
             else pts(i) := s.pts(i).pt;
            end if;
          end loop;
          nei := Create(pts);
          s.pts(k).si := nei;
          nei.pts(k).si := s;
    end if;
  end Update;

  procedure Update ( s : in out Simplex; x : in Link_to_Vector; 
                     pos : in Vector ) is
  begin
    for k in pos'first..pos'last-1 loop
      if pos(k)*pos(pos'last) > 0 
       then Update(s,x,k+1);
      end if;
    end loop;
  end Update;

  procedure Update_One ( s : in out Simplex; x : in Link_to_Vector;
                         pos : in Vector ) is

    done : boolean := false;
    nei : Simplex;

  begin
    for k in pos'first..pos'last-1 loop
      if pos(k)*pos(pos'last) > 0
       then nei := s.pts(k+1).si;
            if nei /= Null_Simplex
             then Update_One(nei,x,Position(nei,x.all));
             else Update(s,x,k+1);
            end if;
            done := true;
      end if;
      exit when done;
    end loop;
  end Update_One;

  procedure Update_One ( s : in out Simplex; x : in Link_to_Vector;
                         pos : in Vector; news : out Simplex ) is

    done : boolean := false;
    nei,newsnei : Simplex;

  begin
   -- LOOK FIRST FOR NULL SIMPLEX IN THE DIRECTION TO x :
    for k in pos'first..pos'last-1 loop
      if pos(k)*pos(pos'last) > 0
       then nei := s.pts(k+1).si;
            if nei = Null_Simplex
             then Update(s,x,k+1);
                  news := s.pts(k+1).si;
                  done := true;
            end if;
      end if;
      exit when done;
    end loop;
   -- WALK FURTHER IN THE DIRECTION TO x :
    if not done
     then Update_One(nei,x,Position(nei,x.all),newsnei);
          if newsnei /= Null_Simplex
           then news := newsnei;
                s := nei;
          end if;
    end if;
  end Update_One;

  procedure Update_All ( s : in out Simplex; x : in Link_to_Vector;
                         pos : in Vector; ancestor : in Simplex ) is

    nei : Simplex;
    continue : boolean := true;

  begin
    for k in pos'first..pos'last-1 loop
      if pos(k)*pos(pos'last) > 0
       then nei := s.pts(k+1).si;
            if nei /= Null_Simplex 
             then if not Is_Vertex(nei,x.all) and nei /= ancestor
                   then Update_All(nei,x,Position(nei,x.all),s);
                  end if;
             else Update(s,x,k+1);
                  Process(s.pts(k+1).si,continue);
            end if;
      end if;
      exit when not continue;
    end loop;
  end Update_All;

  procedure Connect ( s1,s2 : in out Simplex ) is

    neighb : boolean;
    index1,index2 : integer32;

  begin
    neighb := true; -- assume they are neighbors
    index1 := 0; 
   -- SEARCH FOR INDEX OF POINT IN s1 THAT DOES NOT BELONG TO s2 :
    for k in s1.pts'range loop
      if not Is_Vertex(s2,s1.pts(k).pt.all) then
        if (index1 = 0) and (s1.pts(k).si = Null_Simplex)
         then index1 := k;      -- kth point is not common
         else neighb := false;  -- more than one point not common
        end if;                 -- or there is already a neighbor
      end if;
      exit when not neighb;
    end loop;  -- either not neighb or index1 -> point in s1, not in s2
   -- SEARCH FOR INDEX OF POINT IN s2 THAT DOES NOT BELONG TO s1 :
    if neighb
     then index2 := 0;
          for k in s2.pts'range loop
            if not Is_Vertex(s1,s2.pts(k).pt.all)
             then if (index2 = 0) and (s2.pts(k).si = Null_Simplex)
                   then index2 := k;      -- kth point is not common
                   else neighb := false;  -- more than one point not common
                                          -- or there is already a neighbor
                  end if;
            end if;
            exit when not neighb;
          end loop;  -- either no neighb or index2 -> point in s2, not in s1
   -- CONNECT THE SIMPLICES WITH EACH OTHER :
          if neighb
           then s1.pts(index1).si := s2;
                s2.pts(index2).si := s1;
          end if;
    end if;
  end Connect;

  procedure Flatten ( s : in out Simplex ) is
  begin
    s.nor := (s.nor'range => 0);
    s.nor(s.n) := 1;
    for k in s.pts'range loop
      s.pts(k).pt(s.n) := 0;
    end loop;
  end Flatten;

-- SELECTORS :

  function Dimension ( s : Simplex ) return natural32 is
  begin
    return natural32(s.n);
  end Dimension;

  function Normal ( s : Simplex ) return Vector is
  begin
    return s.nor;
  end Normal;

  function Is_Flat ( s : Simplex ) return boolean is
  begin
    for i in s.nor'first..(s.nor'last-1) loop
      if s.nor(i) /= 0
       then return false;
      end if;
    end loop;
    return (s.nor(s.nor'last) = 1);
  end Is_Flat;

  function Vertices ( s : Simplex ) return VecVec is

    res : VecVec(s.pts'range);

  begin
    for k in res'range loop
      res(k) := s.pts(k).pt;
    end loop;
    return res;
  end Vertices;

  function Vertex ( s : Simplex; k : integer32 ) return Vector is
  begin
    return s.pts(k).pt.all;
  end Vertex;

  function Is_Vertex ( s : Simplex; x : Vector ) return boolean is
  begin
    for k in s.pts'range loop
      if s.pts(k).pt.all = x
       then return true;
      end if;
    end loop;
    return false;
  end Is_Vertex;

  function Equal ( s1,s2 : Simplex ) return boolean is

    found : boolean;

  begin
    if s1.nor /= s2.nor                -- check if normals are the same
     then return false;
     else for k in s1.pts'range loop   -- check if vertices are the same
            -- check whether s1.pts(k).pt.all occurs in s2.pts
            found := false;
            for l in s2.pts'range loop
              if s1.pts(k).pt.all = s2.pts(l).pt.all
               then found := true;
              end if;
              exit when found;
            end loop;
            if not found
             then return false;
            end if;
          end loop;
          return true;
    end if;
  end Equal;

  function Index ( s : Simplex; x : Vector ) return integer32 is
  begin
    for k in s.pts'range loop
      if s.pts(k).pt.all = x
       then return k;
      end if;
    end loop;
    return 0;
  end Index;

  function Neighbor ( s : Simplex; k : integer32 ) return Simplex is
  begin
    return s.pts(k).si;
  end Neighbor;

  function Neighbor ( s : Simplex; k : integer32; pos : Vector )
                    return Simplex is
  begin
    if pos(k-1)*pos(pos'last) > 0
     then return s.pts(k).si;
     else return Null_Simplex;
    end if;
  end Neighbor;

  function Position ( s : Simplex; x : Vector ) return Vector is

    m : matrix(x'first..x'last-1,x'range);
    pos : Vector(x'range);
    res : Vector(0..pos'last);
 
  begin 
   -- nbpos := nbpos + 1;
   -- put("# position computations : "); put(nbpos,1); new_line;
   -- transform point and simplex
    m := Diagonalize(s,x);
   -- solve the system
    pos := (pos'range => 0);
    Solve0(m,pos);
    res(pos'first..pos'last) := pos;
    res(0) := 0;
    for k in pos'range loop
      res(0) := res(0) + pos(k);
    end loop;
    res(0) := -res(0);
    return res;
  end Position;

  function Is_In ( s : Simplex; x : Vector ) return boolean is

    pos : constant Vector(0..x'last) := Position(s,x);

  begin
    return Is_In(pos);
  end Is_In;
    
  function Is_In ( pos : Vector ) return boolean is
  begin
    for k in pos'first..pos'last-1 loop
      if pos(k)*pos(pos'last) > 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_In;

  function Is_In_All ( s : Simplex; x : Vector ) return boolean is

    pos : constant Vector(0..x'last) := Position(s,x);

  begin
    return Is_In_All(s,x,pos);
  end Is_In_All;

  function Is_In_All ( s : Simplex; x : Vector ) return Simplex is

    pos : constant Vector(0..x'last) := Position(s,x);

  begin
    return Is_In_All(s,x,pos);
  end Is_In_All;

  function Is_In_All ( s : Simplex; x,pos : Vector ) return boolean is

    ins : boolean := true;    -- assumes that x belongs to s

  begin
    for k in pos'first..pos'last-1 loop
      if pos(k)*pos(pos'last) > 0
       then if s.pts(k+1).si /= Null_Simplex
             then return Is_In_All(s.pts(k+1).si,x);
             else ins := false;
            end if;
      end if;
    end loop;
    return ins;
  end Is_In_All;

  function Is_In_All ( s : Simplex; x,pos : Vector ) return Simplex is

    ins : boolean := true;    -- assumes that x belongs to s

  begin
    for k in pos'first..pos'last-1 loop
      if pos(k)*pos(pos'last) > 0
       then if s.pts(k+1).si /= Null_Simplex
             then return Is_In_All(s.pts(k+1).si,x);
             else ins := false;
            end if;
      end if;
    end loop;
    if ins
     then return s;
     else return Null_Simplex;
    end if;
  end Is_In_All;

  procedure Neighbors ( s : in out Simplex; x : in Vector ) is

    cont : boolean;

    procedure Neighbors ( s : in out Simplex; x : in Vector;
                          cont : out boolean ) is

      pos : constant Vector(0..x'last) := Position(s,x);
      continue : boolean := true;

    begin
      for k in pos'first..pos'last-1 loop
        if pos(k)*pos(pos'last) > 0 then
          if s.pts(k+1).si /= Null_Simplex
           then Neighbors(s.pts(k+1).si,x,continue);
           else Process_Neighbor(s,k+1,continue);
          end if;
        end if;
        exit when not continue;
      end loop;
      cont := continue;
    end Neighbors;
  
  begin
    Neighbors(s,x,cont);
  end Neighbors;

  function Volume ( s : Simplex ) return natural32 is

    m : constant matrix := Diagonalize(s);
    vol : integer32 := 1;

  begin
    for k in m'range(1) loop
      vol := vol*m(k,k);
    end loop;
    if vol >= 0
     then return natural32(vol);
     else return natural32(-vol);
    end if;
  end Volume;

-- DESTRUCTORS :

  procedure Destroy_Neighbor ( s : in out Simplex; k : in integer32 ) is
  begin
    s.pts(k).si := Null_Simplex;
  end Destroy_Neighbor;

  procedure Destroy_Neighbors ( s : in out Simplex ) is
  begin
    for k in s.pts'range loop
      Destroy_Neighbor(s,k);
    end loop;
  end Destroy_Neighbors;

  procedure Clear_Neighbor ( s : in out Simplex; k : in integer32 ) is
  begin
    Clear(s.pts(k).si);
  end Clear_Neighbor;

  procedure Clear_Neighbors ( s : in out Simplex ) is
  begin
    for k in s.pts'range loop
      Clear_Neighbor(s,k);
    end loop;
  end Clear_Neighbors;

  procedure Clear ( s : in out Simplex ) is

    procedure free is new unchecked_deallocation(Simplex_Rep,Simplex);

  begin
    Clear(s.tra);
    free(s);
  end Clear;

end Standard_Integer32_Simplices;
