with unchecked_deallocation;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Generic_Lists;
with Standard_Natural_Vectors;           use Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Set_Structure;                      use Set_Structure;
with Permutations,Permute_Operations;    use Permutations,Permute_Operations;
with Templates;                          use Templates;

package body Symmetric_Set_Structure is

-- DATASTRUCTURES :

  type set is array ( natural32 range <> ) of boolean;
  type boolean_array is array ( natural32 range <> ) of boolean;
  type link_to_boolean_array is access boolean_array;
  procedure free is 
    new unchecked_deallocation(boolean_array,link_to_boolean_array);

  type boolean_matrix is array ( natural32 range <> ) of link_to_boolean_array;
  type link_to_boolean_matrix is access boolean_matrix;
  procedure free is 
    new unchecked_deallocation(boolean_matrix,link_to_boolean_matrix);
  type set_coord is record
    k,l : natural32;
  end record;
  type Dependency_Structure is array ( natural32 range <> ) of set_coord;
  type Link_to_Dependency_Structure is access Dependency_Structure;
  procedure free is 
    new unchecked_deallocation
          (Dependency_Structure,Link_to_Dependency_Structure);

  package Lists_of_Dependency_Structures
    is new Generic_Lists (Link_to_Dependency_Structure);
  type Covering is new Lists_of_Dependency_Structures.List;

-- INTERNAL DATA : 

  cov : Covering;  -- covering of the set structure
  lbm : link_to_boolean_matrix;
    -- auxiliary data structure for bookeeping during the construction
    -- of the covering,
    -- to remember which sets have already been treated.

-- AUXILIARY ROUTINES FOR CONSTRUCTING THE COVERING :

  function Give_Set ( n,i,j : natural32 ) return set is
   
  -- DESCRIPTION :
  --   Returns the (i,j)-th set out of the set structure.

    s : set(1..n);

  begin
    for k in 1..n loop
      s(k) := Is_In(i,j,k);
    end loop;
    return s;
  end Give_Set;

  function Equal ( s1,s2 : set ) return boolean is

  -- DESCRIPTION :
  --   Returns true if both sets are equal.

  begin
    for i in s1'range loop
      if s1(i) /= s2(i)
       then return false;
      end if;
    end loop;
    return true;
  end Equal;

  function Find ( i,n : natural32; s : set ) return natural32 is

  -- DESCRIPTION :
  --   Returns the first occurence of the set s in the i-th row
  --   of the set structure;
  --   returns zero if the set does not occur in the i-th row.

  begin
    for j in 1..Number_Of_Sets(i) loop
      if not lbm(i)(j) and then Equal(s,Give_Set(n,i,j))
       then return j;
      end if;
    end loop;
    return 0;
  end Find;

  function Apply ( p : Permutation; s : set ) return set is

  -- DESCRIPTION :
  --   Returns the result after application of p on the set s.

    r : set(s'range);

  begin
    for i in p'range loop
      r(natural32(i)) := s(natural32(p(i)));
    end loop;
    return r;
  end Apply;

  procedure Init_Covering ( n : in natural32 ) is

  -- DESCRIPTION :
  --   Initialization of lbm.

  begin
    lbm := new boolean_matrix(1..n);
    for i in 1..n loop
      lbm(i) := new boolean_array'(1..Number_of_Sets(i) => false);
    end loop;
  end Init_Covering;

  procedure Update ( dps : Dependency_Structure ) is

  -- DESCRIPTION :
  --   All pairs in dps are marked in lbm.

  begin
    for i in dps'range loop
      lbm(dps(i).k)(dps(i).l) := true;
    end loop;
  end Update;

  procedure Search ( n : in natural32; i,j : out natural32;
                     empty : out boolean ) is

  -- DESCRIPTION :
  --   Searches in lbm the first (i,j)-th free set;
  --   returns empty if all sets have already been used.

  begin
    for k in 1..n loop
      for l in lbm(k)'range loop
        if not lbm(k)(l)
         then i := k; j := l; empty := false; return;
        end if;
      end loop;
    end loop;
    empty := true;
  end Search;

-- CONSTRUCTOR FOR DEPENDENCY STRUCTURE AND COVERING :

  procedure Construct_Dependency_Structure
                ( n,m : in natural32; v,w : in List_Of_Permutations;
                  i,j : in natural32; dps : in out Dependency_Structure;
		  fail : out boolean ) is

  -- DESCRIPTION :
  --   A dependency structure will be constructed.

  -- ON ENTRY :
  --   n         the dimension;
  --   m         number of elements in dps,v and w;
  --   v,w       matrix representations;
  --   i,j       coordinates of a set in the dependency structure.

  -- ON RETURN :
  --   dps       the dependency structure;
  --   fail      is true if the set structure is not symmetric.

    s : constant set(1..n) := Give_Set(n,i,j);
    lv,lw : List_Of_Permutations;
    pv,pw : Permutation(1..integer32(n));
    ps : set(1..n);
    res : natural32;

  begin
    lv := v;  lw := w;
    for x in 1..m loop
      pw := Permutation(Head_Of(lw).all);
      dps(x).k := natural32(pw(integer32(i)));
      pv := Permutation(Head_Of(lv).all);
      ps := Apply(pv,s);
      res := Find(dps(x).k,n,ps);
      exit when (res = 0);
      dps(x).l := res;
      lv := Tail_Of(lv);
      lw := Tail_Of(lw);
    end loop;
    fail := (res = 0);
  end Construct_Dependency_Structure;

  procedure Construct_Covering
                ( n,m : in natural32; v,w : in List_Of_Permutations;
                  fail : out boolean ) is

  -- DESCRIPTION :
  --   A covering of the set structure will be constructed.

  -- EFFECT :
  --   Initially, all entries in lbm are false;
  --   at the end, all entries in lbm are true (if not fail).

    dps : Dependency_Structure(1..m);
    ldps : Link_to_Dependency_Structure;
    empty,fl : boolean;
    i,j : natural32;

  begin
    Init_Covering(n);
    Search(n,i,j,empty);
    while not empty loop
      Construct_Dependency_Structure(n,m,v,w,i,j,dps,fl);
      exit when fl;
      Update(dps);
      ldps := new Dependency_Structure(1..m);
      ldps.all := dps;
      Construct(ldps,cov);
      Search(n,i,j,empty);
    end loop;
    fail := fl;
  end Construct_Covering;

-- OUTPUT PROCEDURES FOR COVERING :

  procedure Write_Set ( n,i,j : natural32 ) is

  -- DESCRIPTION :
  --   Writes the (i,j)-th set on the standard output.

  begin
    put('{');
    for k in 1..n loop
      if Is_In(i,j,k)
       then put(' '); put('x'); put(k,1);
      end if;
    end loop;
    put(" }");
  end Write_Set;

  procedure Write_Coord ( k,l : in natural32 ) is
  begin
    put('['); put(k,1); put(' '); put(l,1); put(']');
  end Write_Coord;

  procedure Write_Covering is

    tmp : Covering := cov;
    ldps : Link_to_Dependency_Structure;

  begin
    put_line("The covering :");
    while not Is_Null(tmp) loop
      ldps := Head_Of(tmp);
      declare
	nb : natural32 := 0;
      begin
        for i in ldps'range loop
  	  Write_Coord(ldps(i).k,ldps(i).l);
	  nb := nb+1;
	  if nb > 7
	   then new_line; nb := 0;
          end if;
        end loop;
        new_line;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Covering;

  procedure Write_Coord ( file : in file_type; k,l : in natural32 ) is
  begin
    put(file,'['); put(file,k,1); put(file,' '); put(file,l,1); put(file,']');
  end Write_Coord;

  procedure Write_Covering ( file : in file_type ) is

    tmp : Covering := cov;
    ldps : Link_to_Dependency_Structure;

  begin
    put_line(file,"The covering :");
    while not Is_Null(tmp) loop
      ldps := Head_Of(tmp);
      declare
        nb : natural := 0;
      begin
        for i in ldps'range loop
          Write_Coord(file,ldps(i).k,ldps(i).l);
          nb := nb+1;
          if nb > 7
           then new_line(file); nb := 0;
          end if;
        end loop;
        new_line(file);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Covering;

-- CONSTRUCTION OF TEMPLATES :

  procedure Init_Template ( n : in natural32 ) is

  -- DESCRIPTION :
  --   Initialization of the template.

    h : constant Standard_Natural_Vectors.Vector(0..integer32(n))
      := (0..integer32(n) => 0);

  begin
    Templates.Create(n);
    for i in 1..n loop
      for j in 1..Number_Of_Sets(i) loop
        Templates.Add_Hyperplane(i,h);
      end loop;
    end loop;
  end Init_Template;

  procedure First_Equivariant_Template
                ( n : in natural32; cnt : in out natural32 ) is

  -- DESCRIPTION :
  --   Constructs the first equation of the template, for an equivariant
  --   linear product system system

  -- ON ENTRY :
  --   n          the dimension;
  --   cnt        counts the number of free coefficients.

    h : Standard_Natural_Vectors.Vector(0..integer32(n));

  begin
    for j in 1..Templates.Number_of_Hyperplanes(1) loop
      Templates.Get_Hyperplane(1,j,h);
      cnt := cnt + 1; h(0) := cnt;
      for k in 1..n loop
        if Set_Structure.Is_In(1,j,k) then
          if cnt = h(0)
           then cnt := cnt + 1;
          end if;
          h(integer32(k)) := cnt;
        end if;
      end loop;
      Templates.Change_Hyperplane(1,j,h);
    end loop;
  end First_Equivariant_Template;

  function Action ( i,n : natural32; g : List_of_Permutations )
                  return Permutation is

  -- DESCRIPTION :
  --   Returns the group action from the list g that permutes the first
  --   array of sets into the ith one.

    p : Permutation(1..integer32(n));
    first,second : Standard_Natural_Vectors.Vector(1..integer32(n));
    tmp : List_of_Permutations := g;

  begin
    for k in 1..integer32(n) loop
      if Set_Structure.Is_In(1,1,natural32(k))
       then first(k) := 1;
       else first(k) := 0;
      end if;
      if Set_Structure.Is_In(i,1,natural32(k))
       then second(k) := 1;
       else second(k) := 0;
      end if;
    end loop;
    while not Is_Null(tmp) loop
      p := Permutation(Head_Of(tmp).all);
      if second = p*first
       then return p;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    p := (p'range => 0);
    return p;
  end Action;

  procedure Propagate_Equivariant_Template
                   ( n : in integer32; g : in List_of_Permutations;
                     fail : out boolean ) is

  -- DESCRIPTION :
  --   Given a template whose first equation is already constructed,
  --   the rest of the template will be constructed, with the aid of the
  --   list of generating permutations.

    h : Standard_Natural_Vectors.Vector(0..n);
    p : Permutation(1..n);

  begin
    for i in 2..natural32(n) loop
      p := Action(i,natural32(n),g);
      if p = (p'range => 0)
       then fail := true; return;
      end if;
      for j in 1..Templates.Number_of_Hyperplanes(i) loop
        Templates.Get_Hyperplane(1,j,h);
        h(1..n) := p*h(1..n);
        Templates.Change_Hyperplane(i,j,h);
      end loop;
    end loop;
    fail := false;
  end Propagate_Equivariant_Template;

  procedure Construct_Part_of_Template
                ( n,m : in natural32; v : in List_Of_Permutations;
                  dps : in Dependency_Structure; invpv1 : in Permutation;
                  cnt : in out natural32 ) is

  -- DESCRIPTION :
  --   This procedure constructs the coefficients of the hyperplanes
  --   associated with the sets in the dependency structure dps.
  --   cnt counts the number of free coefficients.

    lv : List_Of_Permutations;
    pv : Permutation(1..integer32(n));
    h : Standard_Natural_Vectors.Vector(0..integer32(n));
    indi : integer32;

  begin
   -- GENERATE CONSTANT COEFFICIENT :
    cnt := cnt+1;
    for j in 1..m loop
      Templates.Get_Hyperplane(dps(j).k,dps(j).l,h);
      h(0) := cnt;
      Templates.Change_Hyperplane(dps(j).k,dps(j).l,h);
    end loop;
   -- GENERATE THE OTHER COEFFICIENTS :
    for i in 1..integer32(n) loop
     -- GENERATE :
      if Is_In(dps(1).k,dps(1).l,natural32(i)) then
        Templates.Get_Hyperplane(dps(1).k,dps(1).l,h);
        if h(i) = 0 then
          cnt := cnt + 1;
          -- PROPAGATE :
          --put("PROPAGATING "); put(i,1);
          --put_line("-th coefficient :");
          lv := v;
          for j in 1..m loop
            pv := Permutation(Head_Of(lv).all);
            indi := 0;
            for l in 1..integer32(n) loop
              if pv(l) = invpv1(i)
               then indi := l; exit;
              end if;
            end loop;
	    --Write_Coord(dps(j).k,dps(j).l); put(" : ");
	    --Write_Set(n,dps(j).k,dps(j).l);
            --put(" indi : "); put(indi,1); new_line;
            Templates.Get_Hyperplane(dps(j).k,dps(j).l,h);
            h(indi) := cnt;
            Templates.Change_Hyperplane(dps(j).k,dps(j).l,h);
            lv := Tail_Of(lv);
          end loop;
          --put_line("RANDOM PRODUCT SYSTEM AFTER PROPAGATION :");
          --Write_RPS(n,2,4,3);
          --for l in 1..75 loop put("+"); end loop; new_line;
        end if;
      end if;
    end loop;
  end Construct_Part_of_Template;

  procedure Construct_Template
		( n,m : in natural32; v : in List_Of_Permutations;
		  nbfree : out natural32 ) is

  -- DESCRIPTION :
  --   Given a covering of the set structure,
  --   the data of the package Random_Product_System will be filled.

  -- ON ENTRY :
  --   n          the dimension of the vectors
  --   m          the number of entries in v
  --   v          matrix representations of the group

  -- ON RETURN :
  --   nbfree     the number of free coefficients

    tmp : Covering := cov;
    ldps : Link_to_Dependency_Structure;
    invpv1 : Permutation(1..integer32(n));
    cnt : natural32;

  begin
    Init_Template(n);
    cnt := 0;
    -- CONSTRUCT THE BASE SET OF dps :
    invpv1 := inv(Permutation(Head_Of(v).all));
      -- then for each pv in v: permutation of the base set
      -- is defined as pv*invpv1.
     --put("invpv1 : "); Put(invpv1); new_line;
    while not Is_Null(tmp) loop
      ldps := Head_Of(tmp);
      Construct_Part_of_Template(n,m,v,ldps.all,invpv1,cnt);
      tmp := Tail_Of(tmp);
    end loop;
    nbfree := cnt;
  end Construct_Template;

  procedure Construct_Equivariant_Template
                 ( n : in natural32; g : in List_of_Permutations;
                   cntfree : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Constructs a template for an equivariant system.  The list g contains
  --   the generating elements of the group.  The variable cntfree counts the
  --   number of free coefficients.

  begin
    Init_Template(n);
    First_Equivariant_Template(n,cntfree);
    Propagate_Equivariant_Template(integer32(n),g,fail);
  end Construct_Equivariant_Template;

  procedure Write_Templates ( n : in natural32 ) is
  begin
    Write_Templates(Standard_Output,n);
  end Write_Templates;

  procedure Write_Templates ( file : in file_type; n : in natural32 ) is

    h : Standard_Natural_Vectors.Vector(0..integer32(n));

  begin
    put_line(file,"The templates :");
    for i in 1..n loop
      for j in 1..Number_of_Hyperplanes(i) loop
        put(file,"("); put(file,i,1); put(file,","); put(file,j,1);
        put(file,") : "); Get_Hyperplane(i,j,h); put(file,h); new_line(file);
      end loop;
    end loop;
  end Write_Templates;

-- CONSTRUCTION OF START SYSTEMS :

  procedure Equivariant_Start_System
                  ( n : in natural32; g : in List_of_Permutations;
                    fail : out boolean ) is

    nbfree : natural32 := 0;
    fl : boolean := false;

  begin
    Construct_Equivariant_Template(n,g,nbfree,fl);
    if not fl
     then Templates.Polynomial_System(n,nbfree);
    end if;
    fail := fl;
  end Equivariant_Start_System;

  procedure Symmetric_Start_System 
               ( n,bb : in natural32; lp : in List;
                 v,w : in List_Of_Permutations;
                 notsymmetric,degenerate : out boolean ) is

    m : constant natural32 := Number(v);
    fl : boolean;
    nbfree : natural32;

  begin
    Construct_Covering(n,m,v,w,fl);
   -- Write_Covering;
    for i in lbm'range loop
      free(lbm(i));
    end loop;
    free(lbm);
    if fl then
      notsymmetric := true;
      -- put_line("The set structure is not (G,V,W)-symmetric.");
    else
      notsymmetric := false;
      -- put_line("The set structure is (G,V,W)-symmetric.");
      -- Templates.Create(n);
      Construct_Template(n,m,v,nbfree);
      -- Write_Templates(n);
      -- vb := Templates.Verify(n,lp);
      -- put("The bound of Templates.Verify : "); put(vb,1); new_line;
      -- if bb /= vb
      --  then degenerate := true;
      --       put_line("The set structure is degenerate.");
      --  else 
      degenerate := false;
      --       put_line("The set structure is not degenerate.");
      Templates.Polynomial_System(n,nbfree);
      -- end if;
    end if;
  end Symmetric_Start_System;

-- DESTRUCTOR :

  procedure Clear is

    use Lists_of_Dependency_Structures;
    tmp : Covering := cov;
    elem : Link_to_Dependency_Structure;

  begin
    while not Is_Null(tmp) loop
      elem := Head_Of(tmp);
      free(elem);
      tmp := Tail_Of(tmp);
    end loop;
    Clear(cov);
    Templates.Clear;
  end Clear;

end Symmetric_Set_Structure;
