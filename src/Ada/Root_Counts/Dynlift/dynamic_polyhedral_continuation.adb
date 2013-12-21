with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Laurentials;       use Standard_Complex_Laurentials;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Standard_Simpomial_Solvers;
with Transforming_Laurent_Systems;
with Standard_Integer32_Simplices;       use Standard_Integer32_Simplices;
with Standard_Dynamic32_Triangulations;  use Standard_Dynamic32_Triangulations;
with Cayley_Trick;                       use Cayley_Trick;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Unfolding_Subdivisions;             use Unfolding_Subdivisions;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;

package body Dynamic_Polyhedral_Continuation is

-- CAUTION : PATCH FOR FEWNOMIALS => NO SHIFT !!!

-- AUXILIAIRIES :

  procedure Flatten ( t : in out Term ) is

  -- DESCRIPTION :
  --   Flattens the Laurent term, i.e., the last exponent of t becomes zero.

  begin
    t.dg(t.dg'last) := 0;
  end Flatten;

  procedure Flatten ( p : in out Poly ) is

  -- DESCRIPTION :
  --   Flattens the Laurent polynomial,
  --   i.e., the last exponent of every monomial becomes zero.

    procedure Flatten_Term ( t : in out Term; cont : out boolean ) is
    begin
      Flatten(t); cont := true;
    end Flatten_Term;
    procedure Flatten_Terms is new Changing_Iterator(Flatten_Term);

  begin
    Flatten_Terms(p);
  end Flatten;

  procedure Flatten ( p : in out Laur_Sys ) is

  -- DESCRIPTION :
  --   Flattens the Laurent polynomial system,
  --   i.e., the last exponent of every monomial becomes zero.

  begin
    for i in p'range loop
      Flatten(p(i));
    end loop;
  end Flatten;

  function Non_Flattened_Points ( l : List ) return List is

  -- DESCRIPTION :
  --   Returns the list of points with last coordinate /= 0.
   
    tmp,res,res_last : List;
    pt : Link_to_Vector;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt(pt'last) /= 0
       then Append(res,res_last,pt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Non_Flattened_Points;

  function Non_Flattened_Points ( l : Array_of_Lists )
                                return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the points with last coordinate /= 0.

    res : Array_of_Lists(l'range);

  begin
    for i in l'range loop
      res(i) := Non_Flattened_Points(l(i));
    end loop;
    return res;
  end Non_Flattened_Points;

  function Non_Flattened_Points ( fs : Face_Structures )
                                return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the points with last coordinate /= 0.

    res : Array_of_Lists(fs'range);

  begin
    for i in fs'range loop
      res(i) := Non_Flattened_Points(fs(i).l);
    end loop;
    return res;
  end Non_Flattened_Points;

  function Non_Flattened_Points
               ( mix : Vector; mixsub : Mixed_Subdivision )
               return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the points with last coordinate /= 0.

    res,res_last : Array_of_Lists(mix'range);
    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : constant Mixed_Cell := Head_Of(tmp);
      begin
        for i in mic.pts'range loop
          Deep_Concat_Diff(res(i),res_last(i),Non_Flattened_Points(mic.pts(i)));
        end loop;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Non_Flattened_Points;

  function Convert ( s : Simplex ) return Mixed_Cell is

    res : Mixed_Cell;

  begin
    res.nor := new vector'(Normal(s));
    res.pts := new Array_of_Lists(1..1);
    res.pts(1) := Shallow_Create(Vertices(s));
    return res;
  end Convert; 

  function Convert ( normal : in vector; s : Simplex ) return Mixed_Cell is

    res : Mixed_Cell;

  begin
    res.nor := new vector'(normal);
    res.pts := new Array_of_Lists(1..1);
    res.pts(1) := Shallow_Create(Vertices(s));
    return res;
  end Convert;

  function Non_Flattened_Cells
                 ( flatnor : Link_to_Vector; mixsub : Mixed_Subdivision )
                 return Mixed_Subdivision is

  -- DESCRIPTION :
  --   Returns the cells which are not flattened.

    tmp,res,res_last : Mixed_Subdivision;
    mic : Mixed_Cell;

  begin
    tmp := mixsub;
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      if mic.nor.all /= flatnor.all
       then Append(res,res_last,mic);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Non_Flattened_Cells;

  function Convert ( t : Triangulation ) return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : constant Mixed_Cell := Convert(Head_Of(tmp));
      begin
        Append(res,res_last,mic);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Non_Flattened_Cells
                 ( flatnor : Link_to_Vector; t : Triangulation )
                 return Mixed_Subdivision is

    tmp : Triangulation;
    res,res_last : Mixed_Subdivision;
   
  begin
    tmp := t;
    while not Is_Null(tmp) loop
      declare
        s : constant Simplex := Head_Of(tmp);
        nor : constant vector := Normal(s);
      begin
        if nor /= flatnor.all then
          declare
            mic : constant Mixed_Cell := Convert(nor,s);
          begin
            Append(res,res_last,mic);
          end;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Non_Flattened_Cells;

-- UTILITIES :

  function Pointer_to_Last ( l : Solution_List ) return Solution_List is

  -- DESCRIPTION :
  --   Returns a pointer to the last element of l.

  begin
    if Is_Null(l)
     then return l;
     elsif Is_Null(Tail_Of(l))
         then return l;
         else return Pointer_to_Last(Tail_Of(l));
    end if;
  end Pointer_to_Last;

  function Initialize_Polyhedral_Homotopy
                ( n :  integer32; mix : Vector; fs : Face_Structures;
                  p : Laur_Sys ) return Laur_Sys is

  -- DESCRIPTION :
  --   Initializes the polyhedral homotopy w.r.t. to the already lifted
  --   points in the face structure.

    res : Laur_Sys(p'range);

  begin
    if not Is_Null(fs(fs'first).l) then
      declare
        lifted : Array_of_Lists(fs'range);
      begin
        for i in lifted'range loop
          lifted(i) := fs(i).l;
        end loop;
        res := Perform_Lifting(n,mix,lifted,p);
      end;
    end if;
    return res;
  end Initialize_Polyhedral_Homotopy;

  function Initialize_Polyhedral_Homotopy
             ( n : integer32; L : List; p : Laur_Sys ) return Laur_Sys is

  -- DESCRIPTION :
  --   Initializes the polyhedral homotopy w.r.t. to the already lifted
  --   points in the list.

    res : Laur_Sys(p'range);

  begin
    if not Is_Null(L) then
      declare
        lifted : Array_of_Lists(p'range);
        mix : constant Vector(1..1) := (1..1 => n);
      begin
        for i in lifted'range loop
          lifted(i) := L;
        end loop;
        res := Perform_Lifting(n,mix,lifted,p);
      end;
    end if;
    return res;
  end Initialize_Polyhedral_Homotopy;

  function Select_Terms ( p : Poly; L : List ) return Poly is

  -- DESCRIPTION :
  --   Given a tuple of lifted points, selected terms of p, whose exponents
  --   occur in l, will be returned.

    res : Poly := Null_Poly;
    tmp : List := L;
    point : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      point := Head_Of(tmp);
      declare
        d : degrees := new Vector(point'first..point'last-1);
        t : Term;
      begin
        d.all := point(point'first..point'last-1);
        t.cf := Coeff(p,d);
        t.dg := d;
        Add(res,t);
        Standard_Integer_Vectors.Clear
          (Standard_Integer_Vectors.Link_to_Vector(d));
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : Laur_Sys; mix : Vector; L : Array_of_Lists )
                        return Laur_Sys is

  -- DESCRIPTION :
  --   Given a tuple of lifted points, selected terms of p, whose exponents
  --   occur in l, will be returned.

    res : Laur_Sys(p'range);
    cnt : integer32 := p'first;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        res(cnt) := Select_Terms(p(cnt),L(i));
        cnt := cnt + 1;
      end loop;
    end loop;
    return res;
  end Select_Terms;

  procedure Update_Polyhedral_Homotopy
              ( p : in Laur_Sys; point : in Vector; i : in integer32;
                hom : in out Laur_Sys ) is

  -- DESCRIPTION :
  --   Given the lifted point of the ith support list,
  --   the ith polynomial of hom will be updated with a new term.

    d : degrees
      := new Standard_Integer_Vectors.Vector(point'first..point'last-1);
    t : Term;

  begin
    d.all := point(point'first..point'last-1);
    t.cf := Coeff(p(i),d);
    t.dg := new Standard_Integer_Vectors.Vector'(point);
    Add(hom(i),t);
    Standard_Integer_Vectors.Clear(Standard_Integer_Vectors.Link_to_Vector(d));
    Clear(t);
  end Update_Polyhedral_homotopy;

  procedure Update_Polyhedral_Homotopy
              ( p : in Laur_Sys; L : in List; i : in integer32;
                hom : in out Laur_Sys ) is

  -- DESCRIPTION :
  --   Given a list of lifted points of the ith support list,
  --   the ith polynomial of hom will be updated with new terms.

    tmp : List := L;

  begin
    while not Is_Null(tmp) loop
      Update_Polyhedral_Homotopy(p,Head_Of(tmp).all,i,hom);
      tmp := Tail_Of(tmp);
    end loop;
  end Update_Polyhedral_Homotopy;

  procedure Update_Polyhedral_Homotopy
               ( p : in Laur_Sys; lifted : in Array_of_Lists;
                 mix : in Vector; hom : in out Laur_Sys ) is

  -- DESCRIPTION :
  --   Given a lists of lifted points, the polyhedral homotopy hom
  --   will be updated with new terms.

    cnt : integer32 := hom'first;

  begin
    for i in mix'range loop
      for j in 1..mix(i) loop
        Update_Polyhedral_Homotopy(p,lifted(i),cnt,hom);
        cnt := cnt + 1;
      end loop;
    end loop;
  end Update_Polyhedral_Homotopy;

  procedure Update_Polyhedral_Homotopy
               ( p : in Laur_Sys; lifted : in List; hom : in out Laur_Sys ) is

  -- DESCRIPTION :
  --   Given a list of lifted points, the unmixed polyhedral homotopy hom
  --   will be updated with new terms.

  begin
    for i in hom'range loop
      Update_Polyhedral_Homotopy(p,lifted,i,hom);
    end loop;
  end Update_Polyhedral_Homotopy;

  procedure Solve_by_Unfolding
                ( file : in file_type; p,hom : in Laur_Sys; n : in integer32;
                  mix : in Vector; nor : in Link_to_Vector;
                  mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   All cells in the given mixed subdivision have the same normal.
  --   The solutions which correspond to these cells will be computed
  --   by means of the given polyhedral homotopy.

  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   p          randomized system, without lifting;
  --   hom        the polyhedral homotopy;
  --   n          dimension before lifting;
  --   mix        type of mixture;
  --   nor        the same normal to all cells in the subdivision;
  --   mixsub     the mixed subdivision all with the same normal nor.

  -- ON RETURN :
  --   sols       the solutions of p, w.r.t. the cells in mixsub.

    work : Mixed_Subdivision;
    mv : natural32;
    homsub : Laur_Sys(p'range);
    first : boolean := true;
    flatnor : Link_to_Vector;
    sols_last : Solution_List;

    procedure Process ( mic : in Mixed_Cell; newpts : in Array_of_Lists ) is

      subsys : Laur_Sys(p'range) := Select_Terms(p,mix,mic.pts.all);
      subsols : Solution_List;
      tol_zero : constant double_float := 1.0E-12;
      fail,zero_y : boolean;

    begin
      Transforming_Laurent_Systems.Shift(subsys);  -- patch for Simpomials !!
      Standard_Simpomial_Solvers.Solve(subsys,tol_zero,subsols,fail,zero_y);
      put(file,"Number of solutions of subsystem : ");
      put(file,Length_Of(subsols),1); new_line(file);
      if first then
        homsub := Perform_Lifting(n,mix,mic.pts.all,p);
        sols := subsols;
        first := false;
      else
        Update_Polyhedral_Homotopy(p,newpts,mix,homsub);
        Mixed_Continuation(file,homsub,mic.nor.all,subsols);
        Set_Continuation_Parameter(sols,Create(0.0));
        Mixed_Continuation(file,homsub,flatnor.all,sols);
        Flatten(homsub);
        sols_last := Pointer_to_Last(sols);
        Concat(sols,sols_last,subsols);
        Shallow_Clear(subsols);
      end if;
      Clear(subsys);
    end Process;
    procedure Solve_Subsystem is new Unfolding(Process);

  begin
    new_line(file);
    put_line(file,"****  SOLVING BY UNFOLDING  ****");
    new_line(file);
    flatnor := new vector(1..n+1);
    flatnor(1..n) := (1..n => 0);
    flatnor(n+1) := 1;
    Copy(mixsub,work);
    put(file,natural32(n),mix,work,mv);
    Solve_Subsystem(work);
    --Deep_Clear(work);
    Clear(flatnor); Clear(homsub);
    Set_Continuation_Parameter(sols,Create(0.0));
    Mixed_Continuation(file,hom,nor.all,sols);
  end Solve_by_Unfolding;

  procedure Solve_with_Unfolding
              ( file : in file_type; p,hom : in Laur_Sys; n : in integer32; 
                mix : in Vector; mixsub : in Mixed_Subdivision;
                sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Computes the new solutions corresponding to the cells in the
  --   mixed subdivision.

    sols_last : Solution_List;
    newnor : List := Different_Normals(mixsub);
    tmp : List := newnor;
    normal : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      normal := Head_Of(tmp);
      declare
        newcells : Mixed_Subdivision := Extract(normal.all,mixsub);
        bigcell : Mixed_Subdivision;
        newsols : Solution_List;
      begin
       put_line(file,"THE LIST OF NEW CELLS");
       put(file,natural32(n),mix,newcells);
        if Length_Of(newcells) = 1
         then Mixed_Solve(file,hom,mix,newcells,newsols);
      --   else Solve_by_Unfolding(file,p,hom,n,mix,normal,newcells,newsols);
         else bigcell := Merge_Same_Normal(newcells);
              put_line(file,"THE BIG CELL");
              put(file,natural32(n),mix,bigcell);
              Mixed_Solve(file,hom,mix,bigcell,newsols);
             -- Deep_Clear(bigcell);
              Shallow_Clear(bigcell);
        end if;
        sols_last := Pointer_to_Last(sols);
        Concat(sols,sols_last,newsols);
        Shallow_Clear(newsols);
        Shallow_Clear(newcells);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Deep_Clear(newnor);
  end Solve_with_Unfolding;

-- DYNAMIC LIFTING FOR UNMIXED SYSTEMS :

  procedure Dynamic_Unmixed_Solve
              ( file : in file_type; n : in integer32;
                L : in List; order,inter : in boolean; maxli : in integer32;
                lifted,lifted_last : in out List; t : in out Triangulation;
                q : in Poly_Sys; qsols : in out Solution_List ) is

    p : Laur_Sys(q'range) := Polynomial_to_Laurent_System(q);
    qt : Laur_Sys(q'range);         -- the polyhedral homotopy
    firstflat : boolean := true;    -- first time flattening
    flatnor : Link_to_Vector;       -- the flat normal
    mix : Vector(1..1) := (1..1 => n);

    procedure Solve_Before_Flattening
                 ( t : in Triangulation; lifted : in List ) is

    -- DESCRIPTION :
    --   Computes the new solutions, right before the subdivision
    --   is flattened.

      newpts : List;
      newcells : Mixed_Subdivision;

    begin
      if firstflat
       then newpts := lifted;
            newcells := Convert(t);
       else newcells := Non_Flattened_Cells(flatnor,t);
            newpts := Non_Flattened_Points(lifted);
      end if;
      Update_Polyhedral_Homotopy(p,newpts,qt);
      if not firstflat
       then new_line(file);
            put_line(file,"***  EXTENDING THE SOLUTIONS  ***");
            new_line(file);
            Set_Continuation_Parameter(qsols,Create(0.0));
            Mixed_Continuation(file,qt,flatnor.all,qsols);
      end if;
      Solve_with_Unfolding(file,p,qt,n,mix,newcells,qsols);
      if not firstflat
       then Shallow_Clear(newpts); Shallow_Clear(newcells);
       else firstflat := false;
      end if;
      Flatten(qt);
    end Solve_Before_Flattening;

    procedure S_Dynamic_Lifting is
      new Standard_Dynamic32_Triangulations.Dynamic_Lifting_with_Flat
                 ( Before_Flattening => Solve_Before_Flattening );

  begin
    qt := Initialize_Polyhedral_Homotopy(n,lifted,p);
    flatnor := new vector(1..n+1);
    flatnor(1..n) := (1..n => 0);
    flatnor(n+1) := 1;
    S_Dynamic_Lifting(l,order,inter,maxli,lifted,lifted_last,t);
    Solve_Before_Flattening(t,lifted);
    Clear(flatnor); Clear(mix);
    Clear(qt); Clear(p);
  end Dynamic_Unmixed_Solve;

-- DYNAMIC LIFTING FOR THE CAYLEY TRICK :

  procedure Dynamic_Cayley_Solve
              ( file : in file_type; n : in integer32; mix : in Vector;
                supports : in Array_of_Lists; order,inter : in boolean;
                maxli : in integer32; lifted : in out Array_of_Lists;
                mixsub : in out Mixed_Subdivision; numtri : out natural32;
                q : in Poly_Sys; qsols : in out Solution_List ) is

    p : constant Laur_Sys(q'range) := Polynomial_to_Laurent_System(q);
    qt : Laur_Sys(q'range);         -- the polyhedral homotopy
    firstflat : boolean := true;    -- first time flattening
    flatnor : Link_to_Vector;       -- the flat normal

    procedure Solve_Before_Flattening
                ( mixsub : in out Mixed_Subdivision;
                  lifted : in Array_of_Lists ) is

    -- DESCRIPTION :
    --   Computes the new solutions, right before the subdivision
    --   is flattened.

      newpts : Array_of_Lists(lifted'range);
      newcells : Mixed_Subdivision;

    begin
      if firstflat
       then newpts := lifted;
            newcells := mixsub;
       else newcells := Non_Flattened_Cells(flatnor,mixsub);
            newpts := Non_Flattened_Points(lifted);
      end if;
      Update_Polyhedral_Homotopy(p,newpts,mix,qt);
      if not Is_Null(qsols)
       then new_line(file);
            put_line(file,"***  EXTENDING THE SOLUTIONS  ***");
            new_line(file);
            Set_Continuation_Parameter(qsols,Create(0.0));
            Mixed_Continuation(file,qt,flatnor.all,qsols);
      end if;
      if not Is_Null(newcells)
       then
         Solve_with_Unfolding(file,p,qt,n,mix,newcells,qsols);
         if not firstflat
          then Shallow_Clear(newpts); Shallow_Clear(newcells);
          else firstflat := false;
         end if;
      end if;
      Flatten(qt);
    end Solve_Before_Flattening;

    procedure S_Dynamic_Lifting is
      new Cayley_Trick.Dynamic_Cayley_with_Flat
                ( Before_Flattening => Solve_Before_Flattening );

  begin
    flatnor := new vector(1..n+1);
    flatnor(1..n) := (1..n => 0);
    flatnor(n+1) := 1;
    S_Dynamic_Lifting(n,mix,supports,order,inter,maxli,lifted,mixsub,numtri);
    Solve_Before_Flattening(mixsub,lifted);
    Clear(flatnor);
  end Dynamic_Cayley_Solve;

-- DYNAMIC LIFTING FOR SEMI-MIXED SYSTEMS :

  procedure Dynamic_Mixed_Solve
                ( file : in file_type; n : in integer32;
                  mix : in Vector; supports : in Array_of_Lists;
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                  q : in Poly_Sys; qsols : in out Solution_List ) is

    p : Laur_Sys(q'range) := Polynomial_to_Laurent_System(q);
    qt : Laur_Sys(q'range);         -- the polyhedral homotopy
    firstflat : boolean := true;    -- first time flattening 
    flatnor : Link_to_Vector;       -- the flat normal
    qsols_last : Solution_List;

    procedure Solve_Before_Flattening
                ( mixsub : in out Mixed_Subdivision;
                  fs : in Face_Structures ) is

    -- DESCRIPTION :
    --   Computes the new solutions, right before the subdivision 
    --   is flattened.

      newpts : Array_of_Lists(fs'range);
      newcells : Mixed_Subdivision;
      newsols : Solution_List;

    begin
      if firstflat
       then for i in fs'range loop
              newpts(i) := fs(i).l;
            end loop;
            newcells := mixsub;
       else newcells := Non_Flattened_Cells(flatnor,mixsub);
            newpts := Non_Flattened_Points(fs);
      end if;
      Update_Polyhedral_Homotopy(p,newpts,mix,qt);
      if not firstflat
       then new_line(file);
            put_line(file,"***  EXTENDING THE SOLUTIONS  ***");
            new_line(file);
            Set_Continuation_Parameter(qsols,Create(0.0));
            Mixed_Continuation(file,qt,flatnor.all,qsols);
      end if;
      Mixed_Solve(file,qt,mix,newcells,newsols);
      qsols_last := Pointer_to_Last(qsols);
      Concat(qsols,qsols_last,newsols);
      Shallow_Clear(newsols);
      if not firstflat
       then Shallow_Clear(newpts); Shallow_Clear(newcells);
       else firstflat := false;
      end if;
      Flatten(qt);
    end Solve_Before_Flattening;

    procedure S_Dynamic_Lifting is
      new Dynamic_Mixed_Subdivisions.Dynamic_Lifting_with_Flat
                ( Before_Flattening => Solve_Before_Flattening );

  begin
    qt := Initialize_Polyhedral_Homotopy(n,mix,fs,p);
    flatnor := new vector(1..n+1);
    flatnor(1..n) := (1..n => 0);
    flatnor(n+1) := 1;
    S_Dynamic_Lifting(n,mix,supports,order,inter,conmv,maxli,mixsub,fs,
                      nbsucc,nbfail);
    Solve_Before_Flattening(mixsub,fs);
    Clear(flatnor);
    Clear(qt); Clear(p);
  end Dynamic_Mixed_Solve;

end Dynamic_Polyhedral_Continuation;
