with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Initial_Mixed_Cell;
with Standard_Integer32_Vertices;        use Standard_Integer32_Vertices;
with Standard_Integer32_Simplices;       use Standard_Integer32_Simplices;
with Global_Dynamic32_Triangulation;     use Global_Dynamic32_Triangulation;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Dynamic32_Lifting_Functions;        use Dynamic32_Lifting_Functions;
with Flatten_Mixed_Subdivisions;         use Flatten_Mixed_Subdivisions;
with Enumerate_Faces_of_Polytope;        use Enumerate_Faces_of_Polytope;
with Common_Faces_of_Polytope;           use Common_Faces_of_Polytope;
with Integer_Pruning_Methods;            use Integer_Pruning_Methods;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Contributions_to_Mixed_Volume;      use Contributions_to_Mixed_Volume;

package body Dynamic_Mixed_Subdivisions is

-- UTILITIES :

  function Is_Empty ( points : Array_of_Lists ) return boolean is
  begin
    for i in points'range loop
      if not Is_Null(points(i))
       then return false;
      end if;
    end loop;
    return true;
  end Is_Empty;

  function First_Non_Empty ( points : Array_of_Lists ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of the first non empty set in the points.

    res : integer32 := points'first - 1;

  begin
    for i in points'range loop
      if not Is_Null(points(i))
       then res := i;
      end if;
      exit when res >= points'first;
    end loop;
    return res;
  end First_Non_Empty;

  function Is_In ( fs : Face_Structure; point : vector ) return boolean is
  begin
    if Is_Null(fs.t)
     then return Is_In(fs.l,point);
     else return Is_In(fs.t,point);
    end if;
  end Is_In;

--  procedure put ( n : in integer32; fs : in Face_Structure ) is
--  begin
--    put_line("The list of lifted points : "); put(fs.l);
--    put_line("The list of k-faces : "); put(fs.f);
--    put_line("The triangulation : "); put(n,fs.t);
--  end put;

--  procedure put ( n : in integer32; fs : in Face_Structures ) is
--  begin
--    for i in fs'range loop
--      put("face structure for component "); put(i,1);
--      put_line(" :");
--      put(n,fs(i));
--    end loop;
--  end put;

-- FLATTENING :

  procedure Flatten ( points : in out VecVec ) is

    pt : Link_to_Vector;

  begin
    for i in points'range loop
      pt := points(i);
      if pt(pt'last) /= 0
       then pt(pt'last) := 0;
      end if;
    end loop;
  end Flatten;

  procedure Flatten ( f : in out Face ) is
  begin
    Flatten(f.all);
  end Flatten;

  procedure Flatten ( fs : in out Faces ) is

    tmp : Faces := fs;
    f : Face;

  begin
    while not Is_Null(tmp) loop
      f := Head_Of(tmp);
      Flatten(f);
      Set_Head(tmp,f);
      tmp := Tail_Of(tmp);
    end loop;
  end Flatten;

  procedure Flatten ( fs : in out Face_Structure ) is
  begin
    Flatten(fs.l); Flatten(fs.t); Flatten(fs.f);
  end Flatten;

  procedure Flatten ( fs : in out Face_Structures ) is
  begin
    for i in fs'range loop
      Flatten(fs(i));
    end loop;
  end Flatten;

-- ZERO CONTRIBUTION :

  function Extract_Facet ( n : integer32; s : Simplex ) return Face is

    ver : constant VecVec := Standard_Integer32_Simplices.Vertices(s);
    res : constant Face := new VecVec(ver'range);

  begin
    for i in res'range loop
      res(i) := new Standard_Integer_Vectors.Vector'(ver(i)(1..n));
    end loop;
    return res;
  end Extract_Facet;

  function Extract_Facets ( n : integer32; t : Triangulation; x : Vector )
                          return Faces is

  -- DESCRIPTION :
  --   Returns the list of facets in t that all contain x.
  --   The facets are given in their original coordinates.

  -- REQUIRED : x is lifted, i.e., x'length = n+1.

    tmp : Triangulation := t;
    res,res_last : Faces;

  begin
    while not Is_Null(tmp) loop
      declare
        s : constant Simplex := Head_Of(tmp);
      begin
        if Is_Vertex(s,x)
         then Append(res,res_last,Extract_Facet(n,s));
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Extract_Facets;

  function Zero_Contribution ( n : integer32; fs : Face_Structures;
                               x : Vector; i : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the point x does not contribute to the mixed volume
  --   when considered to the ith component of the already lifted points.

  -- REQUIRED : x'length = n+1, n = dimension before lifting.

    res : boolean := false;
    supp : Array_of_Lists(fs'range);

  begin
    for j in supp'range loop
      supp(j) := Reduce(fs(j).l,n+1);
    end loop;
    if Length_Of(fs(i).t) > 1 then
      declare
        facets : constant Faces := Extract_Facets(n,fs(i).t,x);
      begin
        res := Simple_Zero_Contribution(supp,facets,x(1..n),i);
      end;
    else
      res := Simple_Zero_Contribution(supp,x(1..n),i);
    end if;
    Deep_Clear(supp);
    return res;
  end Zero_Contribution;

-- INITIALIZATION :

  procedure Initialize 
               ( n : in integer32; mix : in Vector; points : in Array_of_Lists;
                 mixsub : in out Mixed_Subdivision;
                 fs : in out Face_Structures; rest : in out Array_of_Lists;
                 done : in out boolean ) is

  -- DESCRIPTION :
  --   Performs initialization of the dynamic lifting algorithm.

  -- ON ENTRY :
  --   n         length of the vectors in points;
  --   mix       type of mixture;
  --   mixsub    mixed subdivision of the lifted points;
  --   fs        face structures of the lifted points.

  -- ON RETURN :
  --   mixsub    if empty on entry, then it contains the initial mixed
  --             cell, in case the problem is nondegenerate;
  --   rest      rest of point list to be processed;
  --   done      if true, then either the mixed volume equals zero,
  --             or there are no more points to process.

    null_points : Array_of_Lists(points'range);

  begin
    if Is_Null(mixsub) then -- compute an initial mixed cell
      declare
        mic : Mixed_Cell;
      begin
        Initial_Mixed_Cell(n,mix,points,mic,rest);
        if Mixed_Volume(n,mix,mic) > 0 then -- check on degeneracy
         -- initialize the mixed subdivision :
          Construct(mic,mixsub);
         -- initialize the face structures :
          for i in mic.pts'range loop
            declare
              tmp : List := mic.pts(i);
              fc : Face;
            begin
              while not Is_Null(tmp) loop
                Append(fs(i).l,fs(i).last,Head_Of(tmp).all);
                tmp := Tail_of(tmp);
              end loop;
              fc := new VecVec'(Shallow_Create(mic.pts(i)));
              Construct(fc,fs(i).f);
            end;
          end loop;
          if mix'last = mix'first then
            declare
              v : constant VecVec := Shallow_Create(fs(fs'first).l);
              s : constant Simplex := Create(v);
            begin
              fs(fs'first).t := Create(s);
            end;
          end if;
          done := Is_Empty(rest);
        else -- degenerate problem => mixed volume = 0
          rest := null_points; done := true;
        end if;
      end;
    else
      rest := points;
      done := Is_Empty(rest);
    end if;
  end Initialize;

  function Initial_Triangulation
                ( n : in integer32; l : in List; point : in Link_to_Vector )
                return Triangulation is

  -- DESCRIPTION :
  --   Given a list of lifted points with Length_Of(l) >= n and another
  --   lifted point, a nonempty triangulation will be returned, when the
  --   volume of \conv(l U {point}) > 0.

    res : Triangulation;
    del,span : List;
    delpoint : Link_to_Vector := new vector(1..n);

    function Deflate ( lifted : in List ) return List is

    -- DESCRIPTION :
    --   Discards the lifting value of the points in the list l.

      res,res_last : List;
      tmp : List := lifted;
      pt : Link_to_Vector;

    begin
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        declare
          dpt : constant Link_to_Vector := new vector(1..n);
        begin
          dpt.all := pt(dpt'range);
          Append(res,res_last,dpt);
        end;
        tmp := Tail_Of(tmp);
      end loop;
      return res;
    end Deflate;

    function Lift ( pt,liftpt : Link_to_Vector; lifted : List )
                  return Link_to_Vector is

    -- DESCRIPTION :
    --   Compares the value of pt with that of liftpt and
    --   looks up the lifting value of the point pt in the list lifted.
    --   The lifted point will be returned.

    -- REQUIRED :
    --   Either the point pt should be equal to Deflate(liftpt)
    --   or it should belong to Deflate(lifted).

      tmp : List := lifted;
      res : Link_to_Vector;

    begin
      if pt.all = liftpt(pt'range)
       then return liftpt;
       else while not Is_Null(tmp) loop
              res := Head_Of(tmp);
              if pt.all = res(pt'range)
               then return res;
               else tmp := Tail_Of(tmp);
              end if;
            end loop;
            return res;
      end if;
    end Lift;

  begin
    del := Deflate(l);
    delpoint.all := point(delpoint'range);
    Construct(delpoint,del);
    span := Extremal_Points(n,del);
    if Length_Of(span) = natural32(n)+1
     then
       declare
         verpts : VecVec(1..n+1);
         tmp : List := span;
         pt : Link_to_Vector;
         s : Simplex;
       begin
         for i in verpts'range loop
           verpts(i) := Lift(Head_Of(tmp),point,l);
           tmp := Tail_Of(tmp);
         end loop;
         s := Create(verpts);
         res := Create(s);
         tmp := l;
         while not Is_Null(tmp) loop
           pt := Head_Of(tmp);
           if not Is_In(span,pt)
            then Update(res,pt);
           end if;
           tmp := Tail_Of(tmp);
         end loop;
       end;
    end if;
    Clear(del); Clear(span); Clear(delpoint);
    return res;
  end Initial_Triangulation;

-- CHOOSE NEXT POINT :

  procedure Next_Point
                ( n : in integer32; points : in out Array_of_Lists;
                  order : in boolean;
                  index : in out integer32; point : out Link_to_Vector ) is

  -- DESCRIPTION :
  --   Chooses a point out of the lists of points.

  -- ON ENTRY :
  --   n          length of the elements in the point lists;
  --   points     nonempty lists of points;
  --   order      if true, then the first element out of the next first
  --              nonempty list will be chosen;
  --   index      indicates component where not Is_Null(points(index)).

  -- ON RETURN :
  --   points     lists of points with the point removed;
  --   index      points to the next nonempty list,
  --              if index < points'first, then Is_Empty(points);
  --   point      the next point to be processed.

    newindex : integer32 := points'first-1;
    pt : Link_to_Vector;

  begin
    if order
     then pt := Head_Of(points(index));
     else pt := Max_Extreme(points(index),n,-3,3);
          Swap_to_Front(points(index),pt);
    end if;
    points(index) := Tail_Of(points(index));
    for i in (index+1)..points'last loop
      if not Is_Null(points(i))
       then newindex := i;
      end if;
      exit when newindex >= points'first;
    end loop;
    if newindex < points'first
     then for i in points'first..index loop
            if not Is_Null(points(i))
             then newindex := i;
            end if;
            exit when newindex >= points'first;
          end loop;
    end if;
    index := newindex;
    point := pt;
  end Next_Point;

-- UPDATE ROUTINES :

  procedure Merge ( mic : in out Mixed_Cell;
                    mixsub : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;
    done : boolean := false;

  begin
    while not Is_Null(tmp) loop
      declare
        mic1 : constant Mixed_Cell := Head_Of(tmp);
        last : List;
      begin
        if Equal(mic.nor,mic1.nor) then
          for k in mic1.pts'range loop
            last := mic1.pts(k);
            while not Is_Null(Tail_Of(last)) loop
              last := Tail_Of(last);
            end loop;
            Deep_Concat_Diff(mic.pts(k),mic1.pts(k),last);
          end loop;
          done := true;
        else
          tmp := Tail_Of(tmp);
        end if;
      end;
      exit when done;
    end loop;
    if done
     then Deep_Clear(mic);
     else Construct(mic,mixsub);
    end if;
  end Merge;

  procedure Merge ( cells : in Mixed_Subdivision;
                    mixsub : in out Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := cells;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        Merge(mic,mixsub);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Merge;

  procedure Compute_New_Faces
                 ( fs : in out Face_Structure; k,n : in integer32;
                   point : in out Link_to_Vector; newfs : out Faces ) is

  -- DESCRIPTION :
  --   Given a point, the new faces will be computed and returned.

  -- ON ENTRY :
  --   fs          a face structure;
  --   k           dimension of the faces to generate;
  --   n           dimension without the lifting;
  --   point       point that has to be in all the faces.

  -- ON RETURN :
  --   fs          face structure with updated triangulation,
  --               the new faces are not added, also the list is not updated;
  --   point       lifted conservatively w.r.t. fs.t;
  --   newfs       k-faces, which all contain the given point.

    res,res_last : Faces;

    procedure Append ( fc : in VecVec ) is

      f : Face;

    begin
      f := new VecVec'(fc);
      Append(res,res_last,f);
    end Append;
    procedure EnumLis is new Enumerate_Lower_Faces_in_List(Append);
    procedure EnumTri is new Enumerate_Faces_in_Triangulation(Append);

  begin
   -- COMPUTE THE NEW FACES AND UPDATE fs :
    if Is_Null(fs.t) then
      if Length_Of(fs.l) >= natural32(n) then
        fs.t := Initial_Triangulation(n,fs.l,point);
        if Is_Null(fs.t)
         then EnumLis(fs.l,point,k);
         else EnumTri(fs.t,point,k);
        end if;
      else 
        EnumLis(fs.l,point,k);
      end if;
    else 
      declare
        newt : Triangulation;
      begin
        point(point'last) := Lift_to_Place(fs.t,point.all);
        Update(fs.t,point,newt);
        Enumtri(newt,point,k);
      end;
    end if;
    Append(fs.l,fs.last,point);
    newfs := res;
  end Compute_New_Faces;

  procedure Swap ( index : in integer32; points : in out Array_of_Lists ) is

  -- DESCRIPTION :
  --   Swaps the first component of points with the component index.

    tmp : constant List := points(index);

  begin
    points(index) := points(points'first);
    points(points'first) := tmp;
  end Swap;

  procedure Swap ( index : in integer32; mixsub : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Swaps the first component of each cell with the component index.

    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : constant Mixed_Cell := Head_Of(tmp);
      begin
        Swap(index,mic.pts.all);
       -- Set_Head(tmp,mic);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Swap;

  procedure Create_Cells
                ( index,n : in integer32; mix : in Vector;
                  afs : in Array_of_Faces; lif : in Array_of_Lists;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                  mixsub : out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Creates a mixed subdivision by considering the faces 
  --   of component index first.

    res : Mixed_Subdivision;
    wrkmix : Vector(mix'range) := mix;
    wrkafs : Array_of_Faces(afs'range) := afs;
    wrklif : Array_of_Lists(lif'range) := lif;

  begin
    wrkmix(wrkmix'first) := mix(index); wrkmix(index) := mix(mix'first);
    wrkafs(wrkafs'first) := afs(index); wrkafs(index) := afs(afs'first);
    wrklif(wrklif'first) := lif(index); wrklif(index) := lif(lif'first);
    Create_CS(n,wrkmix,wrkafs,wrklif,nbsucc,nbfail,res);
    Swap(index,res);
    mixsub := res;
  end Create_Cells;

  procedure Compute_New_Cells 
                ( n : in integer32; mix : in Vector; mic : in Mixed_Cell;
                  afs : in Array_of_Faces; index : in integer32;
                  lifted : in Array_of_Lists; mixsub : out Mixed_Subdivision;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Computes the new cells by considering a given mixed cell with
  --   a tuple of faces.

  -- ON ENTRY :
  --   n          dimension before lifting;
  --   mix        type of mixture;
  --   mic        a given mixed cell;
  --   afs        type of faces;
  --   index      component where new faces are computed for;
  --   lifted     tuple of lifted points;
  --   nbsucc     number of succesful tests of face-face combinations;
  --   nbfail     number of unsuccesful tests of face-face combinations.

  -- ON RETURN :
  --   mixsub     the new mixed cells;
  --   nbsucc     updated number of succesful tests;
  --   nbfail     updated number of unsuccesful tests;

    res,res2 : Mixed_Subdivision;
    neiafs : Array_of_Faces(afs'range);

  begin
    neiafs(index) := Neighboring_Faces(mic,afs(index),index);
    if not Is_Null(neiafs(index))
     then
       for i in neiafs'range loop
         if i /= index
          then neiafs(i) := Neighboring_Faces(mic,afs(i),i);
         end if;
       end loop;
       if index /= 1
        then Create_Cells(index,n,mix,neiafs,lifted,nbsucc,nbfail,res);
        else Create_CS(n,mix,neiafs,lifted,nbsucc,nbfail,res);
       end if;
       for i in neiafs'range loop
         if i /= index
          then Shallow_Clear(neiafs(i));
         end if;
       end loop;
    end if;
    Shallow_Clear(neiafs(index));
    res2 := Create(lifted,res); -- test on normals
    Shallow_Clear(res);
    mixsub := res2;
  end Compute_New_Cells;

  procedure Compute_New_Cells
               ( n : in integer32; mix : in Vector;
                 mixsub : in Mixed_Subdivision; index : in integer32;
                 newfaces : in Faces; fs : in Face_Structures;
                 newcells : out Mixed_Subdivision; 
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector ) is

    res : Mixed_Subdivision;
    afs : Array_of_Faces(fs'range);
    lifted : Array_of_Lists(fs'range);

  begin
   -- CREATE NEW CELLS ONLY WITH THE NEW FACES :
    for i in afs'range loop
      if i = index
       then afs(index) := newfaces;
       else afs(i) := fs(i).f;
      end if;
      lifted(i) := fs(i).l;
    end loop;
    if not Is_Null(fs(index).t)
     then 
      -- ENUMERATE ALL CELLS IN mixsub :
       declare
         tmp : Mixed_Subdivision := mixsub;
       begin
         while not Is_Null(tmp) loop
           declare
             newcells2 : Mixed_Subdivision;
           begin
             Compute_New_Cells(n,mix,Head_Of(tmp),afs,index,lifted,
                               newcells2,nbsucc,nbfail);
             Merge(newcells2,res); Shallow_Clear(newcells2);
           end;
           tmp := Tail_Of(tmp);
         end loop;
       end;
     else
      -- COMPUTE ALL NEW CELLS AT ONCE :
       declare
         res2 : Mixed_Subdivision;
       begin
         if index /= 1
          then Create_Cells(index,n,mix,afs,lifted,nbsucc,nbfail,res2);
          else Create_CS(n,mix,afs,lifted,nbsucc,nbfail,res2);
         end if;
         res := Create(lifted,res2);
         Shallow_Clear(res2);
       end;
    end if;
    newcells := res;
  end Compute_New_Cells;

-- BASIC VERSION : WITHOUT OUTPUT GENERICS :

  procedure Dynamic_Lifting
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists; 
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector ) is

    rest,inner : Array_of_Lists(points'range);
    index,newindex : integer32;
    finished : boolean := false;
    pt : Link_to_Vector;
    nexli : integer32 := 1;

  begin
    Initialize(n,mix,points,mixsub,fs,rest,finished);
    if not finished
     then 
       index := First_Non_Empty(rest); newindex := index;
       while not finished loop  -- ITERATE UNTIL NO MORE POINTS :
         Next_Point(n,rest,order,newindex,pt);
         declare
           point : Link_to_Vector := new vector(pt'first..pt'last+1);
           newfa : Faces;
           newcells : Mixed_Subdivision;
         begin -- LIFT THE POINT CONSERVATIVELY :
           point(pt'range) := pt.all;
           point(point'last) := 1;
           if inter and then Is_In(fs(index),point.all)
            then
              Clear(point); Construct(pt,inner(index));
            else
              nexli := Conservative_Lifting(mixsub,index,point.all);
              if (maxli > 0) and then nexli > maxli
               then Flatten(mixsub); Flatten(fs);
                    nexli := 1;
              end if;
              point(point'last) := nexli;
             -- COMPUTE NEW FACES AND NEW CELLS :
              Compute_New_Faces(fs(index),mix(index),n,point,newfa);
              if not conmv or else not Zero_Contribution(n,fs,point.all,index)
               then Compute_New_Cells(n,mix,mixsub,index,newfa,fs,
                                      newcells,nbsucc,nbfail);
             -- UPDATE THE MIXED SUBDIVISION AND THE FACE STRUCTURES :
                    Construct(newcells,mixsub);   Shallow_Clear(newcells);
              end if;
              Construct(fs(index).f,newfa); Shallow_Clear(newfa);
           end if;
         end;
         index := newindex;
         finished := (index < rest'first);
       end loop;
    end if;
    if inter then -- lift out the interior points
      for i in inner'range loop
        Constant_Lifting(inner(i),nexli,fs(i).l,fs(i).last);
        Shallow_Clear(inner(i));
      end loop;
    end if;
  exception
    when constraint_error    -- probably due to a too high lifting
       => Flatten(mixsub); Flatten(fs);
          Dynamic_Lifting(n,mix,rest,order,inter,conmv,maxli,mixsub,fs,
                          nbsucc,nbfail);
  end Dynamic_Lifting;

-- EXTENDED VERSIONS : WITH OUTPUT GENERICS :

  procedure Dynamic_Lifting_with_Flat
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists; 
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector ) is

    rest,inner : Array_of_Lists(points'range);
    index,newindex : integer32;
    finished : boolean := false;
    pt : Link_to_Vector;
    nexli : integer32 := 1;

  begin
    Initialize(n,mix,points,mixsub,fs,rest,finished);
    if not finished then
      index := First_Non_Empty(rest); newindex := index;
      while not finished loop -- ITERATE UNTIL NO MORE POINTS :
        Next_Point(n,rest,order,newindex,pt);
        declare
          point : Link_to_Vector := new vector(pt'first..pt'last+1);
          newfa : Faces;
          newcells : Mixed_Subdivision;
        begin -- LIFT THE POINT CONSERVATIVELY :
          point(pt'range) := pt.all;
          point(point'last) := 1;
          if inter and then Is_In(fs(index),point.all) then
            Clear(point); Construct(pt,inner(index));
          else
            nexli := Conservative_Lifting(mixsub,index,point.all);
            if (maxli > 0) and then nexli > maxli then
              Before_Flattening(mixsub,fs); 
              Flatten(mixsub); Flatten(fs);
              nexli := 1;
            end if;
            point(point'last) := nexli;
           -- COMPUTE NEW FACES AND NEW CELLS :
            Compute_New_Faces(fs(index),mix(index),n,point,newfa);
            if not conmv or else not Zero_Contribution(n,fs,point.all,index)
             then Compute_New_Cells(n,mix,mixsub,index,newfa,fs,
                                    newcells,nbsucc,nbfail);
           -- UPDATE THE MIXED SUBDIVISION AND THE FACE STRUCTURES :
                  Construct(newcells,mixsub);   Shallow_Clear(newcells);
            end if;
            Construct(fs(index).f,newfa); Shallow_Clear(newfa);
          end if;
        end;
        index := newindex;
        finished := (index < rest'first);
      end loop;
    end if;
    if inter then -- lift out the interior points
      for i in inner'range loop
        Constant_Lifting(inner(i),nexli,fs(i).l,fs(i).last);
        Shallow_Clear(inner(i));
      end loop;
    end if;
  exception
    when constraint_error    -- probably due to a too high lifting
       => Before_Flattening(mixsub,fs);
          Flatten(mixsub); Flatten(fs);
          Dynamic_Lifting_with_Flat(n,mix,rest,order,inter,conmv,maxli,mixsub,
                                    fs,nbsucc,nbfail);
  end Dynamic_Lifting_with_Flat;

  procedure Dynamic_Lifting_with_New
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists; 
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector ) is

    rest,inner : Array_of_Lists(points'range);
    index,newindex : integer32;
    finished : boolean := false;
    pt : Link_to_Vector;
    nexli : integer32 := 1;

  begin
    Initialize(n,mix,points,mixsub,fs,rest,finished);
    Process_New_Cells(mixsub,0,Head_Of(fs(fs'last).last).all);
    if not finished then
      index := First_Non_Empty(rest); newindex := index;
      while not finished loop -- ITERATE UNTIL NO MORE POINTS :
        Next_Point(n,rest,order,newindex,pt);
        declare
          point : Link_to_Vector := new vector(pt'first..pt'last+1);
          newfa : Faces;
          newcells : Mixed_Subdivision;
        begin -- LIFT THE POINT CONSERVATIVELY :
          point(pt'range) := pt.all;
          point(point'last) := 1;
          if inter and then Is_In(fs(index),point.all) then
            Clear(point); Construct(pt,inner(index));
          else
            nexli := Conservative_Lifting(mixsub,index,point.all);
            if (maxli > 0) and then nexli > maxli then
              Flatten(mixsub); Flatten(fs);
              nexli := 1;
            end if;
            point(point'last) := nexli;
           -- COMPUTE NEW FACES AND NEW CELLS :
            Compute_New_Faces(fs(index),mix(index),n,point,newfa);
            if not conmv or else not Zero_Contribution(n,fs,point.all,index)
             then Compute_New_Cells(n,mix,mixsub,index,newfa,fs,newcells,
                                    nbsucc,nbfail);
                  Process_New_Cells(newcells,index,point.all);
           -- UPDATE THE MIXED SUBDIVISION AND THE FACE STRUCTURES :
                  Construct(newcells,mixsub);   Shallow_Clear(newcells);
            end if;
            Construct(fs(index).f,newfa); Shallow_Clear(newfa);
          end if;
        end;
       index := newindex;
       finished := (index < rest'first);
     end loop;
    end if;
    if inter then -- lift out the interior points
      for i in inner'range loop
        Constant_Lifting(inner(i),nexli,fs(i).l,fs(i).last);
        Shallow_Clear(inner(i));
      end loop;
    end if;
  exception
    when constraint_error    -- probably due to a too high lifting
       => Flatten(mixsub); Flatten(fs);
          Dynamic_Lifting_with_New(n,mix,rest,order,inter,conmv,maxli,mixsub,fs,
                                   nbsucc,nbfail);
  end Dynamic_Lifting_with_New;

  procedure Dynamic_Lifting_with_Flat_and_New
                ( n : in integer32; mix : in Vector;
                  points : in Array_of_Lists; 
                  order,inter,conmv : in boolean; maxli : in integer32;
                  mixsub : in out Mixed_Subdivision;
                  fs : in out Face_Structures;
                  nbsucc,nbfail : in out Standard_Floating_Vectors.Vector ) is

    rest,inner : Array_of_Lists(points'range);
    index,newindex : integer32;
    finished : boolean := false;
    pt : Link_to_Vector;
    nexli : integer32 := 1;

  begin
    Initialize(n,mix,points,mixsub,fs,rest,finished);
    Process_New_Cells(mixsub,0,Head_Of(fs(fs'last).last).all);
    if not finished then
      index := First_Non_Empty(rest); newindex := index;
      while not finished loop -- ITERATE UNTIL NO MORE POINTS 
        Next_Point(n,rest,order,newindex,pt);
        declare
          point : Link_to_Vector := new vector(pt'first..pt'last+1);
          newfa : Faces;
          newcells : Mixed_Subdivision;
        begin -- LIFT THE POINT CONSERVATIVELY :
          point(pt'range) := pt.all;
          point(point'last) := 1;
          if inter and then Is_In(fs(index),point.all) then
            Clear(point); Construct(pt,inner(index));
          else
            nexli := Conservative_Lifting(mixsub,index,point.all);
            if (maxli > 0) and then nexli > maxli then
              Before_Flattening(mixsub,fs);
              Flatten(mixsub); Flatten(fs);
              nexli := 1;
            end if;
            point(point'last) := nexli;
           -- COMPUTE NEW FACES AND NEW CELLS :
            Compute_New_Faces(fs(index),mix(index),n,point,newfa);
            if not conmv or else not Zero_Contribution(n,fs,point.all,index)
             then Compute_New_Cells(n,mix,mixsub,index,newfa,fs,newcells,
                                    nbsucc,nbfail);
                  Process_New_Cells(newcells,index,point.all);
           -- UPDATE THE MIXED SUBDIVISION AND THE FACE STRUCTURES :
                  Construct(newcells,mixsub);   Shallow_Clear(newcells);
            end if;
            Construct(fs(index).f,newfa); Shallow_Clear(newfa);
          end if;
        end;
        index := newindex;
        finished := (index < rest'first);
      end loop;
    end if;
    if inter then -- lift out the interior points
      for i in inner'range loop
        Constant_Lifting(inner(i),nexli,fs(i).l,fs(i).last);
        Shallow_Clear(inner(i));
      end loop;
    end if;
  exception
    when constraint_error    -- probably due to a too high lifting
       => Before_Flattening(mixsub,fs);
          Flatten(mixsub); Flatten(fs);
          Dynamic_Lifting_with_Flat_and_New
              (n,mix,rest,order,inter,conmv,maxli,mixsub,fs,nbsucc,nbfail);
  end Dynamic_Lifting_with_Flat_and_New;

end Dynamic_Mixed_Subdivisions;
