with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;
with Standard_Integer_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Planes_and_Polynomials;             use Planes_and_Polynomials;
with Witness_Sets;                       use Witness_Sets;
with Sampling_Machine;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Sample_Point_Grids;                 use Sample_Point_Grids;
with Standard_Linear_Projections;        use Standard_Linear_Projections;
with Standard_Central_Projections;       use Standard_Central_Projections;
with Standard_Subspace_Restrictions;     use Standard_Subspace_Restrictions;
with Standard_Polynomial_Interpolators;  use Standard_Polynomial_Interpolators;
with Standard_Membership_Tests;          use Standard_Membership_Tests;

package body Standard_Breakup_Components is

-- AUXILIARY DATA STRUCTURE :

  type Boolean_Array is array ( integer32 range <> ) of boolean;

-- AUXILIARY TO CREATE GRID STRUCTURE :

  function Create ( len,dim : integer32 )
                  return Standard_Complex_VecVecs.Array_of_VecVecs is

    use Standard_Complex_VecVecs;

    res : Array_of_VecVecs(1..len);

  begin
    for i in 1..len loop
      res(i) := new VecVec(1..dim);
    end loop;
    return res;
  end Create;

-- TOOLS TO MANAGE INTERPOLATORS :

  procedure Count_Components
                ( file : in file_type; p : in out Link_to_Poly_Sys;
                  tol : in double_float ) is

    dis : double_float;
    found : array(p'range) of boolean := (p'range => false);
    cnt : integer32;

  begin
    for i in p'range loop                        -- compute mutual distances
      if not found(i) then
        for j in i+1..p'last loop
          if not found(j) then
            if Degree(p(i)) = Degree(p(j)) then
              dis := Distance(p(i),p(j));
              put(file,"Distance between p("); put(file,i,1);
              put(file,") and p("); put(file,j,1);
              put(file,") is "); put(file,dis,3); new_line(file);
              if dis < tol
               then found(j) := true;
              end if;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    cnt := 0;                               -- count the number of components
    for i in found'range loop
      if not found(i)
       then cnt := cnt+1;
      end if;
    end loop;
    put(file,"Found "); put(file,cnt,1);
    put(file," component");
    if cnt = 1
     then put_line(file,".");
     else put_line(file,"s.");
    end if;
    declare                                         -- collect the equations
      equ : Poly_Sys(1..cnt);
      ind : integer32 := 0;
    begin
      for i in found'range loop
        if not found(i) then
          ind := ind+1;
          Copy(p(i),equ(ind));
        end if;
      end loop;
      Clear(p);
      p := new Poly_Sys'(equ);
    end;
  end Count_Components;

  function Collect_Components ( p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);
    cnt : integer32 := 0;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        cnt := cnt+1;
        res(cnt) := p(i);
      end if;
    end loop;
    return res(1..cnt);
  end Collect_Components;

-- STOP TEST IN MASSIVE INTERPOLATE :

  procedure Residual_Stop_Test
                ( file : in file_type; intpols : in Poly_Sys;
                  points : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  tol : in double_float;
                  complete : in out Boolean_Array; full_stop : out boolean ) is

    y : Standard_Complex_Vectors.Vector(intpols'range);
    tmp,sum : double_float;
    no_eval : boolean;
    use Standard_Complex_VecVecs;

  begin
    put_line(file,"The interpolating polynomials : ");
    put_line(file,intpols);
    put_line(file,"The results of the evaluation : ");
    for i in points'range loop
      if points(i) /= null then
        put(file,"Evaluating sequence "); put(file,i,1);
        put_line(file," :");
        sum := 0.0;
        no_eval := true;
        for j in points(i)'range loop
          if points(i)(j) /= null then
            no_eval := false;
            y := Eval(intpols,points(i)(j).all);
            put(file,"log10(|y("); put(file,j,1); put(file,")|) :");
            for k in y'range loop
              tmp := AbsVal(y(k));
              if tmp > 10.0E-100
               then put(file,Log10(tmp),4,1,0);
               else put(file,"-100");
              end if;
            end loop;
            new_line(file);
            sum := sum + AbsVal(y(i));
          end if;
        end loop;
        if not no_eval
         then complete(i) := (sum < double_float(points'last)*tol);
        end if;
      end if;
    end loop;
    full_stop := true;
    for i in complete'range loop
      if not complete(i)
       then full_stop := false;
      end if;
      exit when not full_stop;
    end loop;
  end Residual_Stop_Test;

-- INTERPOLATORS :

  function Massive_Interpolate
             ( file : file_type; embsys : Poly_Sys; sols : Solution_List;
               hyp : Standard_Complex_VecVecs.VecVec; level : natural32 )
             return Poly_Sys is

    n : constant integer32 := embsys'last-integer32(level);
    len : constant integer32 := integer32(Length_Of(sols)); -- #generic points
    res : Poly_Sys(1..len);
    intpoly : Poly;                                -- template for interpolator
    timer : Timing_Widget;
    genpt : Standard_Complex_Vectors.Vector(1..n);
    projpoints : Standard_Complex_VecVecs.Array_of_VecVecs(1..len);
    npts : natural32 := Number_of_Terms(natural32(len),level+1) - 1;
    tol : constant double_float := 1.0E-8;
    stop : Boolean_Array(1..len) := (1..len => false);
    full_stop : boolean := false;
    slihyp : constant Standard_Complex_VecVecs.VecVec := Slices(embsys,level);
    samples : constant Standard_Sample_List := Create(sols,slihyp);
    ptrlist : Standard_Sample_List := Create(sols,slihyp);
    grid,grid_last,ptrgrid : Standard_Sample_Grid;

  begin
    tstart(timer);
    Sampling_Machine.Initialize(embsys);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    Sample(samples,npts,grid,grid_last);
    Sampling_Machine.Clear;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Sampling Generic Points");
    new_line(file);
    tstart(timer);
    ptrgrid := grid;
    projpoints := Create(len,integer32(npts));
    for i in projpoints'range loop
      ptrlist := Head_Of(ptrgrid);
      for j in projpoints(i)'range loop
        genpt := Sample_Point(Head_Of(ptrlist)).v;
        projpoints(i)(j) := new Standard_Complex_Vectors.Vector'
          (Evaluate(hyp,genpt,integer32(level)+1));
        ptrlist := Tail_Of(ptrlist);
      end loop;
      ptrgrid := Tail_Of(ptrgrid);
    end loop;
    for d in 1..len loop                            -- increasing degrees 
      intpoly := Create(natural32(d),level+1,1); 
      npts := Number_of_Terms(intpoly) - 1;
      for i in projpoints'range loop
        if not stop(i)
         then res(i) := Interpolate(intpoly,projpoints(i)(1..integer32(npts)));
        end if;
      end loop;
      Residual_Stop_Test(file,res,projpoints,tol,stop,full_stop);
      Clear(intpoly);
      exit when full_stop;
      for i in res'range loop
        if not stop(i)
         then Clear(res(i));
        end if;
      end loop;
    end loop;
    Standard_Complex_VecVecs.Clear(projpoints);
    Deep_Clear(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Interpolating through Generic Points");
    new_line(file);
    return res;
  end Massive_Interpolate;

  procedure Incremental_Interpolator
                ( file : in file_type; sol : in Solution;
                  slihyp,projhyp : in Standard_Complex_VecVecs.VecVec;
                  maxdeg,level : in natural32; tol : in double_float;
                  filter : out Poly ) is

  -- DESCRIPTION :
  --   This is the internal routine in the main loop of the incremental
  --   interpolation method to break up components of solutions.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results;
  --   sol        a generic point to start the sampling;
  --   projhyp    general hyperplanes to project onto;
  --   maxdeg     maximal degree for the interpolating polynomial;
  --   level      number of slices added;
  --   tol        tolerance for residual in stop test.

  -- ON RETURN :
  --   filter     interpolating polynomial.

    intpoly : Poly;
    resi : double_float;
    stop : boolean;
    maxterms : constant natural32 := Number_of_Terms(level+1,maxdeg);
    projpts : Standard_Complex_VecVecs.VecVec(1..integer32(maxterms)+1);
    genpt : Standard_Complex_Vectors.Vector(1..sol.n);
    spt : constant Standard_Sample := Create(sol,slihyp);
    sps,sps_last,prev,tmp : Standard_Sample_List;
    npts,prevnpts,startind : integer32;

  begin
    prevnpts := 1;
    Append(sps,sps_last,spt);
    prev := sps;
    for d in 1..maxdeg loop                        -- for increasing degrees
      intpoly := Create(d,level+1,1);
      npts := integer32(Number_of_Terms(intpoly)); -- oversampling with 1
      Sample(spt,natural32(npts-prevnpts),sps,sps_last);
      if prevnpts = 1 then
        tmp := prev;
        startind := prevnpts;
      else
        tmp := Tail_Of(prev);
        startind := prevnpts+1;
      end if;
      for j in startind..npts loop
        genpt := Sample_Point(Head_Of(tmp)).v;
        projpts(j) := new Vector'(Evaluate(projhyp,genpt,integer32(level)+1));
        tmp := Tail_Of(tmp);
      end loop;
      filter := Interpolate(intpoly,projpts(1..npts-1));
      resi := AbsVal(Eval(filter,projpts(npts).all));
      put(file,"Residual in stop test of interpolator : ");
      put(file,resi,3);
      stop := (resi < tol);
      if stop
       then put(file,"  succeeded at degree ");
       else put(file,"  failed at degree ");
      end if;
      put(file,d,1); new_line(file);
      Clear(intpoly);
      exit when stop or (d = maxdeg); -- even if stop test fails, keep filter
      Clear(filter);
      prevnpts := npts;
      prev := sps_last;
    end loop;
    for i in 1..npts loop
      Standard_Complex_Vectors.Clear(projpts(i));
    end loop;
    Deep_Clear(sps);
  end Incremental_Interpolator;

  function Incremental_Interpolate
             ( file : file_type; embsys : Poly_Sys; sols : Solution_List;
               hyp : Standard_Complex_VecVecs.VecVec; level : natural32 )
             return Poly_Sys is

    len : constant integer32 := integer32(Length_Of(sols)); -- #generic points
    res : Poly_Sys(1..len);
    timer : Timing_Widget;
    tol : constant double_float := 10.0**(-8);
    tmp : Solution_List := sols;
    slihyp : constant Standard_Complex_VecVecs.VecVec := Slices(embsys,level);

  begin
    tstart(timer);
    Sampling_Machine.Initialize(embsys);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..len loop                  -- sequences start at generic points
      if not On_Component(file,res,i,sols,hyp,integer32(level),tol) then
        Incremental_Interpolator
          (file,Head_Of(tmp).all,slihyp,hyp,natural32(len),level,tol,res(i));
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Sampling + Interpolating through Generic Points");
    new_line(file);
    return Collect_Components(res);
  end Incremental_Interpolate;

  procedure Update_Base_Points
                ( file : in file_type;
                  spt : in Standard_Sample;
                  base : in out Standard_Complex_VecVecs.Link_to_VecVec;
                  basecard : in out integer32; prevnpts : in integer32;
                  restpts : in Standard_Complex_VecVecs.Array_of_VecVecs;
                  hyp : in Standard_Complex_VecVecs.VecVec;
                  projpts : in out Standard_Complex_VecVecs.Array_of_VecVecs;
                  dim,level,nbpt : in integer32;
                  restsols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Updates the list of base points by sampling one point and 
  --   computing the projections again of the restricted points;

    use Standard_Complex_VecVecs;
    newspt : Standard_Sample;

  begin
    Sample(spt,newspt);
    if base = null then
      base := new Standard_Complex_VecVecs.VecVec(1..dim-level+1);
      base(1) := new Standard_Complex_Vectors.Vector'(Sample_Point(newspt).v);
      basecard := 1;
      for j in 1..prevnpts-1 loop
        Standard_Complex_Vectors.Clear(projpts(nbpt)(j));              
        projpts(nbpt)(j)
          := new Standard_Complex_Vectors.Vector'
                  (Intersect(hyp(hyp'first).all,base(1).all,
                             restpts(nbpt)(j).all,level+1));
      end loop;
    else
      basecard := basecard+1;
      base(basecard)
        := new Standard_Complex_Vectors.Vector'(Sample_Point(newspt).v);
      declare
        proj : Standard_Complex_Vectors.Vector(base(basecard)'range);
      begin
        proj := Intersect(hyp(1..basecard-1),base(1..basecard-1),
                          base(basecard).all,proj'last);
        Standard_Complex_Vectors.Clear(base(basecard));
        base(basecard) := new Standard_Complex_Vectors.Vector'(proj);
      end;
      for j in 1..prevnpts-1 loop
        Standard_Complex_Vectors.Clear(projpts(nbpt)(j));              
        projpts(nbpt)(j) := new Standard_Complex_Vectors.Vector'
          (Intersect(hyp(1..basecard),base(1..basecard),
                     restpts(nbpt)(j).all,level+1));
      end loop;
    end if;
  end Update_Base_Points;

  procedure Dynamic_Incremental_Interpolator
                ( file : in file_type; embsys : in Poly_Sys;
                  sols : in Solution_List; sol : in Solution;
                  hyp : in Standard_Complex_VecVecs.VecVec;
                  nbpt,maxdeg,level : in integer32; tol : in double_float;
                  genpts,projpts,restpts : 
                    in out Standard_Complex_VecVecs.Array_of_VecVecs;
                  centproj : in boolean; subspace : out Link_to_Poly_Sys;
                  pivots : out Standard_Integer_Vectors.Link_to_Vector;
                  base : out Standard_Complex_VecVecs.Link_to_VecVec;
                  basecard : in out integer32; filter : out Poly ) is

  -- DESCRIPTION :
  --   This is the internal routine in the main loop of the incremental
  --   interpolation method to break up components of solutions.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results;
  --   embsys     embedded polynomial system;
  --   sols       list of generic points as solutions to embsys;
  --   sol        current generic point;
  --   hyp        general hyperplanes to project onto;
  --   nbpt       index to current generic point;
  --   maxdeg     maximal degree for the interpolating polynomial;
  --   level      number of slices added;
  --   tol        tolerance for residual in stop test;
  --   genpts     structure to contain all generic points;
  --   projpts    grid with the projections of the generic points;
  --   restpts    grid with the restrictions of the generic points;
  --   centproj   indicates whether central projections will be used;
  --   bascard    counts the number of base points.

  -- ON RETURN :
  --   genpts     updated samples;
  --   projpts    updated projected points;
  --   restpts    updated restricted points;
  --   subspace   equations for the container subspace, empty when the
  --              container subspace is the whole space;
  --   pivots     remaining variables, empty when subspace is empty;
  --   base       base points used in the projection;
  --   basecard   updated number of base points;
  --   filter     interpolating polynomial.

    n : constant integer32 := embsys'last - level;
    prevnpts,npts,k : integer32;
    intpoly : Poly;
    stop,restricted,skewline : boolean;
    restembsys,restemb : Link_to_Poly_Sys;
    dim : integer32 := n+1;
    restsols : Solution_List;
    resi : double_float;
    slihyp : Standard_Complex_VecVecs.VecVec(1..level)
           := Slices(embsys,natural32(level));
    spt : Standard_Sample := Create(sol,slihyp);
    sps,sps_last,tmp,prev : Standard_Sample_List;
    startind : integer32 := 1;

  begin
    prevnpts := 1;
    restricted := false;
    skewline := false;
    Sampling_Machine.Initialize(embsys);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    Append(sps,sps_last,spt);
    prev := sps;
    for d in 1..maxdeg loop                        -- for increasing degrees
      if restricted then
        if centproj then
          if not skewline and then (dim > level+1) then
            skewline := true;
            Update_Base_Points(file,spt,base,basecard,prevnpts,restpts,
                               hyp,projpts,dim,level,nbpt,restsols);
          else
            if dim - basecard > level+1 then
              Update_Base_Points
                (file,spt,base,basecard,prevnpts,restpts,
                 hyp,projpts,dim,level,nbpt,restsols);
            end if;
          end if;
        end if;
        if skewline
         then intpoly := Create(natural32(d-basecard),natural32(level)+1,1);
         else intpoly := Create(natural32(d),natural32(level)+1,1);
        end if;
        npts := integer32(Number_of_Terms(intpoly));   -- oversampling with 1
        if prevnpts <= npts then
          Sample(spt,natural32(npts-prevnpts),sps,sps_last);
          tmp := Tail_Of(prev);
          for j in prevnpts+1..npts loop
            restpts(nbpt)(j) := new Vector'(Sample_Point(Head_Of(tmp)).v);
            tmp := Tail_Of(tmp);
          end loop;
        end if;
        if skewline then
          for j in startind..npts loop
            projpts(nbpt)(j) := new Standard_Complex_Vectors.Vector'
              (Intersect(hyp(1..basecard),base(1..basecard),
                         restpts(nbpt)(j).all,level+1));
          end loop;
        else
          for j in startind..npts loop
            projpts(nbpt)(j) := new Vector'
              (Evaluate(hyp,restpts(nbpt)(j).all,level+1));
          end loop;
        end if;
      else
        intpoly := Create(natural32(d),natural32(level)+1,1);
        npts := integer32(Number_of_Terms(intpoly));  -- oversampling with 1
        Sample(spt,natural32(npts-prevnpts),sps,sps_last);
        if prevnpts = 1
         then tmp := prev;          startind := prevnpts;
         else tmp := Tail_Of(prev); startind := prevnpts+1;
        end if;
        for j in startind..npts loop
          genpts(nbpt)(j) := new Vector'(Sample_Point(Head_Of(tmp)).v);
          tmp := Tail_Of(tmp);
        end loop;
        k := npts-1;
        if k > n
         then k := n;
        end if;
        Subspace_Restriction
          (file,embsys,natural32(nbpt),natural32(k),natural32(n),
           natural32(level),genpts,restpts,pivots,subspace,restembsys,
           natural32(dim));
        if dim < k then
          restricted := true;
          restemb := new Poly_Sys'(Collapse_Equations
            (restembsys.all,natural32(dim),natural32(level)));
          -- put_line(file,"The collapsed equations : ");
          -- put_line(file,restemb.all);
          Sampling_Machine.Clear;
          Sampling_Machine.Initialize(restemb.all);
          restsols := Restrict_Solution(sols,nbpt,level,pivots.all);
          spt := Create(Retrieve(restsols,natural32(nbpt)).all,slihyp);
          -- put_line(file,"The restricted solution : ");
          -- put(file,Get(restsols,nbpt)); new_line(file);
          -- declare
          --   s : Solution(restemb'last) := Get(restsols,nbpt);
          --   y : Vector(restemb'range) := Eval(restemb.all,s.v);
          -- begin
          --   put_line(file,"Evaluation of restrictions : ");
          --   put_line(file,y);
          -- end;
          for j in 1..startind-1 loop  -- recompute projections
            Clear(projpts(nbpt)(j));
            projpts(nbpt)(j) := new Vector'
              (Evaluate(hyp,restpts(nbpt)(j).all,level+1));
          end loop;
          for j in startind..npts loop
            projpts(nbpt)(j) := new Vector'
              (Evaluate(hyp,restpts(nbpt)(j).all,level+1));
          end loop;
        else
          for j in startind..npts loop
            projpts(nbpt)(j) := new Vector'
              (Evaluate(hyp,genpts(nbpt)(j).all,level+1));
          end loop;
        end if;
      end if;
      if restricted and dim = 1
       then Copy(intpoly,filter);
            stop := true;
       else filter := Interpolate(intpoly,projpts(nbpt)(1..npts-1));
            resi := AbsVal(Eval(filter,projpts(nbpt)(npts).all));
            put(file,"Residual in stop test of interpolator : ");
            put(file,resi,3);
            stop := (resi < tol);
            if stop
             then put(file,"  succeeded at degree ");
             else put(file,"  failed at degree ");
            end if;
            put(file,d,1); new_line(file);
      end if;
      Clear(intpoly);
      exit when stop or (d = maxdeg); -- even if stop test fails, keep filter
      Clear(filter);
      prevnpts := npts;
      prev := sps_last;
    end loop;
    Deep_Clear(sps);
    Sampling_Machine.Clear;
  end Dynamic_Incremental_Interpolator;

  procedure Dynamic_Interpolate
                ( file : in file_type; embsys : in Poly_Sys;
                  level : in natural32; sols : in Solution_List;
                  hyp : in Standard_Complex_VecVecs.VecVec;
                  centproj : in boolean; subspaces : out Array_of_Poly_Sys;
                  pivots : out Standard_Integer_VecVecs.VecVec;
                  basepts : out Standard_Complex_VecVecs.Array_of_VecVecs;
                  basecard : out Standard_Natural_Vectors.Vector;
                  filters : out Link_to_Poly_Sys ) is

    len : constant integer32 := integer32(Length_Of(sols)); -- #generic points
    res : Poly_Sys(1..len);
    timer : Timing_Widget;
    maxpts : constant natural32
           := Number_of_Terms(natural32(len),level+1); -- oversampling with 1;
    genpts,projpts,restpts : Standard_Complex_VecVecs.Array_of_VecVecs(1..len);
    tol : constant double_float := 1.0E-8;
    tmp : Solution_List := sols;

  begin
    tstart(timer);
    subspaces(subspaces'first) := null;   -- only to avoid warnings ...
    pivots(pivots'first) := null;
    basepts(basepts'first) := null;
    genpts := Create(len,integer32(maxpts));
    projpts := Create(len,integer32(maxpts));
    restpts := Create(len,integer32(maxpts));
    for i in 1..len loop                  -- sequences start at generic points
      basecard(i) := 0;
      if not On_Component(file,res,i,subspaces,pivots,basepts,basecard,
                          sols,hyp,integer32(level),tol)
       then Dynamic_Incremental_Interpolator
              (file,embsys,sols,Head_Of(tmp).all,hyp,i,len,integer32(level),
               tol,genpts,projpts,restpts,centproj,subspaces(i),pivots(i),
               basepts(i),integer32(basecard(i)),res(i));
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    Standard_Complex_VecVecs.Clear(genpts);
    Standard_Complex_VecVecs.Clear(projpts);
    Standard_Complex_VecVecs.Clear(restpts);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Sampling + Interpolating through Generic Points");
    new_line(file);
   -- WARNING : THIS COULD LEAD TO INCONSISTENCIES WITH SUBSPACES
   -- filters := new Poly_Sys'(Collect_Components(res));
    filters := new Poly_Sys'(res);
  end Dynamic_Interpolate;

end Standard_Breakup_Components;
