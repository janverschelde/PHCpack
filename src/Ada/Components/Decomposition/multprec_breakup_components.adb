with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Standard_Integer_Vectors;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vector_Tools;     use Multprec_Complex_Vector_Tools;
with Standard_Complex_VecVecs;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Functions;   use Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Systems_io;  use Multprec_Complex_Poly_Systems_io;
with Multprec_Complex_Solutions;
with Planes_and_Polynomials;            use Planes_and_Polynomials;
with Witness_Sets;                      use Witness_Sets;
with Sampling_Machine;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Sample_Point_Grids;                use Sample_Point_Grids;
with Multprec_Linear_Projections;       use Multprec_Linear_Projections;
with Multprec_Central_Projections;      use Multprec_Central_Projections;
with Multprec_Subspace_Restrictions;    use Multprec_Subspace_Restrictions;
with Multprec_Polynomial_Interpolators; use Multprec_Polynomial_Interpolators;
with Multprec_Membership_Tests;         use Multprec_Membership_Tests;

package body Multprec_Breakup_Components is

-- DATA STRUCTURE :

  type Boolean_Array is array ( integer32 range <> ) of boolean;

-- AUXILIARY for intermediate structures 

  function Create ( len,dim : integer32 )
                  return Multprec_Complex_VecVecs.Array_of_VecVecs is

    use Multprec_Complex_VecVecs;

    res : Array_of_VecVecs(1..len);

  begin
    for i in 1..len loop
      res(i) := new VecVec(1..dim);
    end loop;
    return res;
  end Create;

  procedure Write_Diagnostics
              ( file : in file_type;
                s : in Multprec_Complex_Solutions.Solution ) is

  -- DESCRIPTION :
  --   Writes the diagnostics about the solution on file.

  begin
    put(file,"==");
    put(file," err : "); put(file,s.err,2,3,3); put(file," =");
    put(file," rco : "); put(file,s.rco,2,3,3); put(file," =");
    put(file," res : "); put(file,s.res,2,3,3); put_line(file," =");
  end Write_Diagnostics;

  maxsamples : constant natural32 := 500;
   -- curve of degree 30 needs 496 samples

  function Minimum ( a,b : natural32 ) return natural32 is
  begin
    if a < b
     then return a;
     else return b;
    end if;
  end Minimum;

-- TOOLS TO MANAGE INTERPOLATORS :

  function Collect_Components
              ( p : Multprec_Complex_Poly_Systems.Poly_Sys )
              return Multprec_Complex_Poly_Systems.Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Poly_Sys(p'range);
    cnt : integer32 := 0;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        cnt := cnt + 1;
        res(cnt) := p(i);
      end if;
    end loop;
    return res(1..cnt);
  end Collect_Components;

-- STOP TEST IN MASSIVE INTERPOLATOR :

  procedure Residual_Stop_Test
              ( file : in file_type;
                intpols : in Multprec_Complex_Poly_Systems.Poly_Sys;
                points : in Multprec_Complex_VecVecs.Array_of_VecVecs;
                tol : in double_float; complete : in out Boolean_Array; 
                full_stop : out boolean ) is

    y : Multprec_Complex_Vectors.Vector(intpols'range);
    fltacc : Floating_Number;
    tmp,sum : double_float;
    no_eval : boolean;

    use Multprec_Complex_Vectors;
    use Multprec_Complex_VecVecs;

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
            for k in y'range loop
              if intpols(k) /= Null_Poly
               then y(k) := Eval(intpols(k),points(i)(j).all);
              end if;
            end loop;
            put(file,"|y("); put(file,j,1); put(file,")| :");
            for k in y'range loop
              fltacc := AbsVal(y(k));
              tmp := Trunc(fltacc);
              put(file,tmp,4,1,3);
              Clear(fltacc);
            end loop;
            new_line(file);
            fltacc := AbsVal(y(i));
            sum := sum + Trunc(fltacc);
            Clear(fltacc);
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
    Multprec_Complex_Vectors.Clear(y);
  end Residual_Stop_Test;

-- THE INTERPOLATORS :

  function Massive_Interpolate
             ( file : file_type;
                embsys : Standard_Complex_Poly_Systems.Poly_Sys;
                orgsys : Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : Standard_Complex_Solutions.Solution_List;
                hyp : Multprec_Complex_VecVecs.VecVec;
                level,size : natural32 )
              return Multprec_Complex_Poly_Systems.Poly_Sys is

    n : constant integer32 := embsys'last - integer32(level);
    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
      -- #generic points
    res : Multprec_Complex_Poly_Systems.Poly_Sys(1..len);
    intpoly : Poly;                              -- template for interpolator
    timer : Timing_Widget;
    genpt : Multprec_Complex_Vectors.Vector(1..n);
    projpts : Multprec_Complex_VecVecs.Array_of_VecVecs(1..len);
    npts : natural32 := Number_of_Terms(natural32(len),level+1) - 1;
    tol : constant double_float := 1.0E-8;
    invcond : Floating_Number;
    stop : Boolean_Array(1..len) := (1..len => false);
    full_stop : boolean := false;
    slihyp : constant Standard_Complex_VecVecs.VecVec := Slices(embsys,level);
    samples : constant Standard_Sample_List := Create(sols,slihyp);
    ptrlist : Multprec_Sample_List;
    grid,grid_last,ptrgrid : Multprec_Sample_Grid;

  begin
    tstart(timer);
    Sampling_Machine.Initialize(embsys,orgsys,integer32(level),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner(size);
    Sample(samples,npts,grid,grid_last);
    Sampling_Machine.Clear;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Sampling Generic Points");
    new_line(file);
    tstart(timer);
    ptrgrid := grid;
    projpts := Create(len,integer32(npts));
    for i in projpts'range loop
      ptrlist := Head_Of(ptrgrid);
      for j in projpts(i)'range loop
        genpt := Sample_Point(Head_Of(ptrlist)).v;
        projpts(i)(j) := new Multprec_Complex_Vectors.Vector'
          (Evaluate(hyp,genpt,integer32(level)+1));
        ptrlist := Tail_Of(ptrlist);
      end loop;
      ptrgrid := Tail_Of(ptrgrid);
    end loop;
    for d in 1..len loop                            -- increasing degrees 
      intpoly := Create(natural32(d),level+1,1); 
      npts := Number_of_Terms(intpoly) - 1;
      for i in projpts'range loop
        if not stop(i) then
          Clear(res(i));
          Interpolate(intpoly,projpts(i)(1..integer32(npts)),res(i),invcond);
          put(file,"Estimate for inverse condition number : ");
          put(file,invcond,3); new_line(file);
          Clear(invcond);
        end if;
      end loop;
      Residual_Stop_Test(file,res,projpts,tol,stop,full_stop);
      Clear(intpoly);
      exit when full_stop;
      for i in res'range loop
        if not stop(i)
         then Clear(res(i));
        end if;
      end loop;
    end loop;
    Multprec_Complex_VecVecs.Clear(projpts);
    Deep_Clear(grid);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Interpolating through generic points");
    new_line(file);
    return res;
  end Massive_Interpolate;

  procedure Incremental_Interpolator
                ( file : in file_type; mpspt : in Multprec_Sample;
                  prohyp : in Multprec_Complex_VecVecs.VecVec;
                  maxdeg,level : in natural32;
                  tol : in double_float; filter : out Poly ) is

  -- DESCRIPTION :
  --   This is the internal routine in the main loop of the incremental
  --   interpolation method to break up components of solutions.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results;
  --   mpspt      refinement of the sample point (sol,hyp);
  --   prohyp     general hyperplanes to project onto;
  --   maxdeg     maximal degree for the interpolating polynomial;
  --   level      number of added slices;
  --   tol        tolerance on membership.

  -- ON RETURN :
  --   filter     interpolating polynomial.

    n : constant integer32 := Number_of_Variables(mpspt);
    prevnpts,npts : natural32;
    intpoly : Poly;
    stop : boolean;
    eva : Complex_Number;
    abseva,invcond : Floating_Number;
    maxterms : constant natural32
             := Minimum(Number_of_Terms(level+1,maxdeg),maxsamples);
    projpts : Multprec_Complex_VecVecs.VecVec(1..integer32(maxterms)+1);
    genpt : Multprec_Complex_Vectors.Vector(1..n);
    stspt : constant Standard_Sample := Original(mpspt);
    sps,sps_last,prev,tmp : Multprec_Sample_List;
    startind : natural32;

  begin
    prevnpts := 1;
    Append(sps,sps_last,mpspt);
    prev := sps;
    for d in 1..maxdeg loop                     -- for increasing degrees
      intpoly := Create(d,level+1,1);
      npts := Number_of_Terms(intpoly);         -- oversampling with 1
      Sample(stspt,npts-prevnpts,sps,sps_last);
      if prevnpts = 1
       then tmp := prev;          startind := prevnpts;
       else tmp := Tail_Of(prev); startind := prevnpts+1;
      end if;
      for j in startind..npts loop
        genpt := Sample_Point(Head_Of(tmp)).v;
        Write_Diagnostics(file,Sample_Point(Head_Of(tmp)));
        projpts(integer32(j)) := new Multprec_Complex_Vectors.Vector'
          (Evaluate(prohyp,genpt,integer32(level)+1));
        tmp := Tail_Of(tmp);
      end loop;
      Interpolate(intpoly,projpts(1..integer32(npts)-1),filter,invcond);
      put(file,"Estimate for inverse condition number : ");
      put(file,invcond,3); new_line(file);
      Clear(invcond);
      eva := Eval(filter,projpts(integer32(npts)).all);
      abseva := AbsVal(eva);
      put(file,"Residual in stop test of interpolator : ");
      put(file,abseva,3);
      stop := (abseva < tol);
      if stop
       then put(file,"  success at degree ");
       else put(file,"  failure at degree ");
      end if;
      put(file,d,1); new_line(file);
      Clear(abseva); Clear(eva);
      Clear(intpoly);
      exit when stop or (d = maxdeg); -- even if stop test fails, keep filter
      Clear(filter);
      prevnpts := npts;
      prev := sps_last;
    end loop;
    for i in 1..integer32(npts) loop
      Multprec_Complex_Vectors.Clear(projpts(i));
    end loop;
    Deep_Clear(sps);
  end Incremental_Interpolator;

  function Incremental_Interpolate
               ( file : file_type;
                 embsys : Standard_Complex_Poly_Systems.Poly_Sys;
                 orgsys : Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : Standard_Complex_Solutions.Solution_List;
                 hyp : Multprec_Complex_VecVecs.VecVec;
                 level,size : natural32 )
               return Multprec_Complex_Poly_Systems.Poly_Sys is

    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
    res : Multprec_Complex_Poly_Systems.Poly_Sys(1..len);
    timer : Timing_Widget;
    -- maxpts : natural := Minimum(Number_of_Terms(len,level+1),maxsamples);
    -- oversampling with 1 and no more than curve of degree 30
    membtol : constant double_float := 1.0E-16;
    stoptol : constant double_float := 1.0E-16;
    refsols : Multprec_Complex_Solutions.Solution_List;
    slihyp : constant Standard_Complex_VecVecs.VecVec := Slices(embsys,level);
    sps : Standard_Sample_List := Create(sols,slihyp);
    mpsps,mpsps_last,samples : Multprec_Sample_List;

  begin
    tstart(timer);
    Sampling_Machine.Initialize(embsys,orgsys,integer32(level),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner(size);
    Refine(sps,mpsps,mpsps_last);
    samples := mpsps;
    while not Is_Null(samples) loop
      Write_Diagnostics(file,Sample_Point(Head_Of(samples)));
      samples := Tail_Of(samples);
    end loop;
    refsols := Sample_Point_Lists.Sample_Points(mpsps);
    samples := mpsps;
    for i in 1..len loop                  -- sequences start at generic points
      if not On_Component(file,res,i,refsols,hyp,integer32(level),size,membtol)
       then Incremental_Interpolator
              (file,Head_Of(samples),hyp,natural32(len),level,stoptol,res(i));
      end if;
      samples := Tail_Of(samples);
    end loop;
    Sampling_Machine.Clear;
   -- Deep_Clear(sps);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Sampling + Interpolating through Generic Points");
    new_line(file);
    return Collect_Components(res);
  end Incremental_Interpolate;

  procedure Update_Base_Points
              ( file : in file_type;
                stspt : in Standard_Sample;
                base : in out Multprec_Complex_VecVecs.Link_to_VecVec;
                basecard : in out integer32; prevnpts : in integer32;
                restpts : in Multprec_Complex_VecVecs.Array_of_VecVecs;
                hyp : in Multprec_Complex_VecVecs.VecVec;
                projpts : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
                dim,level,nbpt : in integer32; size : in natural32;
                restsols : in Standard_Complex_Solutions.Solution_List ) is

  -- DESCRIPTION :
  --   Updates the list of base points by sampling one point and 
  --   computing the projections again of the restricted points;

    use Multprec_Complex_VecVecs;
    mpspt : Multprec_Sample;

  begin
    Sample(stspt,mpspt);
    Write_Diagnostics(file,Sample_Point(mpspt));
    if base = null then
      base := new Multprec_Complex_VecVecs.VecVec(1..dim-level+1);
      base(1) := new Multprec_Complex_Vectors.Vector'
                       (Sample_Point(mpspt).v);
      basecard := 1;
      for j in 1..prevnpts-1 loop
        Multprec_Complex_Vectors.Clear(projpts(nbpt)(j));              
        projpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
          (Intersect(hyp(hyp'first).all,base(1).all,
                     restpts(nbpt)(j).all,level+1));
      end loop;
    else
      basecard := basecard+1;
      base(basecard) := new Multprec_Complex_Vectors.Vector'
                             (Sample_Point(mpspt).v);
      declare
        proj : Multprec_Complex_Vectors.Vector(base(basecard)'range);
      begin
        proj := Intersect(hyp(1..basecard-1),base(1..basecard-1),
                          base(basecard).all,proj'last);
        Multprec_Complex_Vectors.Clear(base(basecard));
        base(basecard) := new Multprec_Complex_Vectors.Vector'(proj);
      end;
      for j in 1..prevnpts-1 loop
        Multprec_Complex_Vectors.Clear(projpts(nbpt)(j));              
        projpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
          (Intersect(hyp(1..basecard),base(1..basecard),
                     restpts(nbpt)(j).all,level+1));
      end loop;
    end if;
  end Update_Base_Points;

  procedure Dynamic_Incremental_Interpolator
              ( file : in file_type; mpspt : in Multprec_Sample;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                hyp : in Multprec_Complex_VecVecs.VecVec;
                nbpt,maxdeg,n,level : in integer32; size : in natural32;
                tol : in double_float;
                genpts,projpts,restpts
                   : in out Multprec_Complex_VecVecs.Array_of_VecVecs;
                centproj : in boolean;
                subspace : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                pivots : out Standard_Integer_Vectors.Link_to_Vector;
                base : out Multprec_Complex_VecVecs.Link_to_VecVec;
                basecard : in out integer32; filter : out Poly ) is

  -- DESCRIPTION :
  --   This is the internal routine in the main loop of the incremental
  --   interpolation method to break up components of solutions.

  -- ON ENTRY :
  --   file       to write diagnostics and intermediate results;
  --   mpspt      current multi-precision sample;
  --   embsys     embedded polynomial system;
  --   sols       list of generic points as solutions to embsys;
  --   hyp        general hyperplanes to project onto;
  --   nbpt       index to current generic point;
  --   maxdeg     maximal degree for the interpolating polynomial;
  --   n          original dimension;
  --   level      number of added slices;
  --   size       size of the multi-precision numbers;
  --   genpts     structure to collect samples starting at generic points;
  --   projpts    projected sample points;
  --   restpts    points restricted to the linear container space;
  --   centproj   flag to indicate whether to use central projections;
  --   basecard   counts the number of base points.

  -- ON RETURN :
  --   genpts     updated samples;
  --   projpts    updated projected points;
  --   restpts    updated restricted points;
  --   subspace   equations for the container subspace, empty when this
  --              space is the whole space;
  --   pivots     remaining variables when restricting to subspace;
  --   base       base points used in projections;
  --   basecard   updated number of base points;
  --   filter     interpolating polynomial.

    prevnpts,npts,k : integer32;
    intpoly : Poly;
    stop,restricted,skewline : boolean;
    abseva,invcond : Floating_Number;
    dim : integer32 := n+1;
    restorgsys : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    restembsys,restemb : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    restorg : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    restsols : Standard_Complex_Solutions.Solution_List;
    eva : Complex_Number;
   -- insub : boolean;
    stspt : Standard_Sample := Original(mpspt);
    sps,sps_last,prev,tmp : Multprec_Sample_List;
    startind : integer32 := 1;

    use Multprec_Complex_Vectors;
    use Multprec_Complex_VecVecs;

  begin
    prevnpts := 1;
    Append(sps,sps_last,mpspt);
    prev := sps;
    restricted := false;
    skewline := false;
    Sampling_Machine.Initialize(embsys,orgsys,level,size);
    for d in 1..maxdeg loop                     -- for increasing degrees
      if restricted then
        if centproj then
          if not skewline and then (dim > level+1) then
            skewline := true;
            Update_Base_Points(file,stspt,base,basecard,prevnpts,restpts,
                               hyp,projpts,dim,level,nbpt,size,restsols);
          else 
            if dim - basecard > level+1 then
              Update_Base_Points
                (file,stspt,base,basecard,prevnpts,restpts,
                 hyp,projpts,dim,level,nbpt,size,restsols);
            end if;
          end if;
        end if;
        if skewline
         then intpoly := Create(natural32(d-basecard),natural32(level)+1,1);
         else intpoly := Create(natural32(d),natural32(level)+1,1);
        end if;
        npts := integer32(Number_of_Terms(intpoly));  -- oversampling with 1
        if prevnpts <= npts then
          Sample(stspt,natural32(npts-prevnpts),sps,sps_last);
          tmp := Tail_Of(prev);          -- prevnpts > 1 anyway
          for j in prevnpts+1..npts loop
            restpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
              (Sample_Point(Head_Of(tmp)).v);
            Write_Diagnostics(file,Sample_Point(Head_Of(tmp)));
            tmp := Tail_Of(tmp);
          end loop;
        end if;
        if skewline then
          for j in startind..npts loop
            projpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
              (Intersect(hyp(1..basecard),base(1..basecard),
                         restpts(nbpt)(j).all,level+1));
          end loop;
        else
          for j in startind..npts loop
            projpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
              (Evaluate(hyp,restpts(nbpt)(j).all,level+1));
          end loop;
        end if;
      else 
        intpoly := Create(natural32(d),natural32(level)+1,1); 
        npts := integer32(Number_of_Terms(intpoly));  -- oversampling with 1
        Sample(stspt,natural32(npts-prevnpts),sps,sps_last);
        if prevnpts = 1
         then tmp := prev;          startind := prevnpts;
         else tmp := Tail_Of(prev); startind := prevnpts+1;
        end if;
        for j in startind..npts loop
          genpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
            (Sample_Point(Head_Of(tmp)).v);
          Write_Diagnostics(file,Sample_Point(Head_Of(tmp)));
          tmp := Tail_Of(tmp);
        end loop;
        k := npts-1;
        if k > n
         then k := n;
        end if;
        Subspace_Restriction
          (file,orgsys,embsys,natural32(nbpt),natural32(k),natural32(n),
           natural32(level),size,genpts,restpts,pivots,subspace,
           restorgsys,restembsys,natural32(dim));
        put(file,"Number of vectors : ");
        put(file,k,1); 
        put(file,"  Dimension of subspace : ");
        put(file,dim,1); new_line(file);
         if dim < k then
           restricted := true;
          -- for i in 1..npts loop  -- test on correctness
          --   insub := In_Subspace
          --              (file,subspace.all,nbpt,genpts(nbpt)(i).all,tol);
          -- end loop;
           restorg := new Multprec_Complex_Poly_Systems.Poly_Sys'
             (Collapse_Equations(restorgsys.all,natural32(dim)));
          -- put_line(file,"The collapsed original equations :");
          -- put_line(file,restorg.all);
          -- restemb := new Standard_Complex_Poly_Systems.Poly_Sys'
          --   (Collapse_Equations(restembsys.all,dim,level));
           restemb := new Standard_Complex_Poly_Systems.Poly_Sys'
             (Embed_Collapsed_Equations
                (restorg.all,restembsys.all,natural32(dim),natural32(level)));
           Sampling_Machine.Clear;
           Sampling_Machine.Initialize(restemb.all,restorg.all,level,size);
          -- put_line(file,"The collapsed embedded equations :");
          -- put_line(file,restemb.all);
          -- put(file,"The pivots ");
          -- put(file,pivots.all); new_line(file);
           restsols := Restrict_Solution(sols,nbpt,level,pivots.all);
           stspt := Create
              (Standard_Complex_Solutions.Retrieve
                 (restsols,natural32(nbpt)).all,Hyperplane_Sections(stspt));
          -- put_line(file,"The restricted solution : ");
          -- put(file,Standard_Complex_Solutions.Get(restsols,nbpt));
          -- new_line(file);
          -- declare
          --   s : Standard_Complex_Solutions.Solution(restemb'last)
          --     := Standard_Complex_Solutions.Get(restsols,nbpt);
          -- begin
          --   put_line(file,"Evaluation in restricted equations : ");
          --   put_line(file,Eval(restembsys.all,s.v));
          --   put_line(file,"Evaluation in collapsed equations : ");
          --   put_line(file,Eval(restemb.all,s.v));
          -- end;
           for j in 1..prevnpts-1 loop   -- recompute projections
             Multprec_Complex_Vectors.Clear(projpts(nbpt)(j));
             projpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
               (Evaluate(hyp,restpts(nbpt)(j).all,level+1));
           end loop;
           for j in prevnpts..npts loop
             projpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
               (Evaluate(hyp,restpts(nbpt)(j).all,level+1));
           end loop;
         else
           for j in startind..npts loop
             projpts(nbpt)(j) := new Multprec_Complex_Vectors.Vector'
               (Evaluate(hyp,genpts(nbpt)(j).all,level+1));
           end loop;
         end if;
      end if;
      if restricted and then (dim = 1) then
        Copy(intpoly,filter);  -- there is no need to do much else...
        stop := true;
      else
        Interpolate(intpoly,projpts(nbpt)(1..npts-1),filter,invcond);
        put(file,"Estimate for inverse condition number : ");
        put(file,invcond,3); new_line(file);
        Clear(invcond);
        declare
          extpt : Multprec_Complex_Vectors.Vector(1..level+1);
        begin
          Multprec_Complex_Vectors.Copy(projpts(nbpt)(npts).all,extpt);
          Set_Size(extpt,2*size);
          -- eva := Eval(filter,projpts(nbpt)(npts).all);
          eva := Eval(filter,extpt);
          Multprec_Complex_Vectors.Clear(extpt);
        end;
        abseva := AbsVal(eva);
        put(file,"Residual in stop test of interpolator : ");
        put(file,abseva,3);
        stop := (abseva < tol);
        if stop
         then put(file,"  success at degree ");
         else put(file,"  failure at degree ");
        end if;
        put(file,d,1); new_line(file);
        Clear(abseva);
        Clear(eva);
      end if;
      Clear(intpoly);
      exit when stop or (d = maxdeg);  -- even if stop test fails, keep filter
      Clear(filter);
      prevnpts := npts;
      prev := sps_last;
    end loop;
    Deep_Clear(sps);
    Sampling_Machine.Clear;
  end Dynamic_Incremental_Interpolator;

  procedure Dynamic_Interpolate
              ( file : in file_type;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List; 
                level,size : in natural32;
                hyp : Multprec_Complex_VecVecs.VecVec;
                centproj : in boolean;
                subspaces : out Multprec_Complex_Poly_Systems.Array_of_Poly_Sys;
                pivots : out Standard_Integer_VecVecs.VecVec;
                basepts : out Multprec_Complex_VecVecs.Array_of_VecVecs;
                basecard : out Standard_Natural_Vectors.Vector;
                filters
                  : out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys ) is

    n : constant integer32 := orgsys'last;
    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
      -- #generic points
    res : Multprec_Complex_Poly_Systems.Poly_Sys(1..len);
    timer : Timing_Widget;
    genpts,projpts,restpts
      : Multprec_Complex_VecVecs.Array_of_VecVecs(1..len);
    maxpts : constant natural32
           := Minimum(Number_of_Terms(natural32(len),level+1),maxsamples);
      -- oversampling with 1 and no more than degree 30 for curves
    membtol : constant double_float := 1.0E-16;
    stoptol : constant double_float := 1.0E-16;
    refsols : Multprec_Complex_Solutions.Solution_List;
    slihyp : constant Standard_Complex_VecVecs.VecVec := Slices(embsys,level);
    sps : Standard_Sample_List := Create(sols,slihyp);
    mpsps,mpsps_last,samples : Multprec_Sample_List;

  begin
    tstart(timer);
    subspaces(subspaces'first) := null;  -- only to avoid warnings...
    pivots(pivots'first) := null;
    basepts(basepts'first) := null;
    genpts := Create(len,integer32(maxpts));
    projpts := Create(len,integer32(maxpts));
    restpts := Create(len,integer32(maxpts));
    Sampling_Machine.Initialize(embsys,orgsys,integer32(level),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner(size);
    Refine(sps,mpsps,mpsps_last);
    samples := mpsps;
    while not Is_Null(samples) loop
      Write_Diagnostics(file,Sample_Point(Head_Of(samples)));
      samples := Tail_Of(samples);
    end loop;
    Sampling_Machine.Clear;
    refsols := Sample_Point_Lists.Sample_Points(mpsps);
    samples := mpsps;
    for i in 1..len loop                  -- sequences start at generic points
      basecard(i) := 0;
      if not On_Component(file,res,i,subspaces,pivots,basepts,basecard,
                          refsols,hyp,integer32(level),size,membtol)
       then Dynamic_Incremental_Interpolator
              (file,Head_Of(samples),embsys,orgsys,sols,hyp,i,len,n,
               integer32(level),size,stoptol,genpts,
               projpts,restpts,centproj,subspaces(i),pivots(i),basepts(i),
               integer32(basecard(i)),res(i));
      end if;
      samples := Tail_Of(samples);
    end loop;
   -- Multprec_Complex_VecVecs.Clear(genpts);
   -- Multprec_Complex_VecVecs.Clear(projpts);
   -- Multprec_Complex_VecVecs.Clear(restpts);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Sampling + Interpolating through Generic Points");
    new_line(file);
  -- WARNING : THIS COULD LEAD TO INCONSISTENCIES WITH SUBSPACES !
  --  filters := new Multprec_Complex_Poly_Systems.Poly_Sys'
  --                   (Collect_Components(res));
    filters := new Multprec_Complex_Poly_Systems.Poly_Sys'(res);
  end Dynamic_Interpolate;

end Multprec_Breakup_Components;
