with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Lists_of_Floating_Vectors;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_JacoMats;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with Exponent_Vectors;
with Floating_Mixed_Subdivisions_io;
with Induced_Permutations;
with Random_Coefficient_Systems;
with MixedVol_Algorithm;                use MixedVol_Algorithm;
with Semaphore;
with Pipelined_Labeled_Cells;           use Pipelined_Labeled_Cells;
with Polyhedral_Start_Systems;          use Polyhedral_Start_Systems;
with Pipelined_Cell_Trackers;           use Pipelined_Cell_Trackers;

package body Pipelined_Polyhedral_Trackers is

  function Lifted_Supports
              ( n,r : integer32;
                mix : Standard_Integer_Vectors.Vector;
                idx : Standard_Integer_Vectors.Link_to_Vector;
                vtx : Standard_Integer_VecVecs.Link_to_VecVec;
                lft : Standard_Floating_Vectors.Link_to_Vector )
              return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

  -- DESCRIPTION :
  --   Joins the vertex points and their lifting values into one
  --   array of lists of lifted supports.

  -- ON ENTRY :
  --   n        ambient dimension of the points, before the lifting;
  --   r        number of different supports;
  --   mix      type of mixture, number of occurrences of each support;
  --   idx      indices to the vertex points;
  --   vtx      coordinates of the vertex points;
  --   lft      lifting values for the vertex points.

  -- ON RETURN :
  --   An array of range 1..r of lifted supports.

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    res_last : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..r);
    ind : integer32 := 0;
    vpt : Standard_Integer_Vectors.Link_to_Vector;
    idxlft : integer32 := lft'first-1;

  begin
   -- put("mix = "); put(mix); new_line;
   -- put("idx = "); put(idx.all); new_line;
    for k in 1..r loop
      ind := ind + 1;
     -- put("support "); put(ind,1); put_line(" :");
      for i in idx(k-1)..(idx(k)-1) loop
        vpt := vtx(i);
       -- put(vpt); new_line;
        declare
          lpt : Standard_Floating_Vectors.Vector(1..n+1);
          ilp : integer32 := 0;
        begin
          for j in vpt'range loop
            ilp := ilp + 1;
            lpt(ilp) := double_float(vpt(j));
          end loop;
          idxlft := idxlft + 1;
          lpt(n+1) := lft(idxlft);
          Lists_of_Floating_Vectors.Append(res(ind),res_last(ind),lpt);
        end;
      end loop;
    end loop;
    return res;
  end Lifted_Supports;

  function Stable_Start_System
             ( nbequ : integer32; stlb : double_float;
               mix : Standard_Integer_Vectors.Vector;
               lifsup : Array_of_Lists )
             return Standard_Complex_Laur_Systems.Laur_Sys is

  -- DESCRIPTION :
  --   For stable polyhedral continuation, artificial origins were
  --   added to the lifted supports.
  --   The artificial origins with lifting bound equal to stlb
  --   are removed in the random coefficient start system on return.

    sup : Array_of_Lists(lifsup'range)
        := Induced_Permutations.Remove_Artificial_Origin(lifsup,stlb);
    res : Standard_Complex_Laur_Systems.Laur_Sys(1..nbequ)
        := Random_Coefficient_Systems.Create(natural32(nbequ),mix,sup);

  begin
    for i in sup'range loop
      Lists_of_Floating_Vectors.Clear(sup(i));
    end loop;
    return res;
  end Stable_Start_System;

  function Stable_Start_System
             ( nbequ : integer32; stlb : double_float;
               mix : Standard_Integer_Vectors.Vector;
               lifsup : Array_of_Lists )
             return DoblDobl_Complex_Laur_Systems.Laur_Sys is

  -- DESCRIPTION :
  --   For stable polyhedral continuation, artificial origins were
  --   added to the lifted supports.
  --   The artificial origins with lifting bound equal to stlb
  --   are removed in the random coefficient start system on return.

    sup : Array_of_Lists(lifsup'range)
        := Induced_Permutations.Remove_Artificial_Origin(lifsup,stlb);
    res : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..nbequ)
        := Random_Coefficient_Systems.Create(natural32(nbequ),mix,sup);

  begin
    for i in sup'range loop
      Lists_of_Floating_Vectors.Clear(sup(i));
    end loop;
    return res;
  end Stable_Start_System;

  function Stable_Start_System
             ( nbequ : integer32; stlb : double_float;
               mix : Standard_Integer_Vectors.Vector;
               lifsup : Array_of_Lists )
             return QuadDobl_Complex_Laur_Systems.Laur_Sys is

  -- DESCRIPTION :
  --   For stable polyhedral continuation, artificial origins were
  --   added to the lifted supports.
  --   The artificial origins with lifting bound equal to stlb
  --   are removed in the random coefficient start system on return.

    sup : Array_of_Lists(lifsup'range)
        := Induced_Permutations.Remove_Artificial_Origin(lifsup,stlb);
    res : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..nbequ)
        := Random_Coefficient_Systems.Create(natural32(nbequ),mix,sup);

  begin
    for i in sup'range loop
      Lists_of_Floating_Vectors.Clear(sup(i));
    end loop;
    return res;
  end Stable_Start_System;

-- AFTER PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, SILENT :

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    sem : Semaphore.Lock;
    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
    permlif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
            := Lifted_Supports(nbequ,r,mix,idx,vtx,lft);
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : Standard_Complex_VecVecs.VecVec(hom'range);
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tmv : Standard_Natural_Vectors.Vector(2..nt) := (2..nt => 0);
    dpw : Standard_Floating_VecVecs.Array_of_VecVecs(2..nt);
    cft : Standard_Complex_VecVecs.Array_of_VecVecs(2..nt);
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);
    permq : Standard_Complex_Laur_Systems.Laur_Sys(1..nbequ);

    procedure Track ( idtask,r : in integer32; 
                      mtype : in Standard_Integer_Vectors.Link_to_Vector;
                      mic : in out Mixed_Cell ) is
    begin
      if( (stlb = 0.0) or else
         ((stlb > 0.0) and then Is_Original(mic,stlb)) ) then
        Standard_Track_Cell(sem,idtask,nbequ,r,mix,mic,lifsup,cff,
          dpw(idtask),cft(idtask),epv,hom,ejf,jmf,q,tmv(idtask),
          tasksols(idtask),lastsols(idtask));
      end if;
    end Track;

  begin
    if r < nbequ then
      if stlb = 0.0 then
        q := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        q := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      lifsup := permlif;
    else
      if stlb = 0.0 then
        permq
          := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        permq := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      for i in perm'range loop
        q(perm(i)+1) := permq(i+1);
        lifsup(perm(i)+1) := permlif(i+1);
      end loop;
    end if;
    lif := new Array_of_Lists'(lifsup);
    cff := Coeff(q);
    epv := Exponent_Vectors.Create(q);
    hom := Create(q);
    Create(q,ejf,jmf);
    Allocate_Workspace_for_Exponents(epv,dpw);
    Allocate_Workspace_for_Coefficients(cff,cft);
    Pipelined_Mixed_Cells
      (nt,nbequ,false,r,mtype,perm,idx,vtx,lft,mcc,mv,Track'access);
    for k in tasksols'range loop
      Standard_Complex_Solutions.Push(tasksols(k),sols);
    end loop;
  end Silent_Multitasking_Tracker;

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;

    sem : Semaphore.Lock;
    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
    permlif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
            := Lifted_Supports(nbequ,r,mix,idx,vtx,lft);
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : DoblDobl_Complex_VecVecs.VecVec(hom'range);
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tmv : Standard_Natural_Vectors.Vector(2..nt) := (2..nt => 0);
    dpw : Standard_Floating_VecVecs.Array_of_VecVecs(2..nt);
    cft : DoblDobl_Complex_VecVecs.Array_of_VecVecs(2..nt);
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);
    permq : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..nbequ);

    procedure Track ( idtask,r : in integer32; 
                      mtype : in Standard_Integer_Vectors.Link_to_Vector;
                      mic : in out Mixed_Cell ) is
    begin
      if( (stlb = 0.0) or else
         ((stlb > 0.0) and then Is_Original(mic,stlb)) ) then
        DoblDobl_Track_Cell(sem,idtask,nbequ,r,mix,mic,lifsup,cff,
          dpw(idtask),cft(idtask),epv,hom,ejf,jmf,q,tmv(idtask),
          tasksols(idtask),lastsols(idtask));
      end if;
    end Track;

  begin
    if r < nbequ then
      if stlb = 0.0 then
        q := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        q := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      lifsup := permlif;
    else
      if stlb = 0.0 then
        permq
          := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        permq := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      for i in perm'range loop
        q(perm(i)+1) := permq(i+1);
        lifsup(perm(i)+1) := permlif(i+1);
      end loop;
    end if;
    lif := new Array_of_Lists'(lifsup);
    cff := Coeff(q);
    epv := Exponent_Vectors.Create(q);
    hom := Create(q);
    Create(q,ejf,jmf);
    Allocate_Workspace_for_Exponents(epv,dpw);
    Allocate_Workspace_for_Coefficients(cff,cft);
    Pipelined_Mixed_Cells
      (nt,nbequ,false,r,mtype,perm,idx,vtx,lft,mcc,mv,Track'access);
    for k in tasksols'range loop
      DoblDobl_Complex_Solutions.Push(tasksols(k),sols);
    end loop;
  end Silent_Multitasking_Tracker;

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;

    sem : Semaphore.Lock;
    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
    permlif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
            := Lifted_Supports(nbequ,r,mix,idx,vtx,lft);
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : QuadDobl_Complex_VecVecs.VecVec(hom'range);
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tmv : Standard_Natural_Vectors.Vector(2..nt) := (2..nt => 0);
    dpw : Standard_Floating_VecVecs.Array_of_VecVecs(2..nt);
    cft : QuadDobl_Complex_VecVecs.Array_of_VecVecs(2..nt);
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);
    permq : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..nbequ);

    procedure Track ( idtask,r : in integer32; 
                      mtype : in Standard_Integer_Vectors.Link_to_Vector;
                      mic : in out Mixed_Cell ) is
    begin
      if( (stlb = 0.0) or else
         ((stlb > 0.0) and then Is_Original(mic,stlb)) ) then
        QuadDobl_Track_Cell(sem,idtask,nbequ,r,mix,mic,lifsup,cff,
          dpw(idtask),cft(idtask),epv,hom,ejf,jmf,q,tmv(idtask),
          tasksols(idtask),lastsols(idtask));
      end if;
    end Track;

  begin
    if r < nbequ then
      if stlb = 0.0 then
        q := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        q := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      lifsup := permlif;
    else
      if stlb = 0.0 then
        permq
          := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        permq := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      for i in perm'range loop
        q(perm(i)+1) := permq(i+1);
        lifsup(perm(i)+1) := permlif(i+1);
      end loop;
    end if;
    lif := new Array_of_Lists'(lifsup);
    cff := Coeff(q);
    epv := Exponent_Vectors.Create(q);
    hom := Create(q);
    Create(q,ejf,jmf);
    Allocate_Workspace_for_Exponents(epv,dpw);
    Allocate_Workspace_for_Coefficients(cff,cft);
    Pipelined_Mixed_Cells
      (nt,nbequ,false,r,mtype,perm,idx,vtx,lft,mcc,mv,Track'access);
    for k in tasksols'range loop
      QuadDobl_Complex_Solutions.Push(tasksols(k),sols);
    end loop;
  end Silent_Multitasking_Tracker;

-- AFTER PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, REPORTING :

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    sem : Semaphore.Lock;
    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
    permlif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
            := Lifted_Supports(nbequ,r,mix,idx,vtx,lft);
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : Standard_Complex_VecVecs.VecVec(hom'range);
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tmv : Standard_Natural_Vectors.Vector(2..nt) := (2..nt => 0);
    dpw : Standard_Floating_VecVecs.Array_of_VecVecs(2..nt);
    cft : Standard_Complex_VecVecs.Array_of_VecVecs(2..nt);
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);
    permq : Standard_Complex_Laur_Systems.Laur_Sys(1..nbequ);

    procedure Track ( idtask,r : in integer32; 
                      mtype : in Standard_Integer_Vectors.Link_to_Vector;
                      mic : in out Mixed_Cell ) is
    begin
      if( (stlb = 0.0) or else
         ((stlb > 0.0) and then Is_Original(mic,stlb)) ) then
        Standard_Track_Cell(sem,idtask,nbequ,r,mix,mic,lifsup,cff,
          dpw(idtask),cft(idtask),epv,hom,ejf,jmf,q,tmv(idtask),
          tasksols(idtask),lastsols(idtask));
      end if;
    end Track;

  begin
    if r < nbequ then
      if stlb = 0.0 then
        q := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        q := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      lifsup := permlif;
    else
      if stlb = 0.0 then
        permq
          := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        permq := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      for i in perm'range loop
        q(perm(i)+1) := permq(i+1);
        lifsup(perm(i)+1) := permlif(i+1);
      end loop;
    end if;
    put_line(file,q);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lifsup);
    lif := new Array_of_Lists'(lifsup);
    cff := Coeff(q);
    epv := Exponent_Vectors.Create(q);
    hom := Create(q);
    Create(q,ejf,jmf);
    Allocate_Workspace_for_Exponents(epv,dpw);
    Allocate_Workspace_for_Coefficients(cff,cft);
    Pipelined_Mixed_Cells
      (nt,nbequ,true,r,mtype,perm,idx,vtx,lft,mcc,mv,Track'access);
    put("Number of paths tracked : "); put(tmv); 
    put(" => "); put(Standard_Natural_Vectors.Sum(tmv),1); new_line;
    for k in tasksols'range loop
      Standard_Complex_Solutions.Push(tasksols(k),sols);
    end loop;
  end Reporting_Multitasking_Tracker;

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Complex_Solutions;

    sem : Semaphore.Lock;
    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
    permlif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
            := Lifted_Supports(nbequ,r,mix,idx,vtx,lft);
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : DoblDobl_Complex_VecVecs.VecVec(hom'range);
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tmv : Standard_Natural_Vectors.Vector(2..nt) := (2..nt => 0);
    dpw : Standard_Floating_VecVecs.Array_of_VecVecs(2..nt);
    cft : DoblDobl_Complex_VecVecs.Array_of_VecVecs(2..nt);
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);
    permq : DoblDobl_Complex_Laur_Systems.Laur_Sys(1..nbequ);

    procedure Track ( idtask,r : in integer32; 
                      mtype : in Standard_Integer_Vectors.Link_to_Vector;
                      mic : in out Mixed_Cell ) is
    begin
      if( (stlb = 0.0) or else
         ((stlb > 0.0) and then Is_Original(mic,stlb)) ) then
        DoblDobl_Track_Cell(sem,idtask,nbequ,r,mix,mic,lifsup,cff,
          dpw(idtask),cft(idtask),epv,hom,ejf,jmf,q,tmv(idtask),
          tasksols(idtask),lastsols(idtask));
      end if;
    end Track;

  begin
    if r < nbequ then
      if stlb = 0.0 then
        q := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        q := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      lifsup := permlif;
    else
      if stlb = 0.0 then
        permq
          := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        permq := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      for i in perm'range loop
        q(perm(i)+1) := permq(i+1);
        lifsup(perm(i)+1) := permlif(i+1);
      end loop;
    end if;
    put_line(file,q);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lifsup);
    lif := new Array_of_Lists'(lifsup);
    cff := Coeff(q);
    epv := Exponent_Vectors.Create(q);
    hom := Create(q);
    Create(q,ejf,jmf);
    Allocate_Workspace_for_Exponents(epv,dpw);
    Allocate_Workspace_for_Coefficients(cff,cft);
    Pipelined_Mixed_Cells
      (nt,nbequ,true,r,mtype,perm,idx,vtx,lft,mcc,mv,Track'access);
    put("Number of paths tracked : "); put(tmv); 
    put(" => "); put(Standard_Natural_Vectors.Sum(tmv),1); new_line;
    for k in tasksols'range loop
      DoblDobl_Complex_Solutions.Push(tasksols(k),sols);
    end loop;
  end Reporting_Multitasking_Tracker;

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,r : in integer32; stlb : in double_float;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Complex_Solutions;

    sem : Semaphore.Lock;
    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
    permlif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
            := Lifted_Supports(nbequ,r,mix,idx,vtx,lft);
    lifsup : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : QuadDobl_Complex_VecVecs.VecVec(hom'range);
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tmv : Standard_Natural_Vectors.Vector(2..nt) := (2..nt => 0);
    dpw : Standard_Floating_VecVecs.Array_of_VecVecs(2..nt);
    cft : QuadDobl_Complex_VecVecs.Array_of_VecVecs(2..nt);
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);
    permq : QuadDobl_Complex_Laur_Systems.Laur_Sys(1..nbequ);

    procedure Track ( idtask,r : in integer32; 
                      mtype : in Standard_Integer_Vectors.Link_to_Vector;
                      mic : in out Mixed_Cell ) is
    begin
      if( (stlb = 0.0) or else
         ((stlb > 0.0) and then Is_Original(mic,stlb)) ) then
        QuadDobl_Track_Cell(sem,idtask,nbequ,r,mix,mic,lifsup,cff,
          dpw(idtask),cft(idtask),epv,hom,ejf,jmf,q,tmv(idtask),
          tasksols(idtask),lastsols(idtask));
      end if;
    end Track;

  begin
    if r < nbequ then
      if stlb = 0.0 then
        q := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        q := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      lifsup := permlif;
    else
      if stlb = 0.0 then
        permq
          := Random_Coefficient_Systems.Create(natural32(nbequ),mix,permlif);
      else
        q := Stable_Start_System(nbequ,stlb,mix,permlif);
      end if;
      for i in perm'range loop
        q(perm(i)+1) := permq(i+1);
        lifsup(perm(i)+1) := permlif(i+1);
      end loop;
    end if;
    put_line(file,q);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lifsup);
    lif := new Array_of_Lists'(lifsup);
    cff := Coeff(q);
    epv := Exponent_Vectors.Create(q);
    hom := Create(q);
    Create(q,ejf,jmf);
    Allocate_Workspace_for_Exponents(epv,dpw);
    Allocate_Workspace_for_Coefficients(cff,cft);
    Pipelined_Mixed_Cells
      (nt,nbequ,true,r,mtype,perm,idx,vtx,lft,mcc,mv,Track'access);
    put("Number of paths tracked : "); put(tmv); 
    put(" => "); put(Standard_Natural_Vectors.Sum(tmv),1); new_line;
    for k in tasksols'range loop
      QuadDobl_Complex_Solutions.Push(tasksols(k),sols);
    end loop;
  end Reporting_Multitasking_Tracker;

-- DOES PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, SILENT :

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    stlb : constant double_float := 0.0;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,stlb,r,idx,vtx,lft);
    Silent_Multitasking_Tracker
      (nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Silent_Multitasking_Tracker;

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    stlb : constant double_float := 0.0;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,stlb,r,idx,vtx,lft);
    Silent_Multitasking_Tracker
      (nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Silent_Multitasking_Tracker;

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    stlb : constant double_float := 0.0;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,stlb,r,idx,vtx,lft);
    Silent_Multitasking_Tracker
      (nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Silent_Multitasking_Tracker;

-- DOES PREPROCESSING AND LIFTING, ON LAURENT SYSTEMS, REPORTING :

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    stlb : constant double_float := 0.0;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,stlb,r,idx,vtx,lft);
    Reporting_Multitasking_Tracker
      (file,nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Reporting_Multitasking_Tracker;

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    stlb : constant double_float := 0.0;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,stlb,r,idx,vtx,lft);
    Reporting_Multitasking_Tracker
      (file,nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Reporting_Multitasking_Tracker;

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    stlb : constant double_float := 0.0;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,stlb,r,idx,vtx,lft);
    Reporting_Multitasking_Tracker
      (file,nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Reporting_Multitasking_Tracker;

-- DOES PREPROCESSING AND LIFTING, ON STABLE MIXED VOLUMES, SILENT :

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    if stable
     then mv_lift(nbequ,stlb,r,idx,vtx,lft);
     else mv_lift(nbequ,0.0,r,idx,vtx,lft);
    end if;
    Silent_Multitasking_Tracker
      (nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Silent_Multitasking_Tracker;

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    if stable
     then mv_lift(nbequ,stlb,r,idx,vtx,lft);
     else mv_lift(nbequ,0.0,r,idx,vtx,lft);
    end if;
    Silent_Multitasking_Tracker
      (nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Silent_Multitasking_Tracker;

  procedure Silent_Multitasking_Tracker
              ( nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    if stable
     then mv_lift(nbequ,stlb,r,idx,vtx,lft);
     else mv_lift(nbequ,0.0,r,idx,vtx,lft);
    end if;
    Silent_Multitasking_Tracker
      (nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Silent_Multitasking_Tracker;

-- DOES PREPROCESSING AND LIFTING, STABLE MIXED VOLUMES, REPORTING :

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    if stable
     then mv_lift(nbequ,stlb,r,idx,vtx,lft);
     else mv_lift(nbequ,0.0,r,idx,vtx,lft);
    end if;
    Reporting_Multitasking_Tracker
      (file,nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Reporting_Multitasking_Tracker;

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out DoblDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out DoblDobl_Complex_Solutions.Solution_List ) is

    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    if stable
     then mv_lift(nbequ,stlb,r,idx,vtx,lft);
     else mv_lift(nbequ,0.0,r,idx,vtx,lft);
    end if;
    Reporting_Multitasking_Tracker
      (file,nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Reporting_Multitasking_Tracker;

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type; nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                stable : in boolean; stlb : in double_float;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                lif : out Link_to_Array_of_Lists;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out QuadDobl_Complex_Laur_Systems.Laur_Sys;
                sols : out QuadDobl_Complex_Solutions.Solution_List ) is

    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    if stable
     then mv_lift(nbequ,stlb,r,idx,vtx,lft);
     else mv_lift(nbequ,0.0,r,idx,vtx,lft);
    end if;
    Reporting_Multitasking_Tracker
      (file,nt,nbequ,r,stlb,mtype,perm,idx,vtx,lft,lif,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Reporting_Multitasking_Tracker;

end Pipelined_Polyhedral_Trackers;
