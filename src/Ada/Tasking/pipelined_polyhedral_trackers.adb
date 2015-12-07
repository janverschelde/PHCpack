with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;

with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Lists_of_Floating_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Floating_Vector_Lists;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_JacoMats;
with Exponent_Vectors;
with Supports_of_Polynomial_Systems;
with Floating_Mixed_Subdivisions_io;
with Random_Coefficient_Systems;
with MixedVol_Algorithm;                use MixedVol_Algorithm;
with Single_Polyhedral_Trackers;        use Single_Polyhedral_Trackers;
with Pipelined_Labeled_Cells;           use Pipelined_Labeled_Cells;

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
    put("mix = "); put(mix); new_line;
    put("idx = "); put(idx.all); new_line;
    for k in 1..idx'last loop
      ind := ind + 1;
      put("support "); put(ind,1); put_line(" :");
      for i in idx(k-1)..(idx(k)-1) loop
        vpt := vtx(i);
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

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,r : in integer32;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Complex_Solutions;

    mix : constant Standard_Integer_Vectors.Vector := Mixture(r,mtype);
    lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range)
        := Lifted_Supports(nbequ,r,mix,idx,vtx,lft);
    sup : Arrays_of_Integer_Vector_Lists.Array_of_Lists(1..nbequ);
    pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists(1..nbequ);
    hom : Eval_Coeff_Laur_Sys(1..nbequ);
    cff : Standard_Complex_VecVecs.VecVec(hom'range);
    dpw : Standard_Floating_VecVecs.Link_to_VecVec;
    cft : Standard_Complex_VecVecs.Link_to_VecVec;
    epv : Exponent_Vectors.Exponent_Vectors_Array(hom'range);
    ejf : Eval_Coeff_Jaco_Mat(hom'range,hom'first..hom'last+1);
    jmf : Mult_Factors(ejf'range(1),ejf'range(2));
    tasksols,lastsols : Array_of_Solution_Lists(2..nt);

  begin
    q := Random_Coefficient_Systems.Create(natural32(nbequ),mix,lif);
    put_line(file,q);
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    Floating_Mixed_Subdivisions_io.put(file,lif);
    for i in cff'range loop
      declare
        coeff_lq : constant Standard_Complex_Vectors.Vector
                 := Standard_Complex_Laur_Functions.Coeff(q(i));
      begin
        cff(i) := new Standard_Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          cff(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    epv := Exponent_Vectors.Create(q);
    Create(q,ejf,jmf);
    declare
      lif : Arrays_of_Floating_Vector_Lists.Array_of_Lists(mix'range);

      procedure Track
        ( idtask,r : in integer32; 
          mtype : in Standard_Integer_Vectors.Link_to_Vector;
          mic : in out Mixed_Cell ) is

        sol : Standard_Complex_Solutions.Link_to_Solution;

      begin
        null;
       -- Track_Path(mix,lif,mic.nor,cff,dpw,cft,epv,hom,ejf,jmf,sol);
       -- Append(tasksols(idtask),lastsols(idtask),sol);
      end Track;

    begin
      Pipelined_Mixed_Cells
        (nt,nbequ,true,r,mtype,perm,idx,vtx,lft,mcc,mv,Track'access);
    end;
  end Reporting_Multitasking_Tracker;

  procedure Reporting_Multitasking_Tracker
              ( file : in file_type;
                nt,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                mcc : out Mixed_Subdivision; mv : out natural32;
                q : out Standard_Complex_Laur_Systems.Laur_Sys;
                sols : out Standard_Complex_Solutions.Solution_List ) is

    stlb : constant double_float := 0.0; -- no stable mv for now...
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    put("ind = "); put(ind); new_line;
    put("cnt = "); put(cnt); new_line;
    put("perm = "); put(perm); new_line;
    put("idx = "); put(idx); new_line;
    put("sdx = "); put(sdx); new_line;
    put("ndx = "); put(ndx); new_line;
    mv_lift(nbequ,nbpts,ind,cnt,support.all,stlb,r,idx,vtx,lft);
    Reporting_Multitasking_Tracker
      (file,nt,nbequ,r,mtype,perm,idx,vtx,lft,mcc,mv,q,sols);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Reporting_Multitasking_Tracker;

end Pipelined_Polyhedral_Trackers;
