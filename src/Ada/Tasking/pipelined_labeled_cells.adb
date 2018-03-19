with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Mixed_Volume;
with MixedVol_Algorithm;                use MixedVol_Algorithm;
with Mixed_Labels_Queue;
with Multitasking;

package body Pipelined_Labeled_Cells is

  function Mixture ( r : integer32;
                     mtype : in Standard_Integer_Vectors.Link_to_Vector )
                   return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..r);

  begin
    for i in 1..r loop
      res(i) := mtype(i-1);
    end loop;
    return res;
  end Mixture;

  procedure Produce_Cells
              ( nbequ,r : in integer32; otp : in boolean;
                mtype,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                cells : out CellStack; nbcells : out integer32;
                mixvol : out natural32 ) is

    cellcnt : natural32 := 0;
    CellSize : constant integer32 := Mixed_Volume.cell_size(r,mtype);

    procedure append ( pts : Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      cellcnt := cellcnt + 1;
      if otp then
        put_line("Appending cell " & Multitasking.to_string(cellcnt)
                                   & " to the queue");
      end if;
      Mixed_Labels_Queue.Append(pts);
    end append;

  begin
    nbcells := 0;
    mixvol := 0;
    if otp
     then put_line("starting the cell production ...");
    end if;
    Cs_Init(cells,CellSize);
    Mixed_Volume.MixedVol_with_Callback
      (nbequ,r,CellSize,mtype,idx,vtx,lft,nbcells,cells,mixvol,false,
       append'access);
    Mixed_Labels_Queue.Stop;
    if otp then
      put_line("The mixed volume is " & Multitasking.to_string(mixvol) & ".");
      put("There are " & Multitasking.to_string(nbcells) & " mixed cells.");
    end if;
  end Produce_Cells;

  procedure Process_Cells 
              ( idtask,nbequ,r : in integer32; otp : in boolean;
                mtype,perm : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                mcc : in out Mixed_Subdivision;
                process : access procedure
                  ( idtask,r : in integer32;
                    mtype : in Standard_Integer_Vectors.Link_to_Vector;
                    mic : in out Mixed_Cell ) := null ) is

    labels : Standard_Integer_Vectors.Link_to_Vector;
    mcc_last : Mixed_Subdivision;
    stop : boolean := false;
    cnt : natural32 := 0;

    use Standard_Integer_Vectors;

  begin
    while not stop loop
      if otp
       then Mixed_Labels_Queue.Next(labels,cnt);
       else labels := Mixed_Labels_Queue.Next;
      end if;
      if labels = null then -- check if all labels are produced
        if Mixed_Labels_Queue.Stopped     -- production stopped
         then stop := true;   -- all labels have been processed
        end if;
      else
        if otp then
          put_line("Task " & Multitasking.to_string(idtask)
                           & " processes cell "
                           & Multitasking.to_string(cnt));
        end if;
        if r < nbequ then
          declare
            mic : Mixed_Cell
                := Labels_to_Mixed_Cell(nbequ,r,mtype,labels,Vtx,lft);
          begin
            Append(mcc,mcc_last,mic);
            if process /= null
             then process(idtask,r,mtype,mic);
            end if;
          end;
        else
          declare
            mic : Mixed_Cell
                := Labels_to_Mixed_Cell(nbequ,r,mtype,perm,labels,Vtx,lft);
          begin
            Append(mcc,mcc_last,mic);
            if process /= null
             then process(idtask,r,mtype,mic);
            end if;
          end;
        end if;
      end if;
    end loop;
    if otp then
      put_line("Task " & Multitasking.to_string(idtask) & " processed "
                       & Multitasking.to_string(Length_Of(mcc))
                       & " mixed cells.");
    end if;
  end Process_Cells;

  procedure Pipelined_Mixed_Cells
              ( ntasks,nbequ : in integer32; otp : in boolean;
                r : in integer32;
                mtype,perm,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                sub : out Mixed_Subdivision; mv : out natural32;
                process : access procedure
                  ( idtask,r : in integer32;
                    mtype : in Standard_Integer_Vectors.Link_to_Vector;
                    mic : in out Mixed_Cell ) := null ) is

    mcc : array(2..ntasks) of Mixed_Subdivision;
    cells : CellStack;
    nbc : integer32;
    sub_last : Mixed_Subdivision;

    procedure do_job ( i,n : in integer32 ) is
    begin
      if otp
       then put_line("In do_job with task " & Multitasking.to_string(i));
      end if;
      if i = 1 then
        Produce_Cells(nbequ,r,otp,mtype,idx,vtx,lft,cells,nbc,mv);
      elsif process /= null then
        Process_Cells(i,nbequ,r,otp,mtype,perm,vtx,lft,mcc(i),process);
      else
        Process_Cells(i,nbequ,r,otp,mtype,perm,vtx,lft,mcc(i));
      end if;
    end do_job;
    procedure rep_do_jobs is new Multitasking.Reporting_Workers(do_job);
    procedure sil_do_jobs is new Multitasking.Silent_Workers(do_job);

  begin
    Mixed_Labels_Queue.Start;
    if otp
     then rep_do_jobs(ntasks);
     else sil_do_jobs(ntasks);
    end if;
    for i in mcc'range loop
      Concat(sub,sub_last,mcc(i));
    end loop;
    Cs_Del(cells);
  end Pipelined_Mixed_Cells;

  procedure Pipelined_Mixed_Cells
              ( ntasks,nbequ,nbpts : in integer32; otp : in boolean;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                sub : out Mixed_Subdivision; mv : out natural32;
                process : access procedure
                  ( idtask,r : in integer32;
                    mtype : in Standard_Integer_Vectors.Link_to_Vector;
                    mic : in out Mixed_Cell ) := null ) is

    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,0.0,r,idx,vtx,lft);
    Pipelined_Mixed_Cells
      (ntasks,nbequ,otp,r,mtype,perm,idx,vtx,lft,sub,mv,process);
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
  end Pipelined_Mixed_Cells;

end Pipelined_Labeled_Cells;
