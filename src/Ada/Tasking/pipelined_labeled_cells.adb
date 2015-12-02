with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Cell_Stack;                        use Cell_Stack;
with Mixed_Volume;
with MixedVol_Algorithm;                use MixedVol_Algorithm;
with Mixed_Labels_Queue;
with Multitasking;

package body Pipelined_Labeled_Cells is

  procedure Produce_Cells
              ( nbequ,nbpts,r : in integer32;
                mtype,idx : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector ) is

    cells : CellStack;
    cellcnt : natural32 := 0;
    CellSize : constant integer32 := Mixed_Volume.cell_size(r,mtype);
    nbcells : integer32 := 0;
    mixvol : natural32;

    procedure write_labels ( pts : Standard_Integer_Vectors.Link_to_Vector ) is
    begin
      cellcnt := cellcnt + 1;
      put_line("Appending cell " & Multitasking.to_string(cellcnt)
                                 & " to the queue");
      Mixed_Labels_Queue.Append(pts);
    end write_labels;

  begin
    put_line("starting the cell production ...");
    Cs_Init(cells,CellSize);
    Mixed_Volume.MixedVol_with_Callback
      (nbequ,r,CellSize,mtype,idx,vtx,lft,nbcells,cells,mixvol,false,
       write_labels'access);
    Mixed_Labels_Queue.Stop;
    put_line("The mixed volume is " & Multitasking.to_string(mixvol) & ".");
    put("There are " & Multitasking.to_string(nbcells) & " mixed cells.");
    Cs_Del(cells);
  end Produce_Cells;

  procedure Process_Cells 
              ( idtask,nbequ,nbpts,r : in integer32;
                mtype,perm : in Standard_Integer_Vectors.Link_to_Vector;
                vtx : in Standard_Integer_VecVecs.Link_to_VecVec;
                lft : in Standard_Floating_Vectors.Link_to_Vector;
                mcc : in out Mixed_Subdivision ) is

    labels : Standard_Integer_Vectors.Link_to_Vector;
    mcc_last : Mixed_Subdivision;
    stop : boolean := false;
    cnt : natural32 := 0;

    use Standard_Integer_Vectors;

  begin
    while not stop loop
      labels := Mixed_Labels_Queue.Next;
      if labels = null then -- check if all labels are produced
        if Mixed_Labels_Queue.Stopped     -- production stopped
         then stop := true;   -- all labels have been processed
        end if;
      else
        cnt := cnt + 1;
        put_line("Task " & Multitasking.to_string(idtask)
                         & " processes cell "
                         & Multitasking.to_string(cnt));
        if r < nbequ then
          declare
            mic : constant Mixed_Cell
                := Labels_to_Mixed_Cell(nbequ,r,mtype,labels,Vtx,lft);
          begin
            Append(mcc,mcc_last,mic);
          end;
        else
          declare
            mic : constant Mixed_Cell
                := Labels_to_Mixed_Cell(nbequ,r,mtype,perm,labels,Vtx,lft);
          begin
            Append(mcc,mcc_last,mic);
          end;
        end if;
      end if;
    end loop;
    put_line("Task " & Multitasking.to_string(idtask) & " processed "
                     & Multitasking.to_string(Length_Of(mcc))
                     & " mixed cells.");
  end Process_Cells;

  procedure Pipelined_Mixed_Cells
              ( ntasks,nbequ,nbpts : in integer32;
                ind,cnt : in Standard_Integer_Vectors.Vector;
                support : in Standard_Integer_Vectors.Link_to_Vector;
                r : out integer32;
                mtype,perm : out Standard_Integer_Vectors.Link_to_Vector;
                sub : out Mixed_Subdivision ) is

    stlb : constant double_float := 0.0; -- no stable mv for now...
    mcc : array(2..ntasks) of Mixed_Subdivision;
    idx,sdx,ndx : Standard_Integer_Vectors.Link_to_Vector;
    vtx,spt : Standard_Integer_VecVecs.Link_to_VecVec;
    lft : Standard_Floating_Vectors.Link_to_Vector;
    sub_last : Mixed_Subdivision;

    procedure do_job ( i,n : in integer32 ) is
    begin
      put_line("In do_job with task " & Multitasking.to_string(i));
      if i = 1
       then Produce_Cells(nbequ,nbpts,r,mtype,idx,vtx,lft);
       else Process_Cells(i,nbequ,nbpts,r,mtype,perm,vtx,lft,mcc(i));
      end if;
    end do_job;
    procedure do_jobs is new Multitasking.Reporting_Workers(do_job);

  begin
    Mixed_Labels_Queue.Start;
    mv_upto_pre4mv
      (nbequ,nbpts,ind,cnt,support.all,r,mtype,perm,idx,vtx,sdx,spt,ndx);
    mv_lift(nbequ,nbpts,ind,cnt,support.all,stlb,r,idx,vtx,lft);
    do_jobs(ntasks);
    for i in mcc'range loop
      Concat(sub,sub_last,mcc(i));
    end loop;
    Standard_Integer_Vectors.Clear(idx);
    Standard_Integer_Vectors.Clear(sdx);
    Standard_Integer_Vectors.Clear(ndx);
    Standard_Floating_Vectors.Clear(lft);
    Standard_Integer_VecVecs.Deep_Clear(vtx);
    Standard_Integer_VecVecs.Deep_Clear(spt);
  end Pipelined_Mixed_Cells;

end Pipelined_Labeled_Cells;
