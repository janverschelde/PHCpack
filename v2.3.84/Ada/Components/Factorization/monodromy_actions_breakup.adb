with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Witness_Sets;                       use Witness_Sets;
with Sampling_Machine;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists_io;              use Sample_Point_Lists_io;
with Monodromy_Group_Actions_io;         use Monodromy_Group_Actions_io;

package body Monodromy_Actions_Breakup is

-- AUXILIARIES :

  function Compare_Labels ( file : in file_type; tol : in double_float;
                            sps1,sps2 : in Standard_Sample_List )
                          return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   For two lists of the same samples, the labels of the solutions
  --   are compared with each other.  The i-th component of the vector
  --   on return indicates the position of the i-th sample of sps1 in
  --   the list sps2.  When the i-th sample of sps1 does not occur,
  --   then an error message is printed, and the bug is "fixed" by
  --   assigning i in the i-th entry on return.

    res : Standard_Natural_Vectors.Vector(1..integer32(Length_Of(sps1)));
    tmp : Standard_Sample_List := sps1;

  begin
    for i in res'range loop
      res(i) := Is_In(sps2,tol,Head_Of(tmp));
      put(file,"Sample point "); put(file,i,1);
      put(file," is mapped to "); put(file,res(i),1);
      if res(i) = 0 then
        put(file,"  Oops, bug fixed");
        res(i) := natural32(i);
      end if;
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Compare_Labels;

  procedure Generate_Monodromy_Action
               ( file : in file_type;
                 sps,newsps : in Standard_Sample_List; 
                 newhyp : in Standard_Complex_VecVecs.VecVec;
                 tol : in double_float;
                 rel : out Standard_Natural_Vectors.Vector ) is 
 
  -- DESCRIPTION :
  --   Generates one new element of the monodromy group and applies
  --   this to the decomposition.

  -- ON ENTRY :
  --   file      for intermediate results and diagnostics;
  --   sps       list of generic points, Length_Of(sps) = rep'last;
  --   newsps    list of target generic points, on slices newhyp;
  --   newhyp    slices for target generic points;
  --   tol       tolerance to decide inclusion.

  -- ON RETURN :
  --   rel       relations between the points.

    mapsps,mapsps_last : Standard_Sample_List;

  begin
    Sample(sps,newhyp,mapsps,mapsps_last);
    rel := Compare_Labels(file,tol,mapsps,newsps);
    Deep_Clear(mapsps);
  end Generate_Monodromy_Action;

  procedure Generate_Monodromy_Action
               ( file : in file_type;
                 ind : in Standard_Natural_Vectors.Vector;
                 sps,newsps : in Standard_Sample_List; 
                 newhyp : in Standard_Complex_VecVecs.VecVec;
                 tol : in double_float;
                 rel : out Standard_Natural_Vectors.Vector ) is 
 
  -- DESCRIPTION :
  --   Generates one new element of the monodromy group and applies
  --   this to the decomposition.

  -- ON ENTRY :
  --   file      for intermediate results and diagnostics;
  --   ind       representatives for active generic points in sps;
  --   sps       list of generic points, Length_Of(sps) = rep'last;
  --   newsps    list of target generic points, on slices newhyp;
  --   newhyp    slices for target generic points;
  --   tol       tolerance to decide inclusion.

  -- ON RETURN :
  --   rel       relations between the points.

    mapsps,mapsps_last,tmp : Standard_Sample_List;

  begin
    Sample(sps,newhyp,mapsps,mapsps_last);
    tmp := mapsps;
    for i in 1..integer32(Length_Of(sps)) loop
      rel(integer32(ind(i))) := Is_In(newsps,tol,Head_Of(tmp));
      put(file,"Sample point "); put(file,ind(i),1);
      put(file," is mapped to "); put(file,rel(integer32(ind(i))),1);
      if rel(integer32(ind(i))) = 0 then
        put(file,"  Oops, bug fixed");
        rel(integer32(ind(i))) := ind(i);
      end if; 
      new_line(file);
      tmp := Tail_Of(tmp);
    end loop; 
    Deep_Clear(mapsps);
  end Generate_Monodromy_Action;

  procedure Data_Management
                ( file : in file_type;
                  rel : in Standard_Natural_Vectors.Vector;
                  ic : in out Irreducible_Components;
                  n1,n2,cnt : in out natural32 ) is

  -- DESCRIPTION :
  --   Manages the result of one monodromy iteration.

  -- ON ENTRY :
  --   file       for intermediate output and diagnostics;
  --   rel        rel(i) is related to point i;
  --   ic         current grouping of generic points;
  --   n1         previous number of representative sets in ic;
  --   n2         current number of sets in ic;
  --   cnt        number of consecutive stabilizations.

  -- ON RETURN :
  --   ic         new grouping of the generic points;
  --   n1         equals the value of the incoming n2;
  --   n2         cardinality of the newly computed ic;
  --   cnt        updated counter.

  begin
    Act(ic,rel);
    n2 := Cardinality(ic);
    put(file,"The decomposition : "); put(file,ic);
    put(file,"Degrees : "); put(file,Degrees(ic)); new_line(file);
    put(file,"Cardinality ");
    put(file,n1,1); put(file," -> "); put(file,n2,1);
    if n1 = n2 then
      put(file," is stable with count = ");
      cnt := cnt+1;
      put(file,cnt,1); 
      put(file,".");
    else
      put(file," drops.");
      cnt := 0;
    end if;
    new_line(file);
    n1 := n2;
  end Data_Management;

-- TARGET ROUTINES :

  procedure Generic_Points
              ( sps : in Standard_Sample_List;
                labels : in Standard_Natural_Vectors.Vector;
                gp,gp_last : in out Standard_Sample_List ) is

    tmp : Standard_Sample_List := sps;
    ind : integer32 := 1;

  begin
    for i in 1..Length_Of(sps) loop
      if labels(ind) = i
       then Append(gp,gp_last,Head_Of(tmp)); ind := ind + 1;
      end if;
      exit when (ind > labels'last);
      tmp := Tail_Of(tmp);
    end loop;
  end Generic_Points;

  procedure Breakup 
              ( file : in file_type; sps : in Standard_Sample_List;
                threshold : in natural32; tol : in double_float;
                ic : out Irreducible_Components; nbit : out natural32 ) is

    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : Standard_Complex_VecVecs.VecVec(1..dim);
    mapsps,mapsps_last,prevsps : Standard_Sample_List;
    n1 : natural32 := Length_Of(sps);
    rel : Standard_Natural_Vectors.Vector(1..integer32(n1));
    cnt,n2,nit : natural32 := 0;
    grid,grid_last,tmp : Standard_Sample_Grid;

  begin
    ic := Create(rel'last);
    loop
      nit := nit+1;
      newhyp := Random_Hyperplanes(natural32(dim),natural32(n));  
      declare                                      -- create new sample list
        newsps,newsps_last : Standard_Sample_List;
      begin
        Sample(sps,newhyp,newsps,newsps_last);
        Sampling_Machine.Change_Slices(newhyp);
        Sample(newsps,hyp,mapsps,mapsps_last);
        rel := Compare_Labels(file,tol,sps,mapsps);
        Deep_Clear(mapsps);
        Data_Management(file,rel,ic,n1,n2,cnt);
        exit when (cnt = threshold);
        tmp := grid;
        Sampling_Machine.Change_Slices(newhyp);
        while not Is_Null(tmp) loop          -- from new to all previous ones
          nit := nit + 1;
          prevsps := Head_Of(tmp);
          prevhyp := Hyperplane_Sections(Head_Of(prevsps));
          Sample(newsps,prevhyp,mapsps,mapsps_last);
          rel := Compare_Labels(file,tol,prevsps,mapsps);
          Deep_Clear(mapsps);
          Data_Management(file,rel,ic,n1,n2,cnt);
          exit when ((cnt = threshold) or (n2 = 1));
          tmp := Tail_Of(tmp);
        end loop;
        Append(grid,grid_last,newsps);
      end;
      Sampling_Machine.Change_Slices(hyp);
      exit when ((cnt = threshold) or (n2 = 1));
    end loop;
    nbit := nit;
  end Breakup;

  procedure Breakup
              ( file : in file_type; sps : in Standard_Sample_List;
                threshold : in natural32; tol : in double_float;
                ic : out Irreducible_Components; nbit : out natural32;
                grid,grid_last : in out Standard_Sample_Grid ) is

    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Hyperplane_Sections(Head_Of(sps));
    n1 : natural32 := Length_Of(sps);
    rel : Standard_Natural_Vectors.Vector(1..integer32(n1));
    cnt,n2,nit : natural32 := 0;

  begin
    ic := Create(rel'last);
    loop
      nit := nit+1;
      declare
        parhyp,newhyp : Standard_Complex_VecVecs.VecVec(1..dim);
        parsps,parsps_last,newsps,newsps_last : Standard_Sample_List;
      begin
        for i in parhyp'range loop
          parhyp(i) := new Standard_Complex_Vectors.Vector'(hyp(i).all);
          parhyp(i)(0) := Random1;
        end loop;
        Sample(sps,parhyp,parsps,parsps_last);
        put_line(file,"Summaries of the samples on parallel slices : ");
        Write_Summaries(file,parsps);
        Append(grid,grid_last,parsps);
        newhyp := Random_Hyperplanes(natural32(dim),natural32(n));
        Sampling_Machine.Change_Slices(parhyp);
        Sample(parsps,newhyp,newsps,newsps_last); 
        put_line(file,"Summaries of the samples on the new slices : ");
        Write_Summaries(file,newsps);
        Sampling_Machine.Change_Slices(hyp);
        Generate_Monodromy_Action(file,sps,newsps,newhyp,tol,rel);
        Deep_Clear(newsps);
      end;
      Data_Management(file,rel,ic,n1,n2,cnt);
      exit when ((cnt = threshold) or (n2 = 1));
    end loop;
    nbit := nit;
  end Breakup;

end Monodromy_Actions_Breakup;
