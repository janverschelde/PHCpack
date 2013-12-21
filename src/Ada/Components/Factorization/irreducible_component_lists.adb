with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Standard_Natural_Vectors;
with Standard_Complex_VecVecs;
with Multprec_Complex_VecVecs;
with Witness_Sets;                       use Witness_Sets;
with Sampling_Machine;
with Sample_Points;                      use Sample_Points;
with Projection_Operators;               use Projection_Operators;
with Interpolation_Filters;              use Interpolation_Filters;
with Span_of_Component;                  use Span_of_Component;
with Irreducible_Component_Creators;     use Irreducible_Component_Creators;
with Monodromy_Group_Actions;            use Monodromy_Group_Actions;
with Monodromy_Actions_Breakup;          use Monodromy_Actions_Breakup;
with Standard_Divided_Differences;
with Multprec_Divided_Differences;
with Standard_Trace_Interpolators;
with Multprec_Trace_Interpolators;

package body Irreducible_Component_Lists is

-- AUXILIARIES FOR LABEL MANAGEMENT :

  procedure Add_Label_and_Point
              ( deco : in out Standard_Irreducible_Component_List;
                i : in integer32; k : in natural32;
                spt : in Standard_Sample ) is

  -- DESCRIPTION :
  --   Adds the label k and the sample point spt to the i-th component
  --   in the list deco.

    tmp : Standard_Irreducible_Component_List := deco;
    c : Standard_Irreducible_Component;
    ind : integer32;

  begin
    for j in 1..(i-1) loop
      exit when Is_Null(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    if not Is_Null(tmp) then
      c := Head_Of(tmp);
      Add_Label(c,k,ind);
      Add_Point(c,spt);
      Set_Head(tmp,c);
    end if;
  end Add_Label_and_Point;

  procedure Add_Label_and_Point
              ( deco : in out Multprec_Irreducible_Component_List;
                i : in integer32; k : in natural32;
                spt : in Standard_Sample ) is

  -- DESCRIPTION :
  --   Adds the label k and the sample point spt to the i-th component
  --   in the list deco.

    tmp : Multprec_Irreducible_Component_List := deco;
    c : Multprec_Irreducible_Component;
    ind : integer32;

  begin
    for j in 1..(i-1) loop
      exit when Is_Null(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    if not Is_Null(tmp) then
      c := Head_Of(tmp);
      Add_Label(c,k,ind);
      Add_Point(c,spt);
      Set_Head(tmp,c);
    end if;
  end Add_Label_and_Point;

-- AUXILIARIES FOR INTERMEDIATE OUTPUT :

  procedure Write_Conclusion
                ( file : in file_type;
                  res : in boolean; tol : in double_float ) is

  -- DESCRIPTION :
  --   This procedure ensures a uniform output format in writing the
  --   conclusion of a membership test.

  begin
    if res
     then put(file,"  <"); put(file,tol,3); put_line(file,"  success.");
     else put(file,"  >"); put(file,tol,3); put_line(file,"  failure.");
    end if;
  end Write_Conclusion;

  procedure Write_Banner ( file : in file_type; i : in natural32 ) is
  begin
    put(file,"Processing Generic Point ");
    put(file,i,1);
    put_line(file," :");
  end Write_Banner;

-- CREATORS WITH INCREMENTAL INTERPOLATION :

  procedure Standard_Massive_Interpolate
              ( file : in file_type; full_output : in boolean;
                sps : in Standard_Sample_List;
                stoptol,membtol : in double_float;
                deco,deco_last : out Standard_Irreducible_Component_List ) is

    p : Standard_Projector;
    maxdeg : natural32;
    tmp : Standard_Sample_List;
    spt : Standard_Sample;

  begin
    if not Is_Null(sps) then
      spt := Head_Of(sps);
      p := Create(natural32(Number_of_Slices(spt))+1,
                  natural32(Number_of_Variables(spt)));
      tmp := sps;
      maxdeg := Length_Of(sps);
      for i in 1..maxdeg loop
        Write_Banner(file,i);
        spt := Head_Of(tmp);
        declare
          f : Standard_Filter;
          c : Standard_Irreducible_Component;
        begin
          Standard_Interpolate(file,full_output,spt,maxdeg,stoptol,p,f);
          c := Create(f);
          Append(deco,deco_last,c);
        end;
        tmp := Tail_Of(tmp);
      end loop;   
    end if;
  end Standard_Massive_Interpolate;

  procedure Standard_Incremental_Interpolate
              ( file : in file_type; full_output : in boolean;
                sps : in Standard_Sample_List;
                stoptol,membtol : in double_float;
                deco,deco_last : out Standard_Irreducible_Component_List ) is

    maxdeg : constant natural32 := Length_Of(sps);
    tmp : Standard_Sample_List := sps;
    spt : Standard_Sample := Head_Of(sps);
    ind : integer32;
    p : constant Standard_Projector
      := Create(natural32(Number_of_Slices(Head_Of(sps)))+1,
                natural32(Number_of_Variables(Head_Of(sps))));

  begin
    deco := Standard_Null_List;
    for i in 1..maxdeg loop
      Write_Banner(file,i);
      spt := Head_Of(tmp);
      ind := integer32(On_Component(file,deco,Sample_Point(spt).v,membtol));
      if ind = 0 then
        declare
          f : Standard_Filter;
          c : Standard_Irreducible_Component;
        begin
          Standard_Interpolate(file,full_output,spt,maxdeg,stoptol,p,f);
          c := Create(f);
          Initialize_Labels(c,integer32(Degree(c)));
          Add_Label(c,i,ind);
          Append(deco,deco_last,c);
        end;
      else
        Add_Label_and_Point(deco,ind,i,spt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;   
  end Standard_Incremental_Interpolate;

  procedure Standard_Incremental_Interpolate_with_Span
              ( file : in file_type; full_output : in boolean;
                sps : in Standard_Sample_List;
                stoptol,membtol : in double_float;
                deco,deco_last : out Standard_Irreducible_Component_List ) is

    maxdeg : constant natural32 := Length_Of(sps);
    tmp : Standard_Sample_List := sps;
    spt : Standard_Sample;
    ind : integer32;
    p : constant Standard_Projector
      := Create(natural32(Number_of_Slices(Head_Of(sps)))+1,
                natural32(Number_of_Variables(Head_Of(sps))));

  begin
    deco := Standard_Null_List;
    for i in 1..maxdeg loop
      Write_Banner(file,i);
      spt := Head_Of(tmp);
      ind := integer32(On_Component(file,deco,Sample_Point(spt).v,membtol));
      if ind = 0 then
        declare
          f : Standard_Filter;
          s : Standard_Span;
          c : Standard_Irreducible_Component;
        begin
          Standard_Interpolate(file,full_output,spt,maxdeg,stoptol,p,f,s);
          c := Create(f,s);
          Initialize_Labels(c,integer32(Degree(c)));
          Add_Label(c,i,ind);
          Add_Point(c,spt);
          Append(deco,deco_last,c);
        end;
      else
        Add_Label_and_Point(deco,ind,i,spt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;   
  end Standard_Incremental_Interpolate_with_Span;

  procedure Standard_Incremental_Central_Interpolate
              ( file : in file_type; full_output : in boolean;
                sps : in Standard_Sample_List;
                stoptol,membtol : in double_float;
                deco,deco_last : out Standard_Irreducible_Component_List ) is

    maxdeg : constant natural32 := Length_Of(sps);
    tmp : Standard_Sample_List := sps;
    spt : Standard_Sample;
    ind : integer32;
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    hyps : Standard_Complex_VecVecs.VecVec(1..dim+1);

  begin
    deco := Standard_Null_List;
    for i in hyps'range loop
      hyps(i) := new Standard_Complex_Vectors.Vector'(Random_Vector(0,n));
    end loop;
    for i in 1..maxdeg loop
      Write_Banner(file,i);
      spt := Head_Of(tmp);
      ind := integer32(On_Component(file,deco,Sample_Point(spt).v,membtol));
      if ind = 0 then
        declare
          f : Standard_Filter;
          p : constant Standard_Projector := Create(hyps);
          s : Standard_Span;
          c : Standard_Irreducible_Component;
        begin
          Standard_Central_Interpolate
            (file,full_output,spt,maxdeg,stoptol,p,f,s);
          c := Create(f,s);
          Initialize_Labels(c,integer32(Degree(c)));
          Add_Label(c,i,ind);
          Add_Point(c,spt);
          Append(deco,deco_last,c);
         end;
      else
        Add_Label_and_Point(deco,ind,i,spt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;   
  end Standard_Incremental_Central_Interpolate;

  procedure Multprec_Massive_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 size : in natural32; stoptol,membtol : in double_float;
                 deco,deco_last : out Multprec_Irreducible_Component_List ) is

    p : Multprec_Projector;
    maxdeg : natural32;
    tmp : Standard_Sample_List;
    spt : Standard_Sample;

  begin
    if not Is_Null(sps) then
      spt := Head_Of(sps);
      p := Create(natural32(Number_of_Slices(spt))+1,
                  natural32(Number_of_Variables(spt)),size);
      maxdeg := Length_Of(sps);
      tmp := sps;
      for i in 1..maxdeg loop
        Write_Banner(file,i);
        spt := Head_Of(tmp);
        declare
          f : Multprec_Filter;
          c : Multprec_Irreducible_Component;
        begin
          Multprec_Interpolate
            (file,full_output,spt,maxdeg,size,stoptol,p,f);
          c := Create(f);
          Append(deco,deco_last,c);
        end;
        tmp := Tail_Of(tmp);
      end loop;   
    end if;
  end Multprec_Massive_Interpolate;

  procedure Multprec_Incremental_Interpolate
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 size : in natural32; stoptol,membtol : in double_float;
                 deco,deco_last : out Multprec_Irreducible_Component_List ) is

    maxdeg : constant natural32 := Length_Of(sps);
    tmp : Standard_Sample_List := sps;
    stspt : Standard_Sample;
    mpspt : Multprec_Sample; 
    ind : integer32;
    p : constant Multprec_Projector
      := Create(natural32(Number_of_Slices(Head_Of(sps)))+1,
                natural32(Number_of_Variables(Head_Of(sps))),size); 

  begin
    deco := Multprec_Null_List;
    for i in 1..maxdeg loop
      Write_Banner(file,i);
      stspt := Head_Of(tmp);
      Refine(file,full_output,stspt,mpspt);
      ind := integer32(On_Component(file,deco,Sample_Point(mpspt).v,membtol));
      if ind = 0 then
        declare
          f : Multprec_Filter;
          c : Multprec_Irreducible_Component;
        begin
          Multprec_Interpolate
            (file,full_output,mpspt,maxdeg,size,stoptol,p,f);
          c := Create(f);
          Initialize_Labels(c,integer32(Degree(c)));
          Add_Label(c,i,ind);
          Add_Point(c,stspt);
          Append(deco,deco_last,c);
        end;
      else
        Add_Label_and_Point(deco,ind,i,stspt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;   
  end Multprec_Incremental_Interpolate;

  procedure Multprec_Incremental_Interpolate_with_Span
               ( file : in file_type; full_output : in boolean;
                 sps : in Standard_Sample_List;
                 size : in natural32; stoptol,membtol : in double_float;
                 deco,deco_last : out Multprec_Irreducible_Component_List ) is

    maxdeg : constant natural32 := Length_Of(sps);
    tmp : Standard_Sample_List := sps;
    stspt : Standard_Sample;
    mpspt : Multprec_Sample;
    ind : integer32;
    p : constant Multprec_Projector
      := Create(natural32(Number_of_Slices(Head_Of(sps)))+1,
                natural32(Number_of_Variables(Head_Of(sps))),size);

  begin
    deco := Multprec_Null_List;
    for i in 1..maxdeg loop
      Write_Banner(file,i);
      stspt := Head_Of(tmp);
      Refine(file,full_output,stspt,mpspt);
      ind := integer32(On_Component(file,deco,Sample_Point(mpspt).v,membtol));
      if ind = 0 then
        declare
          f : Multprec_Filter;
          s : Multprec_Span;
          c : Multprec_Irreducible_Component;
        begin
          Multprec_Interpolate
            (file,full_output,mpspt,maxdeg,size,stoptol,p,f,s);
          c := Create(f,s);
          Initialize_Labels(c,integer32(Degree(c)));
          Add_Label(c,i,ind);
          Add_Point(c,stspt);
          Append(deco,deco_last,c);
        end;
      else
        Add_Label_and_Point(deco,ind,i,stspt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;   
  end Multprec_Incremental_Interpolate_with_Span;

  procedure Multprec_Incremental_Central_Interpolate
              ( file : in file_type; full_output : in boolean;
                sps : in Standard_Sample_List;
                size : in natural32; stoptol,membtol : in double_float;
                deco,deco_last : out Multprec_Irreducible_Component_List ) is

    maxdeg : constant natural32 := Length_Of(sps);
    tmp : Standard_Sample_List := sps;
    stspt : Standard_Sample;
    mpspt : Multprec_Sample;
    ind : integer32;
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    hyps : Multprec_Complex_VecVecs.VecVec(1..dim+1);
   
  begin
    deco := Multprec_Null_List;
    for i in hyps'range loop
      hyps(i) := new Multprec_Complex_Vectors.Vector'
                        (Create(Random_Vector(0,n)));
    end loop;
    for i in 1..maxdeg loop
      Write_Banner(file,i);
      stspt := Head_Of(tmp);
      Refine(file,full_output,stspt,mpspt);
      ind := integer32(On_Component(file,deco,Sample_Point(mpspt).v,membtol));
      if ind = 0 then
        declare
          p : constant Multprec_Projector := Create(hyps);
          f : Multprec_Filter;
          s : Multprec_Span;
          c : Multprec_Irreducible_Component;
        begin
          Multprec_Central_Interpolate
            (file,full_output,mpspt,maxdeg,size,stoptol,p,f,s);
          c := Create(f,s);
          Initialize_Labels(c,integer32(Degree(c)));
          Add_Label(c,i,ind);
          Add_Point(c,stspt);
          Append(deco,deco_last,c);
        end;
      else
        Add_Label_and_Point(deco,ind,i,stspt);
      end if;
      tmp := Tail_Of(tmp);
    end loop;   
  end Multprec_Incremental_Central_Interpolate;

-- AUXILIARY CREATOR TO MONODROMY BREAKUP :

  procedure Create_Component_Lists
               ( ic : in Monodromy_Group_Actions.Irreducible_Components;
                 sps : in Standard_Sample_List;
                 deco,deco_last : out Standard_Irreducible_Component_List ) is

  -- DESCRIPTION :
  --   Converts the representation used in Monodromy_Group_Actions into
  --   a list of irreducible components.

  begin
    for i in 1..Sum_of_Degrees(ic) loop
      if not Empty(ic,i) then
        declare
          lab : constant Standard_Natural_Vectors.Vector := Component(ic,i);
          gp,gp_last : Standard_Sample_List;
          c : Standard_Irreducible_Component;
        begin
          Initialize_Labels(c,lab);
          Generic_Points(sps,lab,gp,gp_last);
          Add_Points(c,gp,gp_last);
          Append(deco,deco_last,c);
        end;
      end if;
    end loop;
  end Create_Component_Lists;

-- CREATORS TO PREDICT BREAKUP WITH MONODROMY LOOPS :

  procedure Monodromy_Breakup
               ( file : in file_type; sps : in Standard_Sample_List;
                 threshold : in natural32; tol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List;
                 nit : out natural32 ) is

    ic : Monodromy_Group_Actions.Irreducible_Components;

  begin
    Breakup(file,sps,threshold,tol,ic,nit);
    Create_Component_Lists(ic,sps,deco,deco_last);
  end Monodromy_Breakup;

  procedure New_Monodromy_Breakup
               ( file : in file_type; sps : in Standard_Sample_List;
                 threshold : in natural32; tol : in double_float;
                 deco,deco_last : out Standard_Irreducible_Component_List;
                 nit : out natural32 ) is

    ic : Monodromy_Group_Actions.Irreducible_Components;
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    dim : constant integer32 := Number_of_Slices(Head_Of(sps));
    hyp : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : Standard_Complex_VecVecs.VecVec(1..dim);
    mapsps,mapsps_last,prevsps : Standard_Sample_List;
    n1 : natural32 := Length_Of(sps);
    rel : Standard_Natural_Vectors.Vector(1..integer32(n1));
    cnt,n2,nbit : natural32 := 0;
    grid,grid_last,tmp : Standard_Sample_Grid;

  begin
   -- Breakup(file,sps,threshold,tol,ic,nit);
    ic := Create(rel'last);
    nit := 0;
    loop
      nbit := nbit+1;
      newhyp := Random_Hyperplanes(natural32(dim),natural32(n));
      declare                          --  create new sample list
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
    nit := nbit;
    Create_Component_Lists(ic,sps,deco,deco_last);
  end New_Monodromy_Breakup;

  procedure Monodromy_Breakup
              ( file : in file_type; sps : in Standard_Sample_List;
                threshold : in natural32; tol : in double_float;
                deco,deco_last : out Standard_Irreducible_Component_List;
                nit : out natural32;
                grid,grid_last : in out Standard_Sample_Grid ) is

    ic : Monodromy_Group_Actions.Irreducible_Components;

  begin
    Breakup(file,sps,threshold,tol,ic,nit,grid,grid_last);
    Create_Component_Lists(ic,sps,deco,deco_last);
  end Monodromy_Breakup;

  procedure Distribute_Points 
               ( deco : in out Standard_Irreducible_Component_List;
                 sps : in Standard_Sample_List ) is

    tmp : Standard_Irreducible_Component_List := deco;

  begin
    while not Is_Null(tmp) loop
      declare
        c : Standard_Irreducible_Component := Head_Of(tmp);
      begin
        Select_Labeled_Points(c,sps);
        Set_Head(tmp,c);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Distribute_Points;

  procedure Distribute_Points 
               ( deco : in out Multprec_Irreducible_Component_List;
                 sps : in Standard_Sample_List ) is

    tmp : Multprec_Irreducible_Component_List := deco;

  begin
    while not Is_Null(tmp) loop
      declare
        c : Multprec_Irreducible_Component := Head_Of(tmp);
      begin
        Select_Labeled_Points(c,sps);
        Set_Head(tmp,c);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Distribute_Points;

-- CREATORS WITH BOOTSTRAPPING NEWTON :

  procedure Standard_Newton_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix ) is

    tmp : Standard_Irreducible_Component_List := deco;

  begin
    for i in 1..integer32(Length_Of(tmp)) loop
      declare
        c : constant Standard_Irreducible_Component := Head_Of(tmp);
        d : constant natural32 := Degree(c);
        sps : constant Standard_Sample_List := Points(c);
        dim : constant integer32 := Number_of_Slices(Head_Of(sps));
        q1 : Standard_Divided_Differences.Newton_Interpolator1;
        qk : Standard_Divided_Differences.Newton_Taylor_Form;
      begin
        put(file,"Interpolating at component ");
        put(file,i,1); put(file," of degree ");
        put(file,d,1); put(file," and dimension ");
        put(file,dim,1); put_line(file," :");
        numres(i,1) := double_float(d);
        if Length_Of(sps) > 1 then
          if dim = 1 then
            Standard_Newton_Interpolate1
              (file,sps,q1,numres(i,2),numres(i,3),numres(i,4),numres(i,5));
          else
            Standard_Newton_Interpolate
              (file,sps,qk,numres(i,2),numres(i,3),numres(i,4),numres(i,5));
          end if;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Newton_Interpolate;

  procedure Multprec_Newton_Interpolate
               ( file : in file_type; size : in natural32;
                 deco : in out Multprec_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix ) is

    tmp : Multprec_Irreducible_Component_List := deco;

  begin
    for i in 1..integer32(Length_Of(tmp)) loop
      declare
        c : constant Multprec_Irreducible_Component := Head_Of(tmp);
        d : constant natural32 := Degree(c);
        sps : constant Standard_Sample_List := Points(c);
        dim : constant integer32 := Number_of_Slices(Head_Of(sps));
        q1 : Multprec_Divided_Differences.Newton_Interpolator1;
        qk : Multprec_Divided_Differences.Newton_Taylor_Form;
      begin
        put(file,"Interpolating at component ");
        put(file,i,1); put(file," of degree ");
        put(file,d,1); put_line(file," :");
        numres(i,1) := double_float(d);
        if Length_Of(sps) > 1 then
          if dim = 1 then
            Multprec_Newton_Interpolate1
              (file,sps,size,q1,numres(i,2),numres(i,3),
               numres(i,4),numres(i,5));
          else
            Multprec_Newton_Interpolate
              (file,sps,size,qk,numres(i,2),numres(i,3),
               numres(i,4),numres(i,5));
          end if;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Multprec_Newton_Interpolate;

-- CREATORS WITH TRACE FORMS :

  procedure Standard_Trace_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix ) is

    tmp : Standard_Irreducible_Component_List := deco;

  begin
    for i in 1..integer32(Length_Of(tmp)) loop
      declare
        c : constant Standard_Irreducible_Component := Head_Of(tmp);
        d : constant natural32 := Degree(c);
        sps : constant Standard_Sample_List := Points(c);
        dim : constant integer32 := Number_of_Slices(Head_Of(sps));
        q1 : Standard_Trace_Interpolators.Trace_Interpolator1;
        qk : Standard_Trace_Interpolators.Trace_Interpolator;
      begin
        put(file,"Interpolating at component ");
        put(file,i,1); put(file," of degree ");
        put(file,d,1); put(file," and dimension ");
        put(file,dim,1); put_line(file," :");
        numres(i,1) := double_float(d);
        if Length_Of(sps) > 1 then
          if dim = 1 then
            Standard_Trace_Interpolate1
              (file,sps,q1,numres(i,2),numres(i,3),numres(i,4),numres(i,5));
          else
            Standard_Trace_Interpolate
              (file,sps,qk,numres(i,2),numres(i,3),numres(i,4),numres(i,5));
          end if;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Trace_Interpolate;

  procedure Multprec_Trace_Interpolate
               ( file : in file_type; size : in natural32;
                 deco : in out Multprec_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix ) is

    tmp : Multprec_Irreducible_Component_List := deco;

  begin
    for i in 1..integer32(Length_Of(tmp)) loop
      declare
        c : constant Multprec_Irreducible_Component := Head_Of(tmp);
        d : constant natural32 := Degree(c);
        sps : constant Standard_Sample_List := Points(c);
        dim : constant integer32 := Number_of_Slices(Head_Of(sps));
        q1 : Multprec_Trace_Interpolators.Trace_Interpolator1;
        qk : Multprec_Trace_Interpolators.Trace_Interpolator;
      begin
        put(file,"Interpolating at component ");
        put(file,i,1); put(file," of degree ");
        put(file,d,1); put_line(file," :");
        numres(i,1) := double_float(d);
        if Length_Of(sps) > 1 then 
          if dim = 1 then
            Multprec_Trace_Interpolate1
              (file,sps,size,q1,numres(i,2),numres(i,3),
               numres(i,4),numres(i,5));
          else
            Multprec_Trace_Interpolate
              (file,sps,size,qk,numres(i,2),numres(i,3),
               numres(i,4),numres(i,5));
          end if;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Multprec_Trace_Interpolate;

  procedure Standard_Power_Trace_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix ) is

    tmp : Standard_Irreducible_Component_List := deco;

  begin
    for i in 1..integer32(Length_Of(tmp)) loop
      declare
        c : constant Standard_Irreducible_Component := Head_Of(tmp);
        d : constant natural32 := Degree(c);
        sps : constant Standard_Sample_List := Points(c);
        dim : constant integer32 := Number_of_Slices(Head_Of(sps));
        q1 : Standard_Trace_Interpolators.Trace_Interpolator1;
        qk : Standard_Trace_Interpolators.Trace_Interpolator;
      begin
        put(file,"Interpolating at component ");
        put(file,i,1); put(file," of degree ");
        put(file,d,1); put(file," and dimension ");
        put(file,dim,1); put_line(file," :");
        numres(i,1) := double_float(d);
        if Length_Of(sps) > 1 then
          if dim = 1 then
            Standard_Power_Trace_Interpolate1
              (file,sps,q1,numres(i,2),numres(i,3),numres(i,4),numres(i,5));
          else
            put_line(file,"Higher dimensional power traces not done yet.");
            put_line(file,"Invoking other full trace interpolation.");
            Standard_Trace_Interpolate
              (file,sps,qk,numres(i,2),numres(i,3),
               numres(i,4),numres(i,5));
          end if;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Power_Trace_Interpolate;

  procedure Multprec_Power_Trace_Interpolate
              ( file : in file_type; size : in natural32;
                deco : in out Multprec_Irreducible_Component_List;
                numres : out Standard_Floating_Matrices.Matrix ) is

    tmp : Multprec_Irreducible_Component_List := deco;

  begin
    for i in 1..integer32(Length_Of(tmp)) loop
      declare
        c : constant Multprec_Irreducible_Component := Head_Of(tmp);
        d : constant natural32 := Degree(c);
        sps : constant Standard_Sample_List := Points(c);
        dim : constant integer32 := Number_of_Slices(Head_Of(sps));
        q1 : Multprec_Trace_Interpolators.Trace_Interpolator1;
        qk : Multprec_Trace_Interpolators.Trace_Interpolator;
      begin
        put(file,"Interpolating at component ");
        put(file,i,1); put(file," of degree ");
        put(file,d,1); put_line(file," :");
        numres(i,1) := double_float(d);
        if Length_Of(sps) > 1 then
          if dim = 1 then
            Multprec_Power_Trace_Interpolate1
             (file,sps,size,q1,numres(i,2),numres(i,3),
              numres(i,4),numres(i,5));
          else
            put_line(file,"Higher dimensional power traces not done yet.");
            put_line(file,"Invoking other full trace interpolation.");
            Multprec_Trace_Interpolate
              (file,sps,size,qk,numres(i,2),numres(i,3),
               numres(i,4),numres(i,5));
          end if;
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Multprec_Power_Trace_Interpolate;

  procedure Standard_Linear_Trace_Interpolate
               ( file : in file_type;
                 deco : in out Standard_Irreducible_Component_List;
                 numres : out Standard_Floating_Matrices.Matrix ) is

    tmp : Standard_Irreducible_Component_List := deco;

  begin
    for i in 1..integer32(Length_Of(tmp)) loop
      declare
        c : constant Standard_Irreducible_Component := Head_Of(tmp);
        d : constant natural32 := Degree(c);
        sps : constant Standard_Sample_List := Points(c);
       -- dim : constant integer32 := Number_of_Slices(Head_Of(sps));
        q1 : Standard_Complex_Vectors.Vector(0..1);
      begin
        put(file,"Interpolating at component ");
        put(file,i,1); put(file," of degree ");
        put(file,d,1); put_line(file," :");
        numres(i,1) := double_float(d);
        if Length_Of(sps) > 1 then
          Standard_Linear_Trace_Interpolate
            (file,sps,q1,numres(i,2),numres(i,3),numres(i,4),numres(i,5));
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Linear_Trace_Interpolate;

-- SELECTORS :

  function On_Component ( L : Standard_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32 is

    tmp : Standard_Irreducible_Component_List := L;

  begin
    for i in 1..Length_Of(L) loop
      if On_Component(Head_Of(tmp),x,tol)
       then return i;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end On_Component;

  function On_Component ( file : in file_type;
                          L : Standard_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32 is

    tmp : Standard_Irreducible_Component_List := L;

  begin
    for i in 1..Length_Of(L) loop
      put(file,"Membership test at component "); put(file,i,1);
      put_line(file," : "); 
      if On_Component(file,Head_Of(tmp),x,tol)
       then return i;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end On_Component;

  function On_Component ( L : Multprec_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32 is

    tmp : Multprec_Irreducible_Component_List := L;

  begin
    for i in 1..Length_Of(L) loop
      if On_Component(Head_Of(tmp),x,tol)
       then return i;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end On_Component;

  function On_Component ( file : in file_type;
                          L : Multprec_Irreducible_Component_List;
                          x : Standard_Complex_Vectors.Vector;
                          tol : double_float ) return natural32 is

    res : natural32;
    mpx : Multprec_Complex_Vectors.Vector(x'range) := Create(x);

  begin
    res := On_Component(file,L,mpx,tol);
    Multprec_Complex_Vectors.Clear(mpx);
    return res;
  end On_Component;

  function On_Component ( L : Multprec_Irreducible_Component_List;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return natural32 is

    tmp : Multprec_Irreducible_Component_List := L;

  begin
    for i in 1..Length_Of(L) loop
      if On_Component(Head_Of(tmp),x,tol)
       then return i;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end On_Component;

  function On_Component ( file : in file_type;
                          L : Multprec_Irreducible_Component_List;
                          x : Multprec_Complex_Vectors.Vector;
                          tol : double_float ) return natural32 is

    tmp : Multprec_Irreducible_Component_List := L;

  begin
    for i in 1..Length_Of(L) loop
      put(file,"Membership test at component "); put(file,i,1);
      put_line(file," : "); 
      if On_Component(file,Head_Of(tmp),x,tol)
       then return i;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end On_Component;

  function Homotopy_Filter ( L : Standard_Irreducible_Component_List;
                             x : Standard_Complex_Vectors.Vector;
                             tol : double_float ) return natural32 is

    tmp : Standard_Irreducible_Component_List := L;

  begin
    for i in 1..Length_Of(L) loop
      if Homotopy_Membership_Test(Head_Of(tmp),x,tol)
       then return i;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end Homotopy_Filter;

  function Homotopy_Filter ( file : file_type;
                             L : Standard_Irreducible_Component_List;
                             x : Standard_Complex_Vectors.Vector;
                             tol : double_float ) return natural32 is

    tmp : Standard_Irreducible_Component_List := L;

  begin
    for i in 1..Length_Of(L) loop
      put(file,"Membership test at component "); put(file,i,1);
      put_line(file," : "); 
      if Homotopy_Membership_Test(file,Head_Of(tmp),x,tol)
       then return i;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end Homotopy_Filter;

end Irreducible_Component_Lists;
