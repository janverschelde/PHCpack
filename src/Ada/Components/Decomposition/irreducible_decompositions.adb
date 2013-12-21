with unchecked_deallocation;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Witness_Sets;                       use Witness_Sets;
with Sampling_Machine;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Irreducible_Components;             use Irreducible_Components;
with Filtered_Points;

package body Irreducible_Decompositions is

  type Array_of_Solution_Lists is
    array ( integer32 range <> ) of Standard_Complex_Solutions.Solution_List;

  type Array_of_Standard_Component_Lists is
    array ( integer32 range <> ) of Standard_Irreducible_Component_List;

  type Array_of_Multprec_Component_Lists is
    array ( integer32 range <> ) of Multprec_Irreducible_Component_List;

  type Standard_Irreducible_Decomposition_Rep ( k : integer32 ) is record
    ep : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..k);
    gp : Array_of_Solution_Lists(0..k);
    cp,cp_last : Array_of_Standard_Component_Lists(0..k);
  end record;

  type Multprec_Irreducible_Decomposition_Rep ( k : integer32 ) is record
    ep : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(0..k);
    gp : Array_of_Solution_Lists(0..k);
    cp,cp_last : Array_of_Multprec_Component_Lists(0..k);
    op : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
  end record;

-- AUXILIARY TO BREAKUP :

  procedure Create_Filter_Data
                ( dc : in Standard_Irreducible_Decomposition;
                  i : in integer32; fp,fp_last : in out List ) is

  -- DESCRIPTION :
  --   For every i dimensional component, a vector will be added to fp
  --   that has as its i-th entry the degree of the co,ponent.

    tmp : Standard_Irreducible_Component_List := dc.cp(i);

  begin
    while not Is_Null(tmp) loop
      declare
        c : constant Standard_Irreducible_Component := Head_Of(tmp);
        v : Standard_Integer_Vectors.Vector(0..dc.k+1);
      begin
        v := Filtered_Points.Create(dc.k+1,i,integer32(Degree(c)));
        Append(fp,fp_last,v);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Create_Filter_Data;

  procedure Create_Filter_Data
                ( dc : in Multprec_Irreducible_Decomposition;
                  i : in integer32; fp,fp_last : in out List ) is

  -- DESCRIPTION :
  --   For every i dimensional component, a vector will be added to fp
  --   that has as its i-th entry the degree of the co,ponent.

    tmp : Multprec_Irreducible_Component_List := dc.cp(i);

  begin
    while not Is_Null(tmp) loop
      declare
        c : constant Multprec_Irreducible_Component := Head_Of(tmp);
        v : Standard_Integer_Vectors.Vector(0..dc.k+1);
      begin
        v := Filtered_Points.Create(dc.k+1,i,integer32(Degree(c)));
        Append(fp,fp_last,v);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Create_Filter_Data;

  procedure Update_Filter_Data ( fp : in out List; i,ind : integer32 ) is

  -- DESCRIPTION :
  --   One point is found to belong to the component that corresponds to
  --   the element at place ind of the list fp.

    lv : Standard_Integer_Vectors.Link_to_Vector;
    tmpfp : List := fp;

  begin
    for j in 1..(ind-1) loop
      exit when Is_Null(tmpfp);
      tmpfp := Tail_Of(tmpfp);
    end loop;
    if not Is_Null(tmpfp) then
      lv := Head_Of(tmpfp);
      lv(i) := lv(i)+1;
      Set_Head(tmpfp,lv);
    end if;
  end Update_Filter_Data;

  procedure Isolated_Filter_Data
              ( k : in integer32; fp,fp_last : in out List ) is

  -- DESCRIPTION :
  --   Processes the finding of an isolated point.

    lv : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if not Is_Null(fp_last) then
      lv := Head_Of(fp_last);
      if lv(k+1) /= 0 then
        Append(fp,fp_last,Filtered_Points.Create(k+1,0,0));
        lv := Head_Of(fp_last);
      end if;
      lv(0) := lv(0) + 1;
      Set_Head(fp_last,lv);
    end if;
  end Isolated_Filter_Data;

-- CREATORS :

  function Create ( k : integer32 )
                  return Standard_Irreducible_Decomposition is

    res : Standard_Irreducible_Decomposition;
    res_rep : Standard_Irreducible_Decomposition_Rep(k);

  begin
    res := new Standard_Irreducible_Decomposition_Rep'(res_rep);
    return res;
  end Create;

  function Create ( p : Standard_Complex_Poly_Systems.Array_of_Poly_Sys )
                  return Standard_Irreducible_Decomposition is

    res : Standard_Irreducible_Decomposition;
    res_rep : Standard_Irreducible_Decomposition_Rep(p'last);

  begin
    res_rep.ep := p;
    res := new Standard_Irreducible_Decomposition_Rep'(res_rep);
    return res;
  end Create;

  function Create ( p : Standard_Complex_Poly_Systems.Array_of_Poly_Sys )
                  return Multprec_Irreducible_Decomposition is

    res : Multprec_Irreducible_Decomposition;
    res_rep : Multprec_Irreducible_Decomposition_Rep(p'last);

  begin
    res_rep.ep := p;
    res := new Multprec_Irreducible_Decomposition_Rep'(res_rep);
    return res;
  end Create;

  function Create ( dc : Standard_Irreducible_Decomposition )
                  return Multprec_Irreducible_Decomposition is

    res : Multprec_Irreducible_Decomposition;
    res_rep : Multprec_Irreducible_Decomposition_Rep(dc.k);

  begin
    res_rep.ep := dc.ep;
    res_rep.gp := dc.gp;
    res := new Multprec_Irreducible_Decomposition_Rep'(res_rep);
    return res;
  end Create;

  procedure Add_Original
              ( dc : in out Multprec_Irreducible_Decomposition;
                p : in Multprec_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    dc.op := p;
  end Add_Original;

  procedure Add_Embedding
              ( dc : in out Standard_Irreducible_Decomposition;
                i : in integer32;
                ep : in Standard_Complex_Poly_Systems.Link_to_Poly_Sys ) is
  begin
    if dc /= null then
      if i <= dc.k
       then dc.ep(i) := ep;
      end if;
    end if;
  end Add_Embedding;

  procedure Add_Generic_Points
              ( dc : in out Standard_Irreducible_Decomposition;
                i : in integer32;
                gp : in Standard_Complex_Solutions.Solution_List ) is
  begin
    if dc /= null then
      if i <= dc.k
       then dc.gp(i) := gp;
      end if;
    end if;
  end Add_Generic_Points;

  procedure Add_Generic_Points
              ( dc : in out Multprec_Irreducible_Decomposition;
                i : in integer32;
                gp : in Standard_Complex_Solutions.Solution_List ) is
  begin
    if dc /= null then
      if i <= dc.k
       then dc.gp(i) := gp;
      end if;
    end if;
  end Add_Generic_Points;

  procedure Breakup ( file : in file_type; full_output : in boolean;
                      dc : in out Standard_Irreducible_Decomposition;
                      i : in integer32; method : in natural32;
                      stoptol,membtol : in double_float;
                      fp,fp_last : in out List ) is

    sli : Standard_Complex_VecVecs.VecVec(1..i);
    sps : Standard_Sample_List;

  begin
    if not Standard_Complex_Solutions.Is_Null(dc.gp(i)) then
      sli := Slices(dc.ep(i).all,natural32(i));
      sps := Create(dc.gp(i),sli);
      Sampling_Machine.Initialize(dc.ep(i).all);
      case method is
        when 0 =>
          Standard_Massive_Interpolate
            (file,full_output,sps,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when 1 => 
          Standard_Incremental_Interpolate
            (file,full_output,sps,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when 2 => 
          Standard_Incremental_Interpolate_with_Span
            (file,full_output,sps,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when 3 => 
          Standard_Incremental_Central_Interpolate
            (file,full_output,sps,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when others => null;
      end case;
      Sampling_Machine.Clear;
      Create_Filter_Data(dc,i,fp,fp_last);
    end if;
  end Breakup;

  procedure Monodromy_Breakup 
              ( file : in file_type;
                 dc : in out Standard_Irreducible_Decomposition;
                 i : in integer32; threshold : in natural32;
                 tol : in double_float;
                 fp,fp_last : in out List ) is

    sli : Standard_Complex_VecVecs.VecVec(1..i);
    sps : Standard_Sample_List;
    nit : natural32;

  begin
    if not Standard_Complex_Solutions.Is_Null(dc.gp(i)) then
      sli := Slices(dc.ep(i).all,natural32(i));
      sps := Create(dc.gp(i),sli);
      Sampling_Machine.Initialize(dc.ep(i).all);
      Monodromy_Breakup(file,sps,threshold,tol,dc.cp(i),dc.cp_last(i),nit);
      Sampling_Machine.Clear;
      Create_Filter_Data(dc,i,fp,fp_last);
    end if;
  end Monodromy_Breakup;

  procedure Breakup ( file : in file_type; full_output : in boolean;
                      dc : in out Multprec_Irreducible_Decomposition;
                      i : in integer32; method,size : in natural32;
                      stoptol,membtol : in double_float;
                      fp,fp_last : in out List ) is

    sli : Standard_Complex_VecVecs.VecVec(1..i);
    sps : Standard_Sample_List;

  begin
    if not Standard_Complex_Solutions.Is_Null(dc.gp(i)) then
      sli := Slices(dc.ep(i).all,natural32(i));
      sps := Create(dc.gp(i),sli);
      Sampling_Machine.Initialize(dc.ep(i).all,dc.op.all,i,size);
      case method is
        when 0 =>
          Multprec_Massive_Interpolate
            (file,full_output,sps,size,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when 1 => 
          Multprec_Incremental_Interpolate
            (file,full_output,sps,size,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when 2 => 
          Multprec_Incremental_Interpolate_with_Span
            (file,full_output,sps,size,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when 3 => 
          Multprec_Incremental_Central_Interpolate
            (file,full_output,sps,size,stoptol,membtol,dc.cp(i),dc.cp_last(i));
        when others => null;
      end case;
      Sampling_Machine.Clear;
      Create_Filter_Data(dc,i,fp,fp_last);
    end if;
  end Breakup;

  procedure Filter ( file : in file_type;
                     dc : in out Standard_Irreducible_Decomposition;
                     i : in integer32; membtol : in double_float;
                     fp,fp_last : in out List; junkcnt : out natural32 ) is

    use Standard_Complex_Solutions;
    tmp,filtgp,filtgp_last : Solution_List;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    ind : integer32;

  begin
    if not Is_Null(dc.gp(i)) then
      put(file,"Classyfying points at dimension "); put(file,i,1);
      put_line(file," with filters of higher dimension :");
      tmp := dc.gp(i);
      for j in 1..Length_Of(dc.gp(i)) loop     -- test solution j
        ls := Head_Of(tmp);
        ind := 0;
        for l in reverse (i+1)..dc.k loop      -- test on component l
          ind := integer32(On_Component(file,dc.cp(l),ls.v,membtol));
          if ind > 0 then
            put(file,"Point "); put(file,j,1);
            put(file," lies on component "); put(file,ind,1);
            put(file," of dimension "); put(file,l,1);
            put_line(file,".");
            Update_Filter_Data(fp,i,ind);
            cnt := cnt+1;
            exit;
          end if;
        end loop;
        if ind = 0 then
          Append(filtgp,filtgp_last,ls.all);
          if i = 0
           then Isolated_Filter_Data(dc.k,fp,fp_last);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      Clear(dc.gp(i));
      dc.gp(i) := filtgp;
    end if;
    junkcnt := cnt;
  end Filter;

  procedure Filter ( file : in file_type;
                     dc : in out Multprec_Irreducible_Decomposition;
                     i : in integer32; membtol : in double_float;
                     fp,fp_last : in out List; junkcnt : out natural32 ) is

    use Standard_Complex_Solutions;
    tmp,filtgp,filtgp_last : Solution_List;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    ind : integer32;

  begin
    if not Is_Null(dc.gp(i)) then
      put(file,"Classyfying points at dimension "); put(file,i,1);
      put_line(file," with filters of higher dimension :");
      tmp := dc.gp(i);
      for j in 1..Length_Of(dc.gp(i)) loop     -- test solution j
        ls := Head_Of(tmp);
        ind := 0;
        for l in reverse (i+1)..dc.k loop      -- test on component l
          ind := integer32(On_Component(file,dc.cp(l),ls.v,membtol));
          if ind > 0 then
            put(file,"Point "); put(file,j,1);
            put(file," lies on component "); put(file,ind,1);
            put(file," of dimension "); put(file,l,1);
            put_line(file,".");
            Update_Filter_Data(fp,i,ind);
            cnt := cnt+1;
            exit;
          end if;
        end loop;
        if ind = 0 then
          Append(filtgp,filtgp_last,ls.all);
          if i = 0
           then Isolated_Filter_Data(dc.k,fp,fp_last);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      Clear(dc.gp(i));
      dc.gp(i) := filtgp;
    end if;
    junkcnt := cnt;
  end Filter;

  procedure Homotopy_Filter
              ( file : in file_type;
                dc : in out Standard_Irreducible_Decomposition;
                i : in integer32; membtol : in double_float;
                fp,fp_last : in out List; junkcnt : out natural32 ) is

    use Standard_Complex_Solutions;
    tmp,filtgp,filtgp_last : Solution_List;
    ls : Link_to_Solution;
    cnt : natural32 := 0;
    ind : integer32;

  begin
    if not Is_Null(dc.gp(i)) then
      put(file,"Classyfying points at dimension "); put(file,i,1);
      put_line(file," with homotopy membership test :");
      tmp := dc.gp(i);
      for j in 1..Length_Of(dc.gp(i)) loop     -- test solution j
        ls := Head_Of(tmp);
        ind := 0;
        for l in reverse (i+1)..dc.k loop      -- test on component l
          Sampling_Machine.Initialize(dc.ep(l).all);
          ind := integer32(Homotopy_Filter(file,dc.cp(l),ls.v,membtol));
          Sampling_Machine.Clear;
          if ind > 0 then
            put(file,"Point "); put(file,j,1);
            put(file," lies on component "); put(file,ind,1);
            put(file," of dimension "); put(file,l,1);
            put_line(file,".");
            Update_Filter_Data(fp,i,ind);
            cnt := cnt+1;
            exit;
          end if;
        end loop;
        if ind = 0 then
          Append(filtgp,filtgp_last,ls.all);
          if i = 0
           then Isolated_Filter_Data(dc.k,fp,fp_last);
          end if;
        end if;
        tmp := Tail_Of(tmp);
      end loop;
      Clear(dc.gp(i));
      dc.gp(i) := filtgp;
    end if;
    junkcnt := cnt;
  end Homotopy_Filter;

-- SELECTORS :

  function Top_Dimension ( dc : Standard_Irreducible_Decomposition )
                         return integer32 is
  begin
    if dc = null
     then return 0;
     else return dc.k;
    end if;
  end Top_Dimension;

  function Top_Dimension ( dc : Multprec_Irreducible_Decomposition )
                         return integer32 is
  begin
    if dc = null
     then return 0;
     else return dc.k;
    end if;
  end Top_Dimension;

  function Embedding ( dc : Standard_Irreducible_Decomposition;
                       i : integer32 )
                     return Standard_Complex_Poly_Systems.Link_to_Poly_Sys is

    res : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if dc /= null then
      if i <= dc.k
       then res := dc.ep(i);
      end if;
    end if;
    return res;
  end Embedding;

  function Original ( dc : Multprec_Irreducible_Decomposition )
                    return Multprec_Complex_Poly_Systems.Link_to_Poly_Sys is

    res : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;

  begin
    if dc /= null
     then res := dc.op;
    end if;
    return res;
  end Original;

  function Generic_Points ( dc : Standard_Irreducible_Decomposition;
                            i : integer32 )
                          return Standard_Complex_Solutions.Solution_List is

    res : Standard_Complex_Solutions.Solution_List;

  begin
    if dc /= null then
      if i <= dc.k
       then res := dc.gp(i);
      end if;
    end if;
    return res;
  end Generic_Points;

  function Generic_Points ( dc : Multprec_Irreducible_Decomposition;
                            i : integer32 )
                          return Standard_Complex_Solutions.Solution_List is

    res : Standard_Complex_Solutions.Solution_List;

  begin
    if dc /= null then
      if i <= dc.k
       then res := dc.gp(i);
      end if;
    end if;
    return res;
  end Generic_Points;

  function Components ( dc : Standard_Irreducible_Decomposition;
                        i : integer32 )
                      return Standard_Irreducible_Component_List is
  begin
    return dc.cp(i);
  end Components;

  function Components ( dc : Multprec_Irreducible_Decomposition;
                        i : integer32 )
                      return Multprec_Irreducible_Component_List is
  begin
    return dc.cp(i);
  end Components;

-- DESTRUCTORS :

  procedure free is
    new unchecked_deallocation(Standard_Irreducible_Decomposition_Rep,
                               Standard_Irreducible_Decomposition);

  procedure free is
    new unchecked_deallocation(Multprec_Irreducible_Decomposition_Rep,
                               Multprec_Irreducible_Decomposition);

  procedure Clear ( dc : in out Standard_Irreducible_Decomposition ) is
  begin
    if dc /= null then
      for i in 0..dc.k loop
        Standard_Complex_Poly_Systems.Clear(dc.ep(i));
        Standard_Complex_Solutions.Clear(dc.gp(i));
        Clear(dc.cp(i));
      end loop;
      free(dc);
    end if;
  end Clear;

  procedure Clear ( dc : in out Multprec_Irreducible_Decomposition ) is
  begin
    if dc /= null then
      for i in 0..dc.k loop
        Standard_Complex_Poly_Systems.Clear(dc.ep(i));
        Standard_Complex_Solutions.Clear(dc.gp(i));
        Clear(dc.cp(i));
      end loop;
      free(dc);
    end if;
  end Clear;

end Irreducible_Decompositions;
