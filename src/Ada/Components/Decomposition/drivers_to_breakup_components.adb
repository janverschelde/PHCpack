with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Multprec_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;        use Standard_Natural_Vectors_io;
with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Standard_Natural_VecVecs;
with Standard_Integer_VecVecs;
with Standard_Complex_VecVecs;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;    use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems_io;   use Multprec_Complex_Poly_Systems_io;
with Standard_Complex_Solutions_io;      use Standard_Complex_Solutions_io;
with Standard_Breakup_Components;        use Standard_Breakup_Components;
with Multprec_Breakup_Components;        use Multprec_Breakup_Components;
with Witness_Sets,Witness_Sets_io;       use Witness_Sets,Witness_Sets_io;
with Sample_Points;                      use Sample_Points;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with Sampling_Machine;
with Monodromy_Group_Actions;            use Monodromy_Group_Actions;
with Monodromy_Group_Actions_io;         use Monodromy_Group_Actions_io;

package body Drivers_to_Breakup_Components is

-- AUXILIARY DATA STRUCTURE FOR MONODROMY BREAKUP :

  type Boolean_Array is array ( integer32 range <> ) of boolean;

-- BREAK UP COMPONENTS WITH INTERPOLATION FILTERS :

  procedure Standard_Breakup_Components
               ( file : in file_type; level,itp : in natural32;
                 skewproj : in boolean;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List ) is

    n : constant integer32 := p'length;
    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
    hyp : Standard_Complex_VecVecs.VecVec(p'range);
    intpols : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    subspaces : Standard_Complex_Poly_Systems.Array_of_Poly_Sys(1..len);
    pivots : Standard_Integer_VecVecs.VecVec(1..len);
    basepts : Standard_Complex_VecVecs.Array_of_VecVecs(1..len);
    basecard : Standard_Natural_Vectors.Vector(1..len) := (1..len => 0);
    deg,sumdeg : natural32;

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

  begin
    for i in hyp'range loop
      hyp(i) := new Standard_Complex_Vectors.Vector'(Random_Vector(0,n));
    end loop;
    case itp is
      when 1 => intpols := new Poly_Sys'
                  (Massive_Interpolate(file,p,sols,hyp,level));
      when 2 => intpols := new Poly_Sys'
                  (Incremental_Interpolate(file,p,sols,hyp,level));
      when 3 => Dynamic_Interpolate
                  (file,p,level,sols,hyp,skewproj,subspaces,pivots,
                   basepts,basecard,intpols);
      when others => null;
    end case;
    if intpols /= null then
      sumdeg := 0;
      new_line(file);
      put_line(file,"THE INTERPOLATING POLYNOMIALS : ");
      for i in intpols'range loop
        if not (intpols(i) = Null_Poly) then
          put(file,"Equation for component "); 
          put(file,i,1); put(file," of degree ");
          deg := natural32(Degree(intpols(i)));
          put(file,deg,1);
          sumdeg := sumdeg + deg + basecard(i);
          put_line(file," : ");
          put_line(file,intpols(i));
          if subspaces(i) /= null then
            put_line(file,"Container subspace : ");
            put_line(file,subspaces(i).all);
            put(file,"with pivots : ");
            put(file,pivots(i)); new_line(file);
          end if;      
        end if;
      end loop;
      new_line(file);
      put(file,"Sum of degrees of polynomials : ");
      put(file,sumdeg,1); put_line(file,".");
      put(file,"Number of central points : ");
      put(file,basecard); new_line(file);
    end if;
  end Standard_Breakup_Components;

  procedure Multprec_Breakup_Components
               ( file : in file_type; level,size,itp : in natural32;
                 skewproj : in boolean;
                 p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                 sols : in Standard_Complex_Solutions.Solution_List ) is

    n : constant integer32 := p'length;
    len : constant integer32
        := integer32(Standard_Complex_Solutions.Length_Of(sols));
    hyp : Multprec_Complex_VecVecs.VecVec(p'range);
    intpols : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    subspaces : Multprec_Complex_Poly_Systems.Array_of_Poly_Sys(1..len);
    pivots : Standard_Integer_VecVecs.VecVec(1..len);
    basepts : Multprec_Complex_VecVecs.Array_of_VecVecs(1..len);
    basecard : Standard_Natural_Vectors.Vector(1..len) := (1..len => 0);
    cnt : natural32 := 0;
    deg,sumdeg : natural32;

    use Multprec_Complex_Vectors;
    use Multprec_Complex_VecVecs;
    use Multprec_Complex_Polynomials;
    use Multprec_Complex_Poly_Systems;

  begin
    for i in hyp'range loop
      hyp(i) := new Multprec_Complex_Vectors.Vector'
                      (Create(Random_Vector(0,n)));
    end loop;
    new_line(file);
    put_line(file,"THE SYSTEM WITH MULTI-PRECISION COEFFICIENTS : ");
    put(file,mp);
    case itp is
      when 1 => intpols := new Multprec_Complex_Poly_Systems.Poly_Sys'
                  (Massive_Interpolate(file,p,mp,sols,hyp,level,size));
      when 2 => intpols := new Multprec_Complex_Poly_Systems.Poly_Sys'
                  (Incremental_Interpolate(file,p,mp,sols,hyp,level,size));
      when 3 => Dynamic_Interpolate
                  (file,p,mp,sols,level,size,hyp,skewproj,subspaces,pivots,
                   basepts,basecard,intpols);
      when others => null;
    end case;
    if intpols /= null then
      sumdeg := 0;
      new_line(file);
      put_line(file,"THE INTERPOLATING POLYNOMIALS : ");
      for i in intpols'range loop
        if not (intpols(i) = Null_Poly) then
          cnt := cnt + 1;
          put(file,"Equation for component "); 
          put(file,cnt,1); put(file," of degree ");
          deg := natural32(Degree(intpols(i)));
          sumdeg := sumdeg + deg;
          put(file,Degree(intpols(i)),1);
          if basepts(i) /= null then
            put(file," + "); put(file,basecard(i),1);
            put(file," (#central points)");
            sumdeg := sumdeg + basecard(i);
          end if;
          put_line(file," : ");
          put_line(file,intpols(i));
          if subspaces(i) /= null then
            put_line(file,"Container subspace : ");
            put_line(file,subspaces(i).all);
            put(file,"with pivots : ");
            put(file,pivots(i)); new_line(file);
          end if;      
        end if;
      end loop;
      new_line(file);
      put(file,"Sum of degrees of polynomials : ");
      put(file,sumdeg,1); new_line(file);
      put(file,"Number of central points : ");
      put(file,basecard); new_line(file);
    end if;
  end Multprec_Breakup_Components;

  procedure Breakup_Menu
              ( itp,size : out natural32; skewproj : out boolean ) is

    ans : character;
    deci : natural32 := 0;

  begin
    new_line;
    put_line("MENU to break up components with interpolation filters : ");
    put_line("  1. massive interpolate using standard numbers;");
    put_line("  2.                     using multi-precision numbers;");
    put_line("  3. incremental interpolate using standard numbers;");
    put_line("  4.                         using multi-precision numbers;");
    put_line("  5. dynamic interpolate with subspace restriction "
                   & "using standard numbers;");
    put_line("  6.                                               "
                   & "using multi-precision;");
    put_line("  7. dynamic interpolate with central projections "
                   & "using standard numbers;");
    put_line("  8.                                              "
                   & "using multi-precision.");
    put("Type 1,2,3,4,5,6,7, or 8 to select : ");
    Ask_Alternative(ans,"12345678");
    case ans is
      when '1' | '2' => itp := 1; skewproj := false;
      when '3' | '4' => itp := 2; skewproj := false;
      when '5' | '6' => itp := 3; skewproj := false;
      when '7' | '8' => itp := 3; skewproj := true;
      when others    => itp := 0;
    end case;
    if ans = '2' or ans = '4' or ans = '6' or ans = '8' then
      new_line;
      put("Give the number of decimal places : "); get(deci);
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    else
      size := 0;
    end if;
  end Breakup_Menu;

  procedure Breakup_with_Interpolation_Filters is

    file : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    lmp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    level,size,itp : natural32;
    skewproj : boolean;

    use Standard_Complex_Solutions;

  begin
    Standard_Read_Embedding(lp,sols,level);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Breakup_Menu(itp,size,skewproj);
    if size > 0 then
      Get_Multprec_System(lp.all,lmp,size,level);
      put_line(file,"The system in multi-precision format : ");
      put_line(file,lmp.all);
    end if;
    new_line;
    put_line("See the output file for results...");
    new_line;
    put_line(file,"THE EMBEDDED SYSTEM : ");
    put_line(file,lp.all);
    put_line(file,"THE GENERIC POINTS : ");
    put(file,Length_Of(sols),natural32(Head_Of(sols).n),sols);
    if size = 0 then
      Standard_Breakup_Components(file,level,itp,skewproj,lp.all,sols);
    else
      Multprec_Breakup_Components
        (file,level,size,itp,skewproj,lp.all,lmp.all,sols);
    end if;
  end Breakup_with_Interpolation_Filters;

-- BREAK UP COMPONENTS USING MONODROMY GROUP ACTIONS :

  procedure Correlate
               ( file : in file_type;
                 sps : in Standard_Sample_List; tol : in double_float;
                 ic : in Irreducible_Components;
                 rel : out Standard_Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Two samples in sps are related if there are paths connecting them.
  --   If on return rel(i) = j, then (i,j) are related.
  -- REQUIRED :
  --   The sampling machine is initialized and tuned;
  --   rel'range = 1..Length_Of(sps).

    tmp : Standard_Sample_List := sps;
    spt : Standard_Sample;
    found : Boolean_Array(1..integer32(Length_Of(sps)));

  begin
    found := (found'range => false);
    for i in 1..integer32(Length_Of(sps)) loop
      if not Empty(ic,i) and not found(i) then
        put(file,"Sampling from point "); put(file,i,1);
        spt := Head_Of(tmp);
        declare
          newspt : Standard_Sample;
          isin : natural32;
          s,s_last : Standard_Sample_List;
        begin
          Sample(spt,newspt);
          Membership_Test(sps,Sample_Point(newspt).v,tol,isin,s,s_last);
          put(file," moves to "); put(file,isin,1);
          if isin = 0 then
            put(file,"  Oops, bug fixed");
            isin := natural32(i);
          end if;
          new_line(file);
          Deep_Clear(s);
          Deep_Clear(newspt);
          rel(i) := isin;
          found(integer32(isin)) := true;
        end;
      else
        rel(i) := natural32(i);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Correlate;

  function Active_Samples
               ( ic : Irreducible_Components; sps : Standard_Sample_List )
               return Standard_Sample_List is

  -- DESCRIPTION :
  --   A sample is active if in the list of irreducible components
  --   it represents its component.
  --   This function returns the list of active samples as a sublist of sps.

  -- WARNING :
  --   The sample list on input must be the master starting set of generic
  --   points, i.e.: Length_Of(sps) = Sum_of_Degrees(ic).

    res,res_last : Standard_Sample_List;
    tmp : Standard_Sample_List := sps; 

  begin
    for i in 1..integer32(Length_Of(sps)) loop
      if not Empty(ic,i)
       then Append(res,res_last,Head_Of(tmp));
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Active_Samples;

  function Active_Samples
               ( sps : Standard_Sample_List;
                 ind : Standard_Natural_Vectors.Vector )
               return Standard_Sample_List is

  -- DESCRIPTION :
  --   Returns the list of samples selected from sps with the
  --   corresponding indices in ind.
  --   This function returns the list of active samples as a sublist of sps.

  -- WARNING :
  --   The sample list on input must be the master starting set of generic
  --   points, i.e.: Length_Of(sps) = Sum_of_Degrees(ic).

    res,res_last : Standard_Sample_List;
    tmp : Standard_Sample_List := sps; 
    ptr : integer32 := ind'first;

  begin
    for i in 1..Length_Of(sps) loop
      if ind(ptr) = i then
        Append(res,res_last,Head_Of(tmp));
        ptr := ptr + 1;
      end if;
      exit when (ptr > ind'last);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Active_Samples;

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

  procedure Monodromy_Decomposition
                   ( file : in file_type;
                     p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     sols : in Standard_Complex_Solutions.Solution_List;
                     dim : in natural32 ) is

  -- DESCRIPTION :
  --   This does some runs on the sampling membership test.

    tol : double_float := 1.0E-8;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,hyp);
    newhyp : Standard_Complex_VecVecs.VecVec(1..integer32(dim));
   --        := Random_Hyperplanes(dim,p'last);
    newsps,newsps_last : Standard_Sample_List;
    rel,lab : Standard_Natural_Vectors.Vector(1..integer32(Length_Of(sps)));
    n1,n2,cnt,threshold,nit,jump : natural32 := 0;
    timer : Timing_Widget;
    ic : Irreducible_Components := Create(rel'last);

  begin
    new_line;
    put("Give threshold for number of stable iterations : "); get(threshold);
   -- put("Give the sample reduction factor : "); get(jump);
    jump := 1;
    cnt := 0;
    tstart(timer);
    for i in lab'range loop
      lab(i) := natural32(i);
    end loop;
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    n1 := Cardinality(ic);
    new_line;
    put_line("See the output file for results...");
    new_line;
    put(file,"The decomposition : "); put(file,ic);
    nit := 0;
   -- Sample(sps,newhyp,newsps,newsps_last);
    loop
     -- Correlate(file,sps,tol,ic,rel);
      newhyp := Random_Hyperplanes(dim,natural32(p'last));
      Sample(sps,newhyp,newsps,newsps_last);
      if jump > 1 then
        declare
        -- rep : constant Standard_Natural_Vectors.Vector 
        --     := Nonempty_Sets(ic);
          rep : constant Standard_Natural_Vectors.Vector
              := Representatives(ic,integer32(nit),integer32(jump));
        -- act : Standard_Sample_List := Active_Samples(ic,sps);
          act : Standard_Sample_List := Active_Samples(sps,rep);
        begin
          Generate_Monodromy_Action(file,rep,act,newsps,newhyp,tol,rel);
        end;
      else
        Generate_Monodromy_Action(file,lab,sps,newsps,newhyp,tol,rel);
      end if;
      Deep_Clear(newsps);
      Act(ic,rel);
      n2 := Cardinality(ic);
      nit := nit+1;
      put(file,"The decomposition : "); put(file,ic); 
      put(file,"Degrees : "); put(file,Degrees(ic)); new_line(file);
      put(file,"cardinality ");
      put(file,n1,1); put(file," -> "); put(file,n2,1);
      if n1 = n2 then
        put(file," is stable with cnt = ");
        cnt := cnt+1;
        put(file,cnt,1); put_line(file,".");
      else
        put_line(file," drops.");
        cnt := 0;
      end if;
      exit when (cnt = threshold);
      n1 := n2;
    end loop;
    Sampling_Machine.Clear;
    tstop(timer);
    put(file,"Number of iterations : "); put(file,nit,1); new_line(file);
    put(file,"Stabilizing threshold : "); put(file,threshold,1);
    new_line(file);
    put_line(file,"The decomposition : "); put_labels(file,ic);
    new_line(file);
   -- put(file,"Sample reduction factor : "); put(file,jump,1); new_line(file);
   -- new_line(file);
    print_times(file,timer,"Decomposition using Monodromy Group Actions");
  end Monodromy_Decomposition;

  procedure Breakup_with_Monodromy_Group_Actions is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List; 
    dim : natural32;
    file : file_type;

  begin
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    Monodromy_Decomposition(file,lp.all,sols,dim);
  end Breakup_with_Monodromy_Group_Actions;

end Drivers_to_Breakup_Components;
