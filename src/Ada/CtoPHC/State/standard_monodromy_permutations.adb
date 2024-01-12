with unchecked_deallocation;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Rectangular_Sample_Grids;
with Monodromy_Component_Breakup;
with PHCpack_Operations;
with Monodromy_Partitions;              use Monodromy_Partitions;
with Standard_Sampling_Operations;

with text_io;                           use text_io;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Solutions_io;     use Standard_Complex_Solutions_io;

package body Standard_Monodromy_Permutations is

-- INTERNAL DATA :

  type Link_to_Array_of_Solution_Lists is access Array_of_Solution_Lists;

  ind : integer32;
  dim : natural32;
  grid : Link_to_Array_of_Solution_Lists;
  deco : Link_to_VecVec;
  trace_grid : Array_of_Standard_Sample_Lists(0..2);
  trace_grid_maximal_error : double_float := 1.0;
    -- maximal error on samples in trace grid
  trace_grid_minimal_distance : double_float := 0.0;
    -- minimal distance between samples in trace grid

  tol : constant double_float := 1.0E-8;

-- OPERATIONS :

  procedure Write_Grid is
  begin
    for i in 0..ind loop
      put("Solutions at grid("); put(i,1);
      put_line(") :");
      put(standard_output,Length_Of(grid(i)),
          natural32(Head_Of(grid(i)).n),grid(i));
    end loop;
  end Write_Grid;

  procedure Assign_Labels ( sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   The multiplicity is the position of the solution in the list.

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      ls.m := integer32(i);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Assign_Labels;

  procedure Initialize ( n,d,k : in integer32 ) is

    n2 : constant integer32 := n+2;

  begin
    ind := -1;
    dim := natural32(k);
    grid := new Array_of_Solution_Lists(0..n2);
    deco := Monodromy_Partitions.Init_Factors(natural32(d));
  end Initialize;

  procedure On_Slice ( sols : in Solution_List;
                       v : in Standard_Complex_Vectors.Vector ) is 

  -- DESCRIPTION :
  --   This routine prints the result of the evaluation of every
  --   solution at the hyperplane with coefficients in v.
  --   Good for checking purposes.
 
    tmp : Solution_List := sols;
    ls : Link_to_Solution;
    res : Complex_Number := Create(0.0);

  begin
    for i in 1..Length_Of(sols) loop
      ls := Head_Of(tmp);
      res := v(0);
      for j in v'range loop
        res := res + ls.v(j)*v(j);
      end loop;
      put("residual on slice : "); put(res); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end On_Slice;

  function Retrieve ( sols : Solution_List;
                      m : integer32 ) return Link_to_Solution is

  -- DESCRIPTION :
  --   Retrieves the solution with the m equal to the given m.

    tmp : Solution_List := sols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      if ls.m = m
       then return ls;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return ls;
  end Retrieve;

  function Sort ( sols : Solution_List ) return Solution_List is

  -- DESCRIPTION :
  --   Returns the same solution list, sorted in increasing order
  --   of the m field.  This is the needed for the trace grid.

    res,res_last : Solution_List;

  begin
    for i in 1..Length_Of(sols) loop
      Append(res,res_last,Retrieve(sols,i));
    end loop;
    return res;
  end Sort;

  procedure On_Slices_Test ( sols : in Solution_List; i : integer32 ) is

  -- DESCRIPTION :
  --   For hyperplane_sections /= null, it makes sense to test whether
  --   the i-th solution list lies on slice i-2, for i >= 3

    v : constant Standard_Complex_VecVecs.Link_to_VecVec
      := Standard_Sampling_Operations.Retrieve_Slices(i-2);

    use Standard_Complex_VecVecs;

  begin
    if v /= null then
      put("Executing on_slice test for i = "); put(i,1);
      put_line(" ...");
      On_Slice(sols,v(1).all);
    else
      put("slice "); put(i-2,1); put_line(" is empty...");
    end if;
  end On_Slices_Test;

  procedure Store ( sols : in Solution_List ) is

    sli : Standard_Complex_VecVecs.Link_to_VecVec;

  begin
   -- put("Storing with ind = "); put(ind,1); new_line;
    ind := ind + 1;
    Copy(sols,grid(ind));
    if ind = 0
     then Assign_Labels(grid(ind));
    end if;
   -- put(standard_output,Length_Of(sols),Head_Of(sols).n,sols);
    if ind < 3 then
     -- put("storing trace grid for slice "); put(ind,1); new_line;
      sli := Standard_Sampling_Operations.Retrieve_Start_Slices;
     -- put_line("the retrieved slices :"); put(sli);
     -- On_Slice(sols,sli(1).all);
      declare
        cp_sli : Standard_Complex_VecVecs.VecVec(sli'range);
      begin
        for i in sli'range loop
          cp_sli(i) := new Standard_Complex_Vectors.Vector'(sli(i).all);
        end loop;
        if ind = 0 then
          trace_grid(ind) := Create(sols,cp_sli);
        else
          trace_grid(ind) := Create(Sort(sols),cp_sli);
          if ind = 2 then
            trace_grid_maximal_error
              := Rectangular_Sample_Grids.Maximal_Error(trace_grid);
            if not stay_silent then
              if PHCpack_Operations.Is_File_Defined then
                put(PHCpack_Operations.output_file,
                    "maximal error on trace grid : ");
                put(PHCpack_Operations.output_file,
                    trace_grid_maximal_error,3);
                new_line(PHCpack_Operations.output_file);
              end if;
            end if;
            trace_grid_minimal_distance
              := Rectangular_Sample_Grids.Minimal_Distance(trace_grid);
            if not stay_silent then
              if PHCpack_Operations.Is_File_Defined then
                put(PHCpack_Operations.output_file,
                    "minimial distance on trace grid : ");
                put(PHCpack_Operations.output_file,
                    trace_grid_minimal_distance,3);
                new_line(PHCpack_Operations.output_file);
              end if;
            end if;
          end if;
        end if;
      end;
    end if;
   -- if ind > 2 then
   --   On_Slices_Test(sols,ind);
   -- end if;
   -- Write_Grid;
  end Store;

  function Retrieve ( k : integer32 ) return Solution_List is

    res : Solution_List;

  begin
    if ((grid /= null) and then (k <= grid'last))
     then res := grid(k);
    end if;
    return res;
  end Retrieve;

  procedure Trace_Grid_Diagnostics
              ( maximal_error,minimal_distance : out double_float ) is
  begin
    maximal_error := trace_grid_maximal_error;
    minimal_distance := trace_grid_minimal_distance;
  end Trace_Grid_Diagnostics;

  function In_Slice ( label,slice : integer32 ) return integer32 is

    tmp : Solution_List;

  begin
    if ((slice > ind) or else Is_Null(grid(slice))) then
      return 0;
    else
      tmp := grid(slice);
      for i in 1..Length_Of(grid(slice)) loop
        if Head_Of(tmp).m = label
         then return integer32(i);
         else tmp := Tail_Of(tmp);
        end if;
      end loop;
      return 0;
    end if;
  end In_Slice;

  function Retrieve ( label,slice : integer32 ) return Link_to_Solution is

    res : Link_to_Solution;

  begin
    if ((slice > ind ) or else Is_Null(grid(slice))) then
      return res;
    else
      res := Retrieve(grid(slice),label);
    end if;
    return res;
  end Retrieve;

  function Match ( s : Link_to_Solution; slice : integer32;
                   tol : double_float ) return integer32 is

    tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    if ((grid = null) or else (slice > grid'last)) then
      return 0;
    else
      tmp := grid(slice);
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        if Equal(ls.v,s.v,tol)
         then return ls.m;
         else tmp := Tail_Of(tmp);
        end if;
      end loop;
     return 0;
    end if;
  end Match;

  procedure Extinguish ( p : in out Vector ) is

  -- DESCRIPTION :
  --   For every i for which p(i) = 0, the corresponding j,
  --   for which p(j) = i will be set to zero as well.

    done : boolean;

  begin
    for i in p'range loop
      done := true;
      if p(i) = 0 then
        for j in p'range loop
          if integer32(p(j)) = i then
            p(j) := 0;
            done := false;
          end if;
        end loop;
      end if;
      exit when done;
    end loop;
  end Extinguish;

  function Permutation ( vrblvl : integer32 := 0 ) return Vector is

    len : constant integer32 := integer32(Length_Of(grid(0)));
    res : Vector(1..len);
    pt1,pt2 : Solution_List;
    ls1,ls2 : Link_to_Solution;

  begin
    if vrblvl > 0 then
      Write_Grid;
      put_line("Solution list 0 :");
      put(standard_output,natural32(Length_Of(grid(0))),
          natural32(Head_Of(grid(0)).n),grid(0));
      put("Solution list "); put(ind,1); put_line(" :");
      put(standard_output,natural32(Length_Of(grid(ind))),
          natural32(Head_Of(grid(ind)).n),grid(ind));
    end if;
    pt1 := grid(0);
    for i in 1..len loop
      ls1 := Head_Of(pt1);
      if ls1.m > 0 then            -- ignore clustered and failed solutions
        pt2 := grid(ind);
        res(i) := 0;
        while not Is_Null(pt2) loop
          ls2 := Head_Of(pt2);
          if ls2.m > 0 then
            if Equal(ls1.v,ls2.v,tol)
             then res(i) := natural32(ls2.m);
            end if;
          end if;
          exit when (res(i) > 0);
          pt2 := Tail_Of(pt2);
        end loop;
      end if;
      pt1 := Tail_Of(pt1);
    end loop;
    Extinguish(res);
    return res;
  end Permutation;

  function Number_of_Irreducible_Factors return natural32 is
  begin
    return Monodromy_Partitions.Number_of_Factors(deco.all);
  end Number_of_Irreducible_Factors;

  function Component ( k : integer32 ) return Link_to_Vector is

    kcnt : integer32 := 0;

  begin
    for i in deco'range loop     -- careful to skip empty entries
      if deco(i) /= null then
        kcnt := kcnt + 1;
        if kcnt = k
         then return deco(i);
        end if;
      end if;
    end loop;
    return null;
  end Component;

  function Decomposition return Link_to_VecVec is
  begin
    return deco;
  end Decomposition;

  procedure Update_Decomposition
              ( p : in Vector; nf0,nf1 : out natural32 ) is
  begin
    if deco /= null then
      nf0 := Monodromy_Partitions.Number_of_Factors(deco.all);
      nf1 := nf0;
      Monodromy_Partitions.Add_Map(deco,nf1,p);
      if not stay_silent then
        if PHCpack_Operations.Is_File_Defined then
          put(PHCpack_Operations.output_file,"the permutation : ");
          put(PHCpack_Operations.output_file,p);
          put(PHCpack_Operations.output_file," : ");
          put(PHCpack_Operations.output_file,nf0,1);
          put(PHCpack_Operations.output_file," -> ");
          put(PHCpack_Operations.output_file,nf1,1);
          new_line(PHCpack_Operations.output_file);
        else
          put(standard_output,"the permutation : ");
          put(standard_output,p);
          put(standard_output," : ");
          put(standard_output,nf0,1);
          put(standard_output," -> ");
          put(standard_output,nf1,1);
          new_line(standard_output);
        end if;
      end if;
    end if;
  end Update_Decomposition;

  function Sort ( v : Vector ) return Vector is

    r : Vector(v'range) := v;
    min_ind : integer32;
    min : natural32;

  begin
    for i in v'first..v'last-1 loop
      min_ind := i;
      for j in i+1..r'last loop  -- find the minimum
        if r(j) < r(min_ind)
         then min_ind := j;
        end if;
      end loop;
      if min_ind /= i then       -- swap if not at i-th spot
        min := r(min_ind);
        r(min_ind) := r(i);
        r(i) := min;
      end if;
    end loop;
    return r;
  end Sort;

  function Trace_Sum_Difference ( f : Vector ) return double_float is
  begin
    if stay_silent then
      return Monodromy_Component_Breakup.Trace_Sum_Difference
               (Sort(f),trace_grid);
    elsif PHCpack_Operations.Is_File_Defined then
      return Monodromy_Component_Breakup.Trace_Sum_Difference
               (PHCpack_Operations.output_file,Sort(f),trace_grid);
    else
      return Monodromy_Component_Breakup.Trace_Sum_Difference
               (standard_output,Sort(f),trace_grid);
    end if;
  end Trace_Sum_Difference;

  function Certify_with_Linear_Trace return boolean is

    done : boolean;

  begin
    if stay_silent then
      done := Monodromy_Component_Breakup.Is_Factorization
                (tol,deco.all,trace_grid);
    elsif PHCpack_Operations.Is_File_Defined then
      put_line(PHCpack_Operations.output_file,
               "Certifying with linear trace test...");
      done := Monodromy_Component_Breakup.Is_Factorization
                (PHCpack_Operations.output_file,tol,deco.all,trace_grid);
    else
      put_line(standard_output,"Certifying with linear trace test...");
      done := Monodromy_Component_Breakup.Is_Factorization
                (standard_output,tol,deco.all,trace_grid);
    end if;
    return done;
  end Certify_with_Linear_Trace;

  procedure Clear is

    procedure free is 
      new unchecked_deallocation(Array_of_Solution_Lists,
                                 Link_to_Array_of_Solution_Lists);

  begin
    if grid /= null then
      for i in 0..ind loop
        Clear(grid(i));
      end loop;
      free(grid);
    end if;
    Standard_Natural_VecVecs.Deep_Clear(deco);
  end Clear;

end Standard_Monodromy_Permutations;
