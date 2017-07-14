with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with QuadDobl_Random_Numbers;           use QuadDobl_Random_Numbers;
with QuadDobl_Complex_Vectors;          use QuadDobl_Complex_Vectors;
with Witness_Sets;                      use Witness_Sets;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with QuadDobl_Solutions_Container;

package body QuadDobl_Sampling_Operations is

-- INTERNAL DATA :

  start_samples,new_samples : Solution_List;
  original_solutions : Solution_List;
  dimension : natural32;
  start_slices,new_slices : Link_to_VecVec;
  gamma : Link_to_Vector;
  hyperplane_sections : Link_to_Array_of_VecVecs;
  nb_sections : natural32 := 0;
  use_laurent : boolean := false;

-- OPERATIONS :

  procedure Initialize ( p : in Poly_Sys; sols : in Solution_List;
                         k : in integer32 ) is

    n : constant integer32 := p'last;

  begin
    use_laurent := false;
    QuadDobl_Sampling_Machine.Initialize(p);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
   -- put("initializing sampling machine with #solutions = ");
   -- put(Length_Of(sols),1); new_line;
    dimension := natural32(k);
    start_slices := new VecVec'(Slices(p,natural32(k)));
    new_slices := new VecVec(1..k);
    for i in 1..k loop
      new_slices(i) := new Vector'(start_slices(i).all);
    end loop;
    gamma := new Vector(1..n);
    for i in 1..n loop
      gamma(i) := Create(integer32(1));
    end loop;
   -- put("  k = "); put(k,1);
   -- put("  d = "); put(Length_Of(sols),1); new_line;
    Copy(sols,original_solutions);
    Copy(original_solutions,start_samples);
           -- start_samples := original_solutions;
   -- put("number of start samples : ");
   -- put(Length_Of(start_samples),1); new_line;
  end Initialize;

  procedure Initialize ( p : in Laur_Sys; sols : in Solution_List;
                         k : in integer32 ) is

    n : constant integer32 := p'last;

  begin
    use_laurent := true;
    QuadDobl_Sampling_Laurent_Machine.Initialize(p);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
   -- put("initializing sampling machine with #solutions = ");
   -- put(Length_Of(sols),1); new_line;
    dimension := natural32(k);
    start_slices := new VecVec'(Slices(p,natural32(k)));
    new_slices := new VecVec(1..k);
    for i in 1..k loop
      new_slices(i) := new Vector'(start_slices(i).all);
    end loop;
    gamma := new Vector(1..n);
    for i in 1..n loop
      gamma(i) := Create(integer32(1));
    end loop;
   -- put("  k = "); put(k,1);
   -- put("  d = "); put(Length_Of(sols),1); new_line;
    Copy(sols,original_solutions);
    Copy(original_solutions,start_samples);
           -- start_samples := original_solutions;
   -- put("number of start samples : ");
   -- put(Length_Of(start_samples),1); new_line;
  end Initialize;

  procedure Initialize_Slices ( n : in integer32 ) is
  begin
    hyperplane_sections := new Array_of_VecVecs(0..n);
    if start_slices /= null then
      declare
        cp : VecVec(start_slices'range); -- copy of the start_slices
      begin
        for i in cp'range loop
          cp(i) := new Vector'(start_slices(i).all);
        end loop;
        hyperplane_sections(0) := new VecVec'(cp);
      end;
      -- put_line("initialized sections #0 with :");
      -- put_line(hyperplane_sections(nb_sections));
    end if;
    nb_sections := 0;
  end Initialize_Slices;

  function Retrieve_Slices ( i : integer32 ) return Link_to_VecVec is
  begin
    if hyperplane_sections = null then
      return null;
    elsif i > hyperplane_sections'last then 
      return null;
    else
      return hyperplane_sections(i);
    end if;
  end Retrieve_Slices;

  procedure Add_Slices ( v : in VecVec ) is
  begin
    if hyperplane_sections /= null then
      if integer32(nb_sections) < hyperplane_sections'last then
        nb_sections := nb_sections + 1;
        hyperplane_sections(integer32(nb_sections)) := new VecVec'(v);
       -- put("Adding sections #"); put(nb_sections,1); put_line(":");
       -- put_line(hyperplane_sections(nb_sections));
      end if;
    end if;
  end Add_Slices;

  procedure Set_Target_Slices ( i : integer32 ) is
  begin
    if hyperplane_sections /= null then
      if i <= hyperplane_sections'last then
       -- put("Setting the target slices to sections "); put(i,1); 
       -- put_line(" ...");
       -- put_line(hyperplane_sections(i));
        for j in new_slices'range loop
          for k in new_slices(j)'range loop
            start_slices(j)(k) := hyperplane_sections(0)(j)(k);
            new_slices(j)(k) := hyperplane_sections(i)(j)(k);
          end loop;
        end loop;
      end if;
    end if;
    if use_laurent
     then QuadDobl_Sampling_Laurent_Machine.Change_Slices(start_slices.all);
     else QuadDobl_Sampling_Machine.Change_Slices(start_slices.all);
    end if;
   -- put_line("changed slices in QuadDobl_Sampling_Machine, leaving ...");
  end Set_Target_Slices;

  function Retrieve_First_Solutions return Solution_List is
  begin
    return original_solutions;
  end Retrieve_First_Solutions;

  function Retrieve_Start_Slices return Link_to_VecVec is
  begin
    return start_slices;
  end Retrieve_Start_Slices;

  function Retrieve_Target_Slices return Link_to_VecVec is
  begin
    return new_slices;
  end Retrieve_Target_Slices;

  procedure Assign_Slice ( c : in Complex_Number; i,j : in integer32 ) is
  begin
    new_slices(i+1)(j) := c;
   -- put("received "); put(c);
   -- put(" to assign at ("); put(i+1,1); put(","); put(j,1);
   -- put_line(")");
  end Assign_Slice;

  procedure Store_Gamma ( c : in Complex_Number; i : in integer32 ) is
  begin
    gamma(i+1) := c;
   -- put("received "); put(c);
   -- put(" as gamma constant "); put(i+1,1); new_line;
  end Store_Gamma;

  procedure Swap_Slices is

    tmp_slices : constant Link_to_VecVec := new_slices;
    tmp_samples : constant Solution_List := new_samples;

  begin
    new_slices := start_slices;
    start_slices := tmp_slices;
    new_samples := start_samples;
    start_samples := tmp_samples;
    Set_Continuation_Parameter(start_samples,Create(integer32(0)));
    if use_laurent
     then QuadDobl_Sampling_Laurent_Machine.Change_Slices(start_slices.all);
     else QuadDobl_Sampling_Machine.Change_Slices(start_slices.all);
    end if;
  end Swap_Slices;

  procedure Sample is
  begin
   -- put_line("Sampling with gamma =");
   -- put_line(gamma.all);
   -- put_line("Sampling with slices :"); put(new_slices.all);
   -- put("Number of start solutions : ");
   -- put(Length_Of(start_samples),1); new_line;
    Set_Continuation_Parameter(start_samples,Create(integer32(0)));
    if use_laurent then
      QuadDobl_Sampling_Laurent_Machine.Sample
        (start_samples,new_slices.all,gamma.all,new_samples);
    else
      QuadDobl_Sampling_Machine.Sample
        (start_samples,new_slices.all,gamma.all,new_samples);
    end if;
   -- put_line("Computed solutions : ");
   -- put(standard_output,
   --     Length_Of(new_samples),Head_Of(new_samples).n,new_samples);
    QuadDobl_Solutions_Container.Clear;
    QuadDobl_Solutions_Container.Initialize(new_samples);
  end Sample;

  function Sample_Loop 
              ( start,target : integer32; s : Link_to_Solution )
              return Link_to_Solution is

    start_sols,target_sols : Solution_List;

  begin
    if ((hyperplane_sections = null)
        or else (start > hyperplane_sections'last)
        or else (target > hyperplane_sections'last)) then
      return s;
    else
      for i in start_slices'range loop
        for j in start_slices(i)'range loop
          start_slices(i)(j) := hyperplane_sections(start)(i)(j);
          new_slices(i)(j) := hyperplane_sections(target)(i)(j);
        end loop;
      end loop;
      -- put_line("The start solution :"); put(s.all);
      Construct(s,start_sols);
      Set_Continuation_Parameter(start_sols,Create(integer32(0)));
      if use_laurent
       then QuadDobl_Sampling_Laurent_Machine.Change_Slices(start_slices.all);
       else QuadDobl_Sampling_Machine.Change_Slices(start_slices.all);
      end if;
      for i in gamma'range loop
        gamma(i) := Random1;
      end loop;
      if use_laurent then
        QuadDobl_Sampling_Laurent_Machine.Sample
          (start_sols,new_slices.all,gamma.all,target_sols);
      else
        QuadDobl_Sampling_Machine.Sample
          (start_sols,new_slices.all,gamma.all,target_sols);
      end if;
      -- put_line("The computed solution :");
      -- put(Head_Of(target_sols).all);
      return Head_Of(target_sols);
    end if;
  end Sample_Loop;

  procedure Clear is
  begin
    if use_laurent
     then QuadDobl_Sampling_Laurent_Machine.Clear;
     else QuadDobl_Sampling_Machine.Clear;
    end if;
    dimension := 0;
  end Clear;

end QuadDobl_Sampling_Operations;
