with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Sampling_Machine;
with Sampling_Laurent_Machine;
with DoblDobl_Sampling_Machine;
with DoblDobl_Sampling_Laurent_Machine;
with QuadDobl_Sampling_Machine;
with QuadDobl_Sampling_Laurent_Machine;
with Witness_Sets;
with Homotopy_Membership_Tests;          use Homotopy_Membership_Tests;
with Multitasking_Membership_Tests;

package body Homotopy_Membership_Filters is

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := totest;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Machine.Initialize(eqs);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      if verbose then
        put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          Homotopy_Membership_Test
            (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
        else
          success := false; -- diverged path is not retained
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(Standard_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    Sampling_Machine.Clear;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := totest;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Laurent_Machine.Initialize(eqs);
    Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          Homotopy_Membership_Test
            (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
        else
          success := false;     -- do not retain diverging path
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(Standard_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    Sampling_Laurent_Machine.Clear;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Machine.Initialize(eqs);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(2);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if (Double_Double_Numbers.hi_part(ls.rco) < rcotol) then
        if (Double_Double_Numbers.hi_part(ls.res) < restol) then
          Homotopy_Membership_Test
            (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(DoblDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    DoblDobl_Sampling_Machine.Clear;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(eqs);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if (Double_Double_Numbers.hi_part(ls.rco) < rcotol) then
        if (Double_Double_Numbers.hi_part(ls.res) < restol) then
          Homotopy_Membership_Test
            (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(DoblDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    DoblDobl_Sampling_Laurent_Machine.Clear;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Machine.Initialize(eqs);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if (Quad_Double_Numbers.hihi_part(ls.rco) < rcotol) then
        if (Quad_Double_Numbers.hihi_part(ls.res) < restol) then
          Homotopy_Membership_Test
            (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(QuadDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    QuadDobl_Sampling_Machine.Clear;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(eqs);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      if verbose
       then put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if (Quad_Double_Numbers.hihi_part(ls.rco) < rcotol) then
        if (Quad_Double_Numbers.hihi_part(ls.res) < restol) then
          Homotopy_Membership_Test
            (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(QuadDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
    QuadDobl_Sampling_Laurent_Machine.Clear;
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : natural32;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := totest;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      if verbose then
        put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          idx := Standard_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(Standard_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : natural32;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := totest;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      if verbose then
        put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          idx := Standard_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(Standard_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      if verbose then
        put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if Double_Double_Numbers.hi_part(ls.rco) < rcotol then
        if Double_Double_Numbers.hi_part(ls.res) < restol then
          idx := DoblDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(DoblDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      if verbose then
        put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if Double_Double_Numbers.hi_part(ls.rco) < rcotol then
        if Double_Double_Numbers.hi_part(ls.res) < restol then
          idx := DoblDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := true; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(DoblDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      if verbose then
        put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if Quad_Double_Numbers.hihi_part(ls.rco) < rcotol then
        if Quad_Double_Numbers.hihi_part(ls.res) < restol then
          idx := QuadDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(QuadDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := totest;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      if verbose then
        put("Testing point "); put(i,1); put_line(" : ");
      end if;
      if Quad_Double_Numbers.hihi_part(ls.rco) < rcotol then
        if Quad_Double_Numbers.hihi_part(ls.res) < restol then
          idx := QuadDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        if verbose then
          put("The point is considered regular for its tolerance ");
          put(ls.rco,3); put(" < "); put(rcotol,3); put_line(".");
        end if;
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    if verbose then
      put("Tested ");
      put(QuadDobl_Complex_Solutions.Length_Of(totest),1); 
      put_line(" candidates for solutions:");
      put("  Found "); put(cnt,1);
      put_line(" solutions on the components.");
    end if;
  end Filter;

  procedure Filter
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := pts;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Machine.Initialize(eqs);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          Homotopy_Membership_Test
            (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,Standard_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    Sampling_Machine.Clear;
  end Filter;

  procedure Filter
              ( file : in file_type;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := pts;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    sli : constant Standard_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    Sampling_Laurent_Machine.Initialize(eqs);
    Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          Homotopy_Membership_Test
            (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
        else
          success := false;
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,Standard_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    Sampling_Laurent_Machine.Clear;
  end Filter;

  procedure Filter
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Machine.Initialize(eqs);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(2);
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Double_Double_Numbers.hi_part(ls.rco) < rcotol then
        if Double_Double_Numbers.hi_part(ls.res) < restol then
          Homotopy_Membership_Test
            (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,DoblDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    DoblDobl_Sampling_Machine.Clear;
  end Filter;

  procedure Filter
              ( file : in file_type;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    sli : constant DoblDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    DoblDobl_Sampling_Laurent_Machine.Initialize(eqs);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    DoblDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Double_Double_Numbers.hi_part(ls.rco) < rcotol then
        if Double_Double_Numbers.hi_part(ls.res) < restol then
          Homotopy_Membership_Test
            (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,DoblDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    DoblDobl_Sampling_Laurent_Machine.Clear;
  end Filter;

  procedure Filter
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Machine.Initialize(eqs);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Quad_Double_Numbers.hihi_part(ls.rco) < rcotol then
        if Quad_Double_Numbers.hihi_part(ls.res) < restol then
          Homotopy_Membership_Test
            (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,QuadDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    QuadDobl_Sampling_Machine.Clear;
  end Filter;

  procedure Filter
              ( file : in file_type;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    sli : constant QuadDobl_Complex_VecVecs.VecVec
        := Witness_Sets.Slices(eqs,dim);
    success,found : boolean;
    cnt : natural32 := 0;

  begin
    QuadDobl_Sampling_Laurent_Machine.Initialize(eqs);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Sampler(2);
    QuadDobl_Sampling_Laurent_Machine.Default_Tune_Refiner;
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Quad_Double_Numbers.hihi_part(ls.rco) < rcotol then
        if Quad_Double_Numbers.hihi_part(ls.res) < restol then
          Homotopy_Membership_Test
            (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,QuadDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
    QuadDobl_Sampling_Laurent_Machine.Clear;
  end Filter;

  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := pts;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          idx := Standard_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,Standard_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
  end Filter;

  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Laur_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in Standard_Complex_Solutions.Solution_List;
                mempts : out Standard_Complex_Solutions.Solution_List;
                outpts : out Standard_Complex_Solutions.Solution_List ) is

    tmp : Standard_Complex_Solutions.Solution_List := pts;
    mempts_last : Standard_Complex_Solutions.Solution_List := mempts;
    outpts_last : Standard_Complex_Solutions.Solution_List := outpts;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..Standard_Complex_Solutions.Length_Of(totest) loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if ls.rco < rcotol then   -- do not test regular points
        if ls.res < restol then -- do not test diverged paths
          idx := Standard_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (ls.res < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          Standard_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          Standard_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,Standard_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
  end Filter;

  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Double_Double_Numbers.hi_part(ls.rco) < rcotol then
        if Double_Double_Numbers.hi_part(ls.res) < restol then
          idx := DoblDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,DoblDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
  end Filter;

  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in DoblDobl_Complex_Solutions.Solution_List;
                mempts : out DoblDobl_Complex_Solutions.Solution_List;
                outpts : out DoblDobl_Complex_Solutions.Solution_List ) is

    tmp : DoblDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : DoblDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : DoblDobl_Complex_Solutions.Solution_List := outpts;
    ls : DoblDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..DoblDobl_Complex_Solutions.Length_Of(totest) loop
      ls := DoblDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Double_Double_Numbers.hi_part(ls.rco) < rcotol then
        if Double_Double_Numbers.hi_part(ls.res) < restol then
          idx := DoblDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Double_Double_Numbers.hi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          DoblDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          DoblDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := DoblDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,DoblDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
  end Filter;

  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Quad_Double_Numbers.hihi_part(ls.rco) < rcotol then
        if Quad_Double_Numbers.hihi_part(ls.res) < restol then
          idx := QuadDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,QuadDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
  end Filter;

  procedure Filter
              ( file : in file_type; nt : in natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                pts : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; rcotol,restol,homtol : in double_float;
                totest : in QuadDobl_Complex_Solutions.Solution_List;
                mempts : out QuadDobl_Complex_Solutions.Solution_List;
                outpts : out QuadDobl_Complex_Solutions.Solution_List ) is

    tmp : QuadDobl_Complex_Solutions.Solution_List := pts;
    mempts_last : QuadDobl_Complex_Solutions.Solution_List := mempts;
    outpts_last : QuadDobl_Complex_Solutions.Solution_List := outpts;
    ls : QuadDobl_Complex_Solutions.Link_to_Solution;
    success,found : boolean;
    cnt,idx : natural32 := 0;

    use Multitasking_Membership_Tests;

  begin
    for i in 1..QuadDobl_Complex_Solutions.Length_Of(totest) loop
      ls := QuadDobl_Complex_Solutions.Head_Of(tmp);
      put(file,"Testing point "); put(file,i,1); put(file," : ");
      if Quad_Double_Numbers.hihi_part(ls.rco) < rcotol then
        if Quad_Double_Numbers.hihi_part(ls.res) < restol then
          idx := QuadDobl_Membership_Test(nt,dim,homtol,eqs,pts,ls.v);
          success := true;
          found := (idx /= 0);
        else
          success := false; -- diverged path
        end if;
      else
        success := (Quad_Double_Numbers.hihi_part(ls.res) < restol);
        found := false;  -- the point does not lie on the solution set
        put(file,"The point is considered regular for its tolerance ");
        put(file,ls.rco,3); put(file," < ");
        put(file,rcotol,3); put_line(file,".");
      end if;
      if success then
        if found then
          cnt := cnt + 1;
          QuadDobl_Complex_Solutions.Append(mempts,mempts_last,ls.all);
        else
          QuadDobl_Complex_Solutions.Append(outpts,outpts_last,ls.all);
        end if;
      end if;
      tmp := QuadDobl_Complex_Solutions.Tail_Of(tmp);
    end loop;
    put(file,"Tested ");
    put(file,QuadDobl_Complex_Solutions.Length_Of(totest),1); 
    put_line(file," candidates for solutions:");
    put(file,"  Found "); put(file,cnt,1);
    put_line(file," solutions on the components.");
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop -- removing junk at dimension dim
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                       rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                       rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                       rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                       rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                       rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                       rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop -- removing junk at dimension dim
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := Standard_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := Standard_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if Standard_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : Standard_Complex_Solutions.Solution_List;
              begin
                Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                Standard_Complex_Solutions.Clear(mempts); -- members are junk
                Standard_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := Standard_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
              begin
                Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
                DoblDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := DoblDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

  procedure Filter
              ( verbose : in boolean; nt : in natural32;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; rcotol,restol,homtol : in double_float;
                filcnt : out Standard_Natural_VecVecs.VecVec;
                times : out Array_of_Duration; alltime : out duration ) is

    total_timer,level_timer : Timing_Widget;
    topcnt : Standard_Natural_Vectors.Vector(0..0);

  begin
    tstart(total_timer);
    topcnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(topdim));
    filcnt(topdim) := new Standard_Natural_Vectors.Vector'(topcnt);
    for dim in reverse 0..topdim-1 loop
      tstart(level_timer);
      declare
        cnt : Standard_Natural_Vectors.Vector(0..topdim-dim);
        idx : integer32 := 0;
      begin
        cnt := (cnt'range => 0);
        cnt(0) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
        for witdim in reverse dim+1..topdim loop
          idx := idx + 1;
          if QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
            cnt(idx) := cnt(idx-1); -- no witness set at witdim, no junk
          else
            if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
              if verbose then
                put("Filtering junk at dimension "); put(dim,1);
                put(" with witness sets at dimension "); put(witdim,1);
                new_line;
              end if;   
              declare
                mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
              begin
                Filter
                  (verbose,nt,eqs(witdim).all,pts(witdim),natural32(witdim),
                   rcotol,restol,homtol,pts(dim),mempts,outpts);
                QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
                QuadDobl_Complex_Solutions.Clear(pts(dim));
                pts(dim) := outpts; -- points not on higher dimensional sets
                cnt(idx) := QuadDobl_Complex_Solutions.Length_Of(pts(dim));
                if verbose then
                  put("Junk removal : "); put(cnt(idx-1),1);
                  put(" -> "); put(cnt(idx),1); new_line;
                end if;
              end;
            end if;
          end if;
        end loop;
        filcnt(dim) := new Standard_Natural_Vectors.Vector'(cnt);
      end;
      tstop(level_timer);
      times(integer(dim)) := Elapsed_User_Time(level_timer);
    end loop;
    tstop(total_timer);
    alltime := Elapsed_User_Time(total_timer);
  end Filter;

end Homotopy_Membership_Filters;
