with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
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

package body Homotopy_Membership_Filters is

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; restol,homtol : in double_float;
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
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      if verbose
       then put("Testing point "); put(i,1); put(" : ");
      end if;
      Homotopy_Membership_Test
        (verbose,eqs,dim,sli,pts,ls.v,restol,homtol,success,found);
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
              ( file : in file_type;
                eqs : in Standard_Complex_Poly_Systems.Poly_Sys;
                pts : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; restol,homtol : in double_float;
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
      Homotopy_Membership_Test
        (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      Homotopy_Membership_Test
        (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      Homotopy_Membership_Test
        (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      Homotopy_Membership_Test
        (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      Homotopy_Membership_Test
        (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
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
                dim : in natural32; restol,homtol : in double_float;
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
      Homotopy_Membership_Test
        (file,eqs,dim,sli,pts,ls.all,restol,homtol,success,found);
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
              ( verbose : in boolean;
                eqs : in Standard_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; restol,homtol : in double_float ) is
  begin
    for dim in reverse 0..topdim-1 loop
      for witdim in reverse dim+1..topdim loop
        if not Standard_Complex_Solutions.Is_Null(pts(witdim)) then
          if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
            if verbose then
              put("Filtering junk at dimension "); put(dim,1);
              put(" with witness sets at dimension "); put(witdim,1); new_line;
            end if;   
            declare
              mempts,outpts : Standard_Complex_Solutions.Solution_List;
            begin
              Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                     restol,homtol,pts(dim),mempts,outpts);
              Standard_Complex_Solutions.Clear(mempts); -- members are junk
              Standard_Complex_Solutions.Clear(pts(witdim));
              pts(witdim) := outpts; -- points not on higher dimensional sets
            end;
          end if;
        end if;
      end loop;
    end loop;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in Standard_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out Standard_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; restol,homtol : in double_float ) is
  begin
    for dim in reverse 0..topdim-1 loop
      for witdim in reverse dim+1..topdim loop
        if not Standard_Complex_Solutions.Is_Null(pts(witdim)) then
          if not Standard_Complex_Solutions.Is_Null(pts(dim)) then
            if verbose then
              put("Filtering junk at dimension "); put(dim,1);
              put(" with witness sets at dimension "); put(witdim,1); new_line;
            end if;   
            declare
              mempts,outpts : Standard_Complex_Solutions.Solution_List;
            begin
              Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                     restol,homtol,pts(dim),mempts,outpts);
              Standard_Complex_Solutions.Clear(mempts); -- members are junk
              Standard_Complex_Solutions.Clear(pts(witdim));
              pts(witdim) := outpts; -- points not on higher dimensional sets
            end;
          end if;
        end if;
      end loop;
    end loop;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; restol,homtol : in double_float ) is
  begin
    for dim in reverse 0..topdim-1 loop
      for witdim in reverse dim+1..topdim loop
        if not DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
          if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
            if verbose then
              put("Filtering junk at dimension "); put(dim,1);
              put(" with witness sets at dimension "); put(witdim,1); new_line;
            end if;   
            declare
              mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
            begin
              Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                     restol,homtol,pts(dim),mempts,outpts);
              DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
              DoblDobl_Complex_Solutions.Clear(pts(witdim));
              pts(witdim) := outpts; -- points not on higher dimensional sets
            end;
          end if;
        end if;
      end loop;
    end loop;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in DoblDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out DoblDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; restol,homtol : in double_float ) is
  begin
    for dim in reverse 0..topdim-1 loop
      for witdim in reverse dim+1..topdim loop
        if not DoblDobl_Complex_Solutions.Is_Null(pts(witdim)) then
          if not DoblDobl_Complex_Solutions.Is_Null(pts(dim)) then
            if verbose then
              put("Filtering junk at dimension "); put(dim,1);
              put(" with witness sets at dimension "); put(witdim,1); new_line;
            end if;   
            declare
              mempts,outpts : DoblDobl_Complex_Solutions.Solution_List;
            begin
              Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                     restol,homtol,pts(dim),mempts,outpts);
              DoblDobl_Complex_Solutions.Clear(mempts); -- members are junk
              DoblDobl_Complex_Solutions.Clear(pts(witdim));
              pts(witdim) := outpts; -- points not on higher dimensional sets
            end;
          end if;
        end if;
      end loop;
    end loop;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Poly_Systems.Array_of_Poly_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; restol,homtol : in double_float ) is
  begin
    for dim in reverse 0..topdim-1 loop
      for witdim in reverse dim+1..topdim loop
        if not QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
          if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
            if verbose then
              put("Filtering junk at dimension "); put(dim,1);
              put(" with witness sets at dimension "); put(witdim,1); new_line;
            end if;   
            declare
              mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
            begin
              Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                     restol,homtol,pts(dim),mempts,outpts);
              QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
              QuadDobl_Complex_Solutions.Clear(pts(witdim));
              pts(witdim) := outpts; -- points not on higher dimensional sets
            end;
          end if;
        end if;
      end loop;
    end loop;
  end Filter;

  procedure Filter
              ( verbose : in boolean;
                eqs : in QuadDobl_Complex_Laur_Systems.Array_of_Laur_Sys;
                pts : in out QuadDobl_Complex_Solutions.Array_of_Solution_Lists;
                topdim : in integer32; restol,homtol : in double_float ) is
  begin
    for dim in reverse 0..topdim-1 loop
      for witdim in reverse dim+1..topdim loop
        if not QuadDobl_Complex_Solutions.Is_Null(pts(witdim)) then
          if not QuadDobl_Complex_Solutions.Is_Null(pts(dim)) then
            if verbose then
              put("Filtering junk at dimension "); put(dim,1);
              put(" with witness sets at dimension "); put(witdim,1); new_line;
            end if;   
            declare
              mempts,outpts : QuadDobl_Complex_Solutions.Solution_List;
            begin
              Filter(verbose,eqs(witdim).all,pts(witdim),natural32(witdim),
                     restol,homtol,pts(dim),mempts,outpts);
              QuadDobl_Complex_Solutions.Clear(mempts); -- members are junk
              QuadDobl_Complex_Solutions.Clear(pts(witdim));
              pts(witdim) := outpts; -- points not on higher dimensional sets
            end;
          end if;
        end if;
      end loop;
    end loop;
  end Filter;

end Homotopy_Membership_Filters;
