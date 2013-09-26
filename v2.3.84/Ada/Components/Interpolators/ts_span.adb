with text_io;                           use text_io;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Standard_to_Multprec_Convertors;   use Standard_to_Multprec_Convertors;
with Multprec_to_Standard_Convertors;   use Multprec_to_Standard_Convertors;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;      use Multprec_Complex_Poly_SysFun;
with Standard_Complex_Solutions;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Standard_Subspace_Restrictions;    use Standard_Subspace_Restrictions;
with Multprec_Subspace_Restrictions;    use Multprec_Subspace_Restrictions;
with Sampling_Machine;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Span_of_Component;                 use Span_of_Component;
with Span_of_Component_io;              use Span_of_Component_io;
--with Span_of_Component_Creators;        use Span_of_Component_Creators;
with Span_of_Witness_Set;               use Span_of_Witness_Set;

procedure ts_span is

-- DESCRIPTION :
--   Test on the creation of the linear subspaces spanned by a component.

  procedure Write_Span ( sp : in Standard_Span ) is

  -- DESCRIPTION :
  --   Prompts the user for a file name and writes the span sp to it.

    outfile : file_type;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    put(outfile,sp);
    Close(outfile);
  end Write_Span;

  procedure Test_Standard_io ( sp : in Standard_Span) is

  -- DESCRIPTION :
  --   Tests whether input and output formats are consistent.

    infile : file_type;
    span_read : Standard_Span;

  begin
    Write_Span(sp);
    new_line;
    put_line("Reading the name of the same file for input.");
    Read_Name_and_Open_File(infile);
    get(infile,span_read);
    Close(infile);
    put_line("The span of the component : "); put(span_read);
  end Test_Standard_io;

  procedure Test_Multprec_io ( sp : in Multprec_Span) is

  -- DESCRIPTION :
  --   Tests whether input and output formats are consistent.

    infile,outfile : file_type;
    span_read : Multprec_Span;

  begin
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    put(outfile,sp);
    Close(outfile);
    new_line;
    put_line("Reading the name of the same file for input.");
    Read_Name_and_Open_File(infile);
    get(infile,span_read);
    Close(infile);
    put_line("The span of the component : "); put(span_read);
  end Test_Multprec_io;

  procedure Test_Standard_Eval
              ( sp : in Standard_Span; tol : in double_float;
                samples : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Tests whether all samples belong to the span or not.

    tmp : Standard_Sample_List := samples;
    spt : Standard_Sample;
    equ : constant Standard_Complex_Poly_Systems.Poly_Sys := Equations(sp);
    eva : double_float;

  begin
    for i in equ'range loop
      put("equ("); put(i,1); put(") has ");
      put(Number_of_Unknowns(equ(i)),1); put_line(" unknowns.");
    end loop;
    while not Is_Null(tmp) loop
      spt := Head_Of(tmp);
      eva := Max_Norm(Eval(equ,Sample_Point(spt).v));
      put("Residual : "); put(eva);
      if In_Span(sp,tol,Sample_Point(spt).v)
       then put_line(" belongs to span.");
       else put_line(" does not belong to span.");
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Standard_Eval;

  procedure Test_Multprec_Eval
              ( sp : in Multprec_Span; tol : in double_float;
                samples : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Tests whether all samples belong to the span or not.

    tmp : Multprec_Sample_List := samples;
    spt : Multprec_Sample;
    equ : constant Multprec_Complex_Poly_Systems.Poly_Sys := Equations(sp);

  begin
    while not Is_Null(tmp) loop
      spt := Head_Of(tmp);
      declare
        eva : Multprec_Complex_Vectors.Vector(equ'range)
            := Eval(equ,Sample_Point(spt).v); 
        val : Floating_Number := Max_Norm(eva);
        fltval : constant double_float := Trunc(val);
      begin
        put("Residual : "); put(fltval);
        if In_Span(sp,tol,Sample_Point(spt).v)
         then put_line(" belongs to span.");
         else put_line(" does not belong to span.");
        end if;
        Multprec_Complex_Vectors.Clear(eva);
        Clear(val);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Multprec_Eval;

  procedure Test_Standard_Create
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_SOlutions.Solution_List;
                dim : in natural32 )  is

  -- DESCRIPTION :
  --   Test the creation of the span from samples taken from sols.

    sp : Standard_Span;
    tol : constant double_float := 1.0E-8;
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    spt : Standard_Sample;
    samples,samples_last : Standard_Sample_List;
    nb : natural32 := 0;
    ans : character;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Length_Of(sols) loop
      spt := Create(Head_Of(tmp).all,hyp);
      new_line;
      put("Sampling from generic point "); put(i,1); put_line(" :");
      loop
        new_line;
        put("Give the number of samples (0 to exit) to determine span : ");
        get(nb);
        exit when (nb = 0);
        Sample(spt,nb,samples,samples_last);
        sp := Create(samples,tol);
        if Empty(sp) then
          put_line("The span of the component is empty.");
        else
          put_line("The span of the component : "); put(sp);
          Test_Standard_Eval(sp,tol,samples);
          skip_line;
          put("Do you want to write the span to file ? ");
          Ask_Yes_or_No(ans);
          if ans = 'y'
           then Write_Span(sp);
          end if;
        end if;
        Clear(samples);
        Clear(sp);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Test_Standard_Create;

  procedure Test_Multprec_Create
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_SOlutions.Solution_List;
                dim : in natural32 )  is

    sp : Multprec_Span;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;
    tol : double_float;
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    spt : Standard_Sample;
    samples,samples_last : Multprec_Sample_List;
    nb : natural32 := 0;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
   -- tol := 10.0**(-deci+8);
    tol := 1.0E-8;
    Get_Multprec_System(p,mp,size,dim);
    Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner(size);
    for i in 1..Length_Of(sols) loop
      spt := Create(Head_Of(tmp).all,hyp);
      new_line;
      put("Sampling from generic point "); put(i,1); put_line(" :");
      loop
        new_line;
        put("Give the number of samples (0 to exit) to determine span : ");
        get(nb);
        exit when (nb = 0);
        Sample(spt,nb,samples,samples_last);
        sp := Create(samples,size,tol);
        if Empty(sp) then
          put_line("The span of the component is empty.");
        else
          put_line("The span of the component : "); put(sp);
          Test_Multprec_Eval(sp,tol,samples);
          skip_line;
          Test_Multprec_io(sp);
        end if;
        Clear(samples);
        Clear(sp);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Test_Multprec_Create;

  procedure Test_Restricted_Solutions
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Tests the residuals of the evaluation of the polynomials in
  --   the given solutions.

    tmp : Standard_Sample_List := sps;
    spt : Standard_Sample;

  begin
    for i in 1..Length_Of(tmp) loop
      spt := Head_Of(tmp);
      put("Evaluating sample point "); put(i,1); put_line(" :");
      put_line(Eval(p,Sample_Point(spt).v)); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Restricted_Solutions;

  procedure Test_Restricted_Solutions
              ( p : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sps : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Tests the residuals of the evaluation of the polynomials in
  --   the given solutions.

    tmp : Multprec_Sample_List := sps;
    spt : Multprec_Sample;

  begin
    for i in 1..Length_Of(tmp) loop
      spt := Head_Of(tmp);
      put("Evaluating sample point "); put(i,1); put_line(" :");
      put_line(Eval(p,Sample_Point(spt).v)); new_line;
      tmp := Tail_Of(tmp);
    end loop;
  end Test_Restricted_Solutions;

  procedure Test_Standard_Sample_in_Span
              ( sp : in Standard_Span; spt : in Standard_Sample;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                samples : in Standard_Sample_List;
                dim : in natural32; tol : in double_float )  is

  -- DESCRIPTION :
  --   Computes the restriction of p to the span and conducts
  --   samples with this restriced polynomial system.

    restp : constant Standard_Complex_Poly_Systems.Poly_Sys(p'range)
          := Restrict(sp,dim,tol,p);
    colrp : constant Standard_Complex_Poly_Systems.Poly_Sys
              (1..integer32(Dimension(sp)+dim))
          := Collapse_Equations(restp,Dimension(sp),dim);
    use Standard_Complex_Solutions;
    restspt : Standard_Sample;
    restsamp,newsamples,newsamples_last : Standard_Sample_List;
    nb : natural32 := 0;

  begin
   -- put_line("Polynomial System restricted to span : "); put(restp);
    restsamp := Restrict(sp,dim,samples);
   -- Test_Restricted_Solutions(restp,restsamp);
    Test_Restricted_Solutions(colrp,restsamp);
    Sampling_Machine.Initialize_Restricted(colrp);
    Sampling_Machine.Initialize_Restricted(colrp);
    restspt := Create(Sample_Point(spt),Slices(colrp,dim));
    put_line("Evaluation of start sample : ");
    put_line(Eval(colrp,Sample_Point(restspt).v));
    new_line;
    put_line("Sampling from restricted solution :");
    loop
      new_line;
      put("Give the number of samples (0 to exit) : "); get(nb);
      exit when (nb = 0);
      Sample(Standard_Output,false,restspt,nb,newsamples,newsamples_last);
    end loop;
    Sampling_Machine.Clear_Restricted;
  end Test_Standard_Sample_in_Span;
             
  procedure Test_Multprec_Sample_in_Span
              ( sp : in Multprec_Span; spt : Standard_Sample;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                samples : in Multprec_Sample_List;
                dim,deci,size : in natural32; tol : in double_float )  is

  -- DESCRIPTION :
  --   Computes the restriction of p to the span and conducts
  --   samples with this restriced polynomial system.

    timer : Timing_Widget;
    convep : constant Multprec_Complex_Poly_Systems.Poly_Sys(ep'range)
           := Convert(ep);
    mprestep : Multprec_Complex_Poly_Systems.Poly_Sys(ep'range);
    strestep : Standard_Complex_Poly_Systems.Poly_Sys(ep'range);
    restmp : Multprec_Complex_Poly_Systems.Poly_Sys(mp'range);
    colrsp : Standard_Complex_Poly_Systems.Poly_Sys
               (1..integer32(Dimension(sp)+dim));
    colrmp : Multprec_Complex_Poly_Systems.Poly_Sys
               (1..integer32(Dimension(sp)));
    restspt : Standard_Sample;
    restsamp,newsamples,newsamples_last : Multprec_Sample_List;
    nb : natural32 := 0;
   -- anykey : character;

  begin
    restsamp := Restrict(sp,dim,samples);
    restmp := Restrict(sp,0,tol,mp);
   -- put_line("Original Polynomial System restricted to Span : ");
   -- put(restmp);
   -- Test_Restricted_Solutions(restmp,restsamp);
   -- put("Give any character to continue "); get(anykey);
    colrmp := Collapse_Equations(restmp,Dimension(sp));
   -- put_line("The collapsed restricted original equations : ");
   -- put(colrmp);
   -- Test_Restricted_Solutions(colrmp,restsamp);
   -- put("Give any character to continue "); get(anykey);
    mprestep := Restrict(sp,dim,tol,convep);
   -- Test_Restricted_Solutions(mprestep,restsamp);
   -- put("Give any character to continue "); get(anykey);
    strestep := Convert(mprestep);
   -- Test_Restricted_Solutions(Convert(strestep),restsamp);
   -- put("Give any character to continue "); get(anykey);
   -- put_line("Embedded Polynomial System restricted to Span : ");
   -- put(strestep);
    colrsp := Collapse_Equations(strestep,Dimension(sp),dim);
   -- Test_Restricted_Solutions(Convert(colrsp),restsamp);
   -- put("Give any character to continue "); get(anykey);
    Sampling_Machine.Initialize_Restricted(colrsp,colrmp,integer32(dim),size);
    put_line("Sampling from restricted solution :");
    restspt := Create(Sample_Point(spt),Slices(strestep,dim));
    loop
      new_line;
      put("Give the number of samples (0 to exit) : "); get(nb);
      exit when (nb = 0);
      tstart(timer);
      Sample(Standard_Output,false,restspt,nb,newsamples,newsamples_last);
      tstop(timer);
      put("Elapsed User CPU Time : ");
      print_hms(Standard_Output,Elapsed_User_Time(timer)); new_line;
    end loop;
    Sampling_Machine.Clear_Restricted;
  end Test_Multprec_Sample_in_Span;

  procedure Evaluate_Hyperplane_Sections ( spt : in Standard_Sample ) is

    hyp : constant Standard_Complex_VecVecs.VecVec
        := Hyperplane_Sections(spt);
    sol : constant Standard_Complex_Vectors.Vector
        := Sample_Point(spt).v;
    eva : Complex_Number;

    use Standard_Complex_Vectors;
  
  begin
    for i in hyp'range loop
      put("At slice "); put(i,1); put(" : ");
      eva := hyp(i)(0) + hyp(i)(sol'range)*sol;
      put(eva); new_line;
    end loop;
  end Evaluate_Hyperplane_Sections;

  procedure Test_Slice_Restriction
                ( sp : in Standard_Span; dim : in natural32;
                  spt : in Standard_Sample ) is

  -- DESCRIPTION :
  --   Tests whether the slices of the restricted sample goes
  --   through the solution.

    res_spt : Standard_Sample := Restrict(sp,dim,spt);

  begin
    put_line("Evaluation of original slices at original sample :");
    Evaluate_Hyperplane_Sections(spt);  
    put_line("Evaluation of restricted slices at restricted sample :");
    Evaluate_Hyperplane_Sections(res_spt);  
    Deep_Clear(res_spt);
  end Test_Slice_Restriction;

  procedure Test_Slice_Restriction
                ( sp : in Multprec_Span; dim : in natural32;
                  spt : in Standard_Sample ) is

  -- DESCRIPTION :
  --   Tests whether the slices of the restricted sample goes
  --   through the solution.

    res_spt : Standard_Sample := Restrict(sp,dim,spt);

  begin
    put_line("Evaluation of original slices at original sample :");
    Evaluate_Hyperplane_Sections(spt);  
    put_line("Evaluation of restricted slices at restricted sample :");
    Evaluate_Hyperplane_Sections(res_spt);  
    Deep_Clear(res_spt);
  end Test_Slice_Restriction;

  procedure Test_Standard_Restrict
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_SOlutions.Solution_List;
                dim : in natural32 )  is

  -- DESCRIPTION :
  --   Test the creation of the span from samples taken from sols.

    sp : Standard_Span;
    tol : constant double_float := 1.0E-8;
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    spt : Standard_Sample;
    samples,samples_last : Standard_Sample_List;
    nb : natural32 := 0;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Length_Of(sols) loop
      spt := Create(Head_Of(tmp).all,hyp);
      new_line;
      put("Sampling from generic point "); put(i,1); put_line(" :");
      loop
        new_line;
        put("Give the number of samples (0 to exit) to determine span : ");
        get(nb);
        exit when (nb = 0);
        Sample(spt,nb,samples,samples_last);
        sp := Create(samples,tol);
        if Empty(sp) then
          put_line("The span of the component is empty.");
        else
          put_line("The span of the component : "); put(sp);
          Test_Slice_Restriction(sp,dim,spt);
          Test_Standard_Sample_in_Span
            (sp,Restrict(sp,dim,spt),p,samples,dim,tol);
        end if;
        Clear(samples);
        Clear(sp);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Test_Standard_Restrict;

  procedure Test_Multprec_Restrict
              ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32 )  is

    sp : Multprec_Span;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    deci,size : natural32 := 0;
    tol : double_float;
    use Standard_Complex_Solutions;
    tmp : Solution_List := sols;
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    spt : Standard_Sample;
    samples,samples_last : Multprec_Sample_List;
    nb : natural32 := 0;

  begin
    new_line;
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
   -- tol := 10.0**(-deci+8);
    tol := 1.0E-8;
    Get_Multprec_System(p,mp,size,dim);
    Sampling_Machine.Initialize(p,mp.all,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(0);
    Sampling_Machine.Default_Tune_Refiner(size);
    Sampling_Machine.Default_Tune_Refiner;
    for i in 1..Length_Of(sols) loop
      spt := Create(Head_Of(tmp).all,hyp);
      new_line;
      put("Sampling from generic point "); put(i,1); put_line(" :");
      loop
        new_line;
        put("Give the number of samples (0 to exit) to determine span : ");
        get(nb);
        exit when (nb = 0);
        Sample(spt,nb,samples,samples_last);
        sp := Create(samples,size,tol);
        if Empty(sp) then
          put_line("The span of the component is empty.");
        else
          put_line("The span of the component : "); put(sp);
          Test_Slice_Restriction(sp,dim,spt);
          Test_Multprec_Sample_in_Span
            (sp,Restrict(sp,dim,spt),p,mp.all,samples,dim,deci,size,tol);
        end if;
        Clear(samples);
        Clear(sp);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    Sampling_Machine.Clear;
  end Test_Multprec_Restrict;

  procedure Main is

    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim : natural32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing the creation of linear subspaces spanned by components");
    Standard_Read_Embedding(lp,sols,dim);
    loop
      new_line;
      put_line("MENU for testing span of components : ");
      put_line("  0. exit this program;");
      put_line("  1. interactive creation with standard arithmetic;");
      put_line("  2.                      with multiprecision arithmetic;");
      put_line("  3. test restriction with standard arithmetic;");
      put_line("  4.                  with multiprecision arithmetic;");
      put_line("  5. enumerate linear spans with standard arithmetic;");
      put_line("  6.                        with multiprecision arithmetic.");
      put("Type 0,1,2,3,4,5, or 6 to select : ");
      Ask_Alternative(ans,"0123456");
      exit when ans = '0';
      case ans is
        when '1' => Test_Standard_Create(lp.all,sols,dim);
        when '2' => Test_Multprec_Create(lp.all,sols,dim);
        when '3' => Test_Standard_Restrict(lp.all,sols,dim);
        when '4' => Test_Multprec_Restrict(lp.all,sols,dim);
        when '5' => Standard_Enumerate_Linear_Spans(lp.all,sols,dim);
        when '6' => Multprec_Enumerate_Linear_Spans(lp.all,sols,dim);
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_span;
