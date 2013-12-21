with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Multprec_Floating_Matrices_io;     use Multprec_Floating_Matrices_io;
with Standard_Complex_NesVecs_io;       use Standard_Complex_NesVecs_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with Sampling_Machine;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Drivers_to_Grid_Creators;          use Drivers_to_Grid_Creators;
with Standard_Stacked_Sample_Grids;
with Multprec_Stacked_Sample_Grids;
with Standard_Nvariate_Interpolators;
with Standard_Trace_Interpolators;
with Multprec_Trace_Interpolators;

procedure ts_traces is

-- DESCRIPTION :
--   Interactive test facility for the multivariate interpolation with traces.

  procedure Standard_Write_Evaluations
              ( file : in file_type;
                t : in Standard_Trace_Interpolators.Trace_Interpolator1;
                sps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.
  --   This version is for the special case of dealing with a curve.

    use Standard_Trace_Interpolators;
    tmp : Standard_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Write_Evaluations;

  procedure Standard_Write_Evaluations
              ( file : in file_type;
                t : in Standard_Trace_Interpolators.Trace_Interpolator;
                sps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.
  --   This version accepts the trace form for general surfaces.

    use Standard_Trace_Interpolators;
    tmp : Standard_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Standard_Write_Evaluations;

  procedure Multprec_Write_Evaluations
              ( file : in file_type;
                t : in Multprec_Trace_Interpolators.Trace_Interpolator1;
                sps : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.

    use Multprec_Trace_Interpolators;
    tmp : Multprec_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Multprec_Write_Evaluations;

  procedure Multprec_Write_Evaluations
              ( file : in file_type;
                t : in Multprec_Trace_Interpolators.Trace_Interpolator;
                sps : in Multprec_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.
  --   This version accepts the trace form for general surfaces.

    use Multprec_Trace_Interpolators;
    tmp : Multprec_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end Multprec_Write_Evaluations;

  procedure Standard_Trace_Interpolation1
               ( file : in file_type;
                 grid : in Array_of_Standard_Sample_Lists;
                 testsps : in Standard_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.
  --   This routine is restricted to the special case of curves.

    use Standard_Complex_Numbers,Standard_Complex_Vectors;
    use Standard_Complex_Polynomials,Standard_Trace_Interpolators;

    t : Trace_Interpolator1;
    p : Poly;

  begin
    if output
     then t := Create(file,grid);
     else t := Create(grid);
    end if;
    Write_Errors(file,t,grid);
    put_line(file,"The matrix with residuals at the grid : ");
    put(file,Errors(t,grid),3);
    put(file,"Maximal residual of evaluation at the grid : ");
    put(file,Maximal_Error(t,grid),3); new_line(file);
    put_line(file,"Evaluation at the test points : ");
    Standard_Write_Evaluations(file,t,testsps);
    put_line(file,"Test on componentwise construction of traces :");
    for i in 1..Degree(t) loop
      put(file,"Trace of degree "); put(file,i,1); put_line(file," :");
      put_line(file,Trace(t,i));
      declare
        tr : constant Vector := Create(grid,i);
        val,eva : Complex_Number;
      begin
        put_line(file,"The constructed trace : "); put_line(file,tr);
        Eval_Trace(tr,i,testsps,val,eva);
        put(file,"value of trace   : "); put(file,val); new_line(file);
        put(file,"power sum result : "); put(file,eva); new_line(file);
      end;
    end loop;
    p := Expand(t);
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,p);
  end Standard_Trace_Interpolation1;

  procedure Standard_Trace_Interpolation
              ( file : in file_type; d : in integer32;
                grid : in Standard_Stacked_Sample_Grids.Stacked_Sample_Grid;
                testsps : in Standard_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.

    use Standard_Complex_Numbers,Standard_Complex_Polynomials;
    use Standard_Nvariate_Interpolators,Standard_Trace_Interpolators;

    t : Trace_Interpolator;
    maxerr : double_float;
    p : Poly;

  begin
    if output
     then t := Create(file,grid,d);
     else t := Create(grid,d);
    end if;
    Write_Errors(file,t,grid,maxerr);
    put(file,"Maximal residual of evaluation at the grid : ");
  -- put(file,Maximal_Error(t,grid),3); new_line(file);
    put(file,maxerr,3); new_line(file);
    put_line(file,"Evaluation at the test points : ");
    Standard_Write_Evaluations(file,t,testsps);
    put_line(file,"Test on construction of linear trace :");
    declare
      seltr : constant Newton_Form(natural32(grid.k),grid.k,1) := Trace(t,1);
      tr : Newton_Form(natural32(grid.k),grid.k,1);
      val,eva : Complex_Number;
      diff : double_float;
    begin
      if output
       then tr := Create(file,grid,d,1);
       else tr := Create(grid,d,1);
      end if;
      put_line(file,"The divided differences : "); put(file,tr.f);
      put_line(file,"The evaluation at the newly computed trace : ");
      Eval_Trace(tr,1,testsps,val,eva);
      put(file,"value of trace   : "); put(file,val); new_line(file);
      put(file,"power sum result : "); put(file,eva); new_line(file);
      diff := AbsVal(val-eva);
      put(file,"-> difference between value and power sum : ");
      put(file,diff);
      new_line(file);
      put_line(file,"The evaluation at the selected linear trace : ");
      Eval_Trace(seltr,1,testsps,val,eva);
      put(file,"value of trace   : "); put(file,val); new_line(file);
      put(file,"power sum result : "); put(file,eva); new_line(file);
      diff := AbsVal(val-eva);
      put(file,"-> difference between value and power sum : ");
      put(file,diff);
      new_line(file);
    end;
    p := Expand(t);
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,p);
  end Standard_Trace_Interpolation;

  procedure Multprec_Trace_Interpolation1
              ( file : in file_type;
                grid : in Array_of_Multprec_Sample_Lists;
                testsps : in Multprec_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.
  --   This routine is restricted to the special case of curves.

    use Multprec_Complex_Numbers,Multprec_Complex_Vectors;
    use Multprec_Complex_Polynomials,Multprec_Trace_Interpolators;
    t : Trace_Interpolator1;
    p : Poly;

  begin
    if output
     then t := Create(file,grid);
     else t := Create(grid);
    end if;
    Write_Errors(file,t,grid);
    put_line(file,"The matrix with residuals at the grid : ");
    put(file,Errors(t,grid),3);
    put(file,"Maximal residual of evaluation at the grid : ");
    put(file,Maximal_Error(t,grid),3); new_line(file);
    put_line(file,"Evaluation at the test points : ");
    Multprec_Write_Evaluations(file,t,testsps);
    put_line(file,"Test on componentwise construction of traces :");
    for i in 1..Degree(t) loop
      put(file,"Trace of degree "); put(file,i,1); put_line(file," :");
      put_line(file,Trace(t,i));
      declare
        tr : Vector(0..i) := Create(grid,i);
        val,eva : Complex_Number;
      begin
        put_line(file,"The constructed trace : "); put_line(file,tr);
        Eval_Trace(tr,i,testsps,val,eva);
        put(file,"value of trace   : "); put(file,val); new_line(file);
        put(file,"power sum result : "); put(file,eva); new_line(file);
        Clear(tr); Clear(val); Clear(eva);
      end;
    end loop;
    p := Expand(t);
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,p);
  end Multprec_Trace_Interpolation1;

  procedure Multprec_Trace_Interpolation
              ( file : in file_type; d : in integer32;
                grid : in Multprec_Stacked_Sample_Grids.Stacked_Sample_Grid;
                testsps : in Multprec_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.

    use Multprec_Complex_Polynomials;
    use Multprec_Trace_Interpolators;

    t : Trace_Interpolator;
    maxerr : Floating_Number;
    p : Poly;

  begin
    if output
     then t := Create(file,grid,d);
     else t := Create(grid,d);
    end if;
    Write_Errors(file,t,grid,maxerr);
    put(file,"Maximal residual of evaluation at the grid : ");
    put(file,maxerr,3); new_line(file);
    put_line(file,"Evaluation at the test points : ");
    Multprec_Write_Evaluations(file,t,testsps);
    p := Expand(t);
    put_line(file,"The trace interpolator expanded as polynomial : ");
    put_line(file,p);                                                    
  end Multprec_Trace_Interpolation;

  procedure Test_Standard_Trace_Interpolation1
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test on the trace interpolation for curves with standard arithmetic.

    sli : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : Standard_Sample_List := Create(sols,sli);
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Standard_Sample_Lists(0..integer32(len));
    eps,dst : double_float;
    testsps : Standard_Sample_List;

  begin
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    Standard_Test_Samples(file,sps,sli,testsps);
    Standard_Rectangular_Grid_Creator(file,sps,len,grid,eps,dst);
    Standard_Trace_Interpolation1(file,grid,testsps,output);
    Sampling_Machine.Clear;
  end Test_Standard_Trace_Interpolation1;

  procedure Test_Standard_Trace_Interpolation
              ( file : in file_type;
                p : in Standard_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test on the trace interpolation for general surfaces,
  --   with standard arithmetic.

    use Standard_Stacked_Sample_Grids;

    deg : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    gsz : constant natural32 := Full_Grid_Size(dim+1,deg);
    sli : Standard_Complex_VecVecs.VecVec(1..integer32(dim)) := Slices(p,dim);
    sps : Standard_Sample_List := Create(sols,sli);
    testsps : Standard_Sample_List;
    grid : Stacked_Sample_Grid(integer32(dim),integer32(deg));

  begin
    put("Number of samples : "); put(gsz,1); new_line;
    new_line;
    Sampling_Machine.Initialize(p);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    new_line;
    new_line;
    put_line("See the output file for results...");
    new_line;
    Standard_Test_Samples(file,sps,sli,testsps);
    Standard_Stacked_Grid_Creator(file,sps,true,grid);
    Standard_Trace_Interpolation(file,integer32(deg),grid,testsps,output);
    Sampling_Machine.Clear;
  end Test_Standard_Trace_Interpolation;

  procedure Test_Multprec_Trace_Interpolation1
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim,size : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test trace interpolation for curves with multi-precision arithmetic.

    sli : Standard_Complex_VecVecs.VecVec(1..integer32(dim)) := Slices(ep,dim);
    sps : Standard_Sample_List := Create(sols,sli);
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_Multprec_Sample_Lists(0..integer32(len));
    eps,dst : double_float;
    testsps : Multprec_Sample_List;

  begin
    Sampling_Machine.Initialize(ep,mp,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Multprec_Test_Samples(file,sps,sli,testsps);
    Multprec_Rectangular_Grid_Creator(file,sps,len,size,grid,eps,dst);
    Multprec_Trace_Interpolation1(file,grid,testsps,output);
    Sampling_Machine.Clear;
  end Test_Multprec_Trace_Interpolation1;

  procedure Test_Multprec_Trace_Interpolation
              ( file : in file_type;
                ep : in Standard_Complex_Poly_Systems.Poly_Sys;
                mp : in Multprec_Complex_Poly_Systems.Poly_Sys;
                sols : in Standard_Complex_Solutions.Solution_List;
                dim,size : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test on the trace interpolation for general surfaces,
  --   with multi-precision arithmetic.

    use Multprec_Stacked_Sample_Grids;

    deg : constant natural32 := Standard_Complex_Solutions.Length_Of(sols);
    gsz : constant natural32 := Grid_Size(dim+1,deg);
    sli : Standard_Complex_VecVecs.VecVec(1..integer32(dim)) := Slices(ep,dim);
    sps : Standard_Sample_List := Create(sols,sli);
    testsps : Multprec_Sample_List;
    grid : Stacked_Sample_Grid(integer32(dim),integer32(deg));

  begin
    put("Number of Samples : "); put(gsz,1); new_line;
    new_line;
    Sampling_Machine.Initialize(ep,mp,integer32(dim),size);
    Sampling_Machine.Default_Tune_Sampler(2);
    Sampling_Machine.Interactive_Tune_Sampler;
    Sampling_Machine.Default_Tune_Refiner;
    Sampling_Machine.Interactive_Tune_Refiner(size);
    new_line;
    put_line("See the output file for results...");
    new_line;
    Multprec_Test_Samples(file,sps,sli,testsps);
    Multprec_Stacked_Grid_Creator(file,sps,true,size,grid);
    Multprec_Trace_Interpolation(file,integer32(deg),grid,testsps,output);
    Sampling_Machine.Clear;
  end Test_Multprec_Trace_Interpolation;

  procedure Main is

    outfile : file_type;
    lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    sols : Standard_Complex_Solutions.Solution_List;
    dim,deci,size : natural32 := 0;
    ans : character;
    output : boolean;

  begin
    new_line;
    put_line("Multivariate interpolation with trace formulas.");
    Standard_Read_Embedding(lp,sols,dim);
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Give the number of decimal places (<= 16 is standard) : ");
    get(deci);
    put("Do you wish intermediate output during creation ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    if deci > 16 then
      size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
      Get_Multprec_System(lp.all,mp,size,dim);
      if dim = 1 then
        Test_Multprec_Trace_Interpolation1
          (outfile,lp.all,mp.all,sols,dim,size,output);
      else
        Test_Multprec_Trace_Interpolation
          (outfile,lp.all,mp.all,sols,dim,size,output);
      end if;
    elsif dim = 1 then
      Test_Standard_Trace_Interpolation1(outfile,lp.all,sols,dim,output);
    else
      Test_Standard_Trace_Interpolation(outfile,lp.all,sols,dim,output);
    end if;
  end Main;

begin
  Main;
end ts_traces;
