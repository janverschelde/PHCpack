with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_NesVecs_io;       use Standard_Complex_NesVecs_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_NesVecs_io;       use DoblDobl_Complex_NesVecs_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_NesVecs_io;       use QuadDobl_Complex_NesVecs_io;
with Standard_Floating_Matrices_io;     use Standard_Floating_Matrices_io;
with Double_Double_Matrices_io;         use Double_Double_Matrices_io;
with Quad_Double_Matrices_io;           use Quad_Double_Matrices_io;
with Multprec_Floating_Matrices_io;     use Multprec_Floating_Matrices_io;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems_io;  use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems_io;  use QuadDobl_Complex_Poly_Systems_io;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Systems;
with Standard_Complex_Solutions;        use Standard_Complex_Solutions;
with DoblDobl_Complex_Solutions;        use DoblDobl_Complex_Solutions;
with QuadDobl_Complex_Solutions;        use QuadDobl_Complex_Solutions;
with Sampling_Machine;
with DoblDobl_Sampling_Machine;
with QuadDobl_Sampling_Machine;
with Sample_Points;                     use Sample_Points;
with Sample_Point_Lists;                use Sample_Point_Lists;
with DoblDobl_Sample_Points;            use DoblDobl_Sample_Points;
with DoblDobl_Sample_Lists;             use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Points;            use QuadDobl_Sample_Points;
with QuadDobl_Sample_Lists;             use QuadDobl_Sample_Lists;
with Witness_Sets,Witness_Sets_io;      use Witness_Sets,Witness_Sets_io;
with Make_Sample_Grids;                 use Make_Sample_Grids;
with Standard_Stacked_Sample_Grids;
with DoblDobl_Stacked_Sample_Grids;
with QuadDobl_Stacked_Sample_Grids;
with Multprec_Stacked_Sample_Grids;
with Standard_Nvariate_Interpolators;
with DoblDobl_Nvariate_Interpolators;
with QuadDobl_Nvariate_Interpolators;
with Standard_Trace_Interpolators;
with DoblDobl_Trace_Interpolators;
with QuadDobl_Trace_Interpolators;
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

  procedure DoblDobl_Write_Evaluations
              ( file : in file_type;
                t : in DoblDobl_Trace_Interpolators.Trace_Interpolator1;
                sps : in DoblDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.
  --   This version is for the special case of dealing with a curve.

    use DoblDobl_Trace_Interpolators;
    tmp : DoblDobl_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end DoblDobl_Write_Evaluations;

  procedure QuadDobl_Write_Evaluations
              ( file : in file_type;
                t : in QuadDobl_Trace_Interpolators.Trace_Interpolator1;
                sps : in QuadDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.
  --   This version is for the special case of dealing with a curve.

    use QuadDobl_Trace_Interpolators;
    tmp : QuadDobl_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end QuadDobl_Write_Evaluations;

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

  procedure DoblDobl_Write_Evaluations
              ( file : in file_type;
                t : in DoblDobl_Trace_Interpolators.Trace_Interpolator;
                sps : in DoblDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.
  --   This version accepts the trace form for general surfaces.

    use DoblDobl_Trace_Interpolators;
    tmp : DoblDobl_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end DoblDobl_Write_Evaluations;

  procedure QuadDobl_Write_Evaluations
              ( file : in file_type;
                t : in QuadDobl_Trace_Interpolators.Trace_Interpolator;
                sps : in QuadDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Writes the result of the evaluation of t at all samples in sps
  --   on file, each time taking a new line for every evaluation.
  --   This version accepts the trace form for general surfaces.

    use QuadDobl_Trace_Interpolators;
    tmp : QuadDobl_Sample_List := sps;

  begin
    while not Is_Null(tmp) loop
      put(file,Eval(t,Sample_Point(Head_Of(tmp)).v)); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end QuadDobl_Write_Evaluations;

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

  procedure DoblDobl_Trace_Interpolation1
               ( file : in file_type;
                 grid : in Array_of_DoblDobl_Sample_Lists;
                 testsps : in DoblDobl_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.
  --   This routine is restricted to the special case of curves.

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Polynomials,DoblDobl_Trace_Interpolators;

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
    DoblDobl_Write_Evaluations(file,t,testsps);
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
  end DoblDobl_Trace_Interpolation1;

  procedure QuadDobl_Trace_Interpolation1
               ( file : in file_type;
                 grid : in Array_of_QuadDobl_Sample_Lists;
                 testsps : in QuadDobl_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.
  --   This routine is restricted to the special case of curves.

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Polynomials,QuadDobl_Trace_Interpolators;

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
    QuadDobl_Write_Evaluations(file,t,testsps);
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
  end QuadDobl_Trace_Interpolation1;

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

  procedure DoblDobl_Trace_Interpolation
              ( file : in file_type; d : in integer32;
                grid : in DoblDobl_Stacked_Sample_Grids.Stacked_Sample_Grid;
                testsps : in DoblDobl_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.

    use DoblDobl_Complex_Numbers,DoblDobl_Complex_Polynomials;
    use DoblDobl_Nvariate_Interpolators,DoblDobl_Trace_Interpolators;

    t : Trace_Interpolator;
    maxerr : double_double;
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
    DoblDobl_Write_Evaluations(file,t,testsps);
    put_line(file,"Test on construction of linear trace :");
    declare
      seltr : constant Newton_Form(natural32(grid.k),grid.k,1) := Trace(t,1);
      tr : Newton_Form(natural32(grid.k),grid.k,1);
      val,eva : Complex_Number;
      diff : double_double;
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
  end DoblDobl_Trace_Interpolation;

  procedure QuadDobl_Trace_Interpolation
              ( file : in file_type; d : in integer32;
                grid : in QuadDobl_Stacked_Sample_Grids.Stacked_Sample_Grid;
                testsps : in QuadDobl_Sample_List; output : in boolean ) is

  -- DESCRIPTION :
  --   Creates the trace form of the interpolating polynomial through
  --   the samples in the grid, computes the residuals at the grid,
  --   and at all samples in the list testsps.

    use QuadDobl_Complex_Numbers,QuadDobl_Complex_Polynomials;
    use QuadDobl_Nvariate_Interpolators,QuadDobl_Trace_Interpolators;

    t : Trace_Interpolator;
    maxerr : quad_double;
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
    QuadDobl_Write_Evaluations(file,t,testsps);
    put_line(file,"Test on construction of linear trace :");
    declare
      seltr : constant Newton_Form(natural32(grid.k),grid.k,1) := Trace(t,1);
      tr : Newton_Form(natural32(grid.k),grid.k,1);
      val,eva : Complex_Number;
      diff : quad_double;
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
  end QuadDobl_Trace_Interpolation;

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

  procedure Test_DoblDobl_Trace_Interpolation1
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test on the trace interpolation for curves with DoblDobl arithmetic.

    sli : constant DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : DoblDobl_Sample_List := Create(sols,sli);
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_DoblDobl_Sample_Lists(0..integer32(len));
    eps,dst : double_double;
    testsps : DoblDobl_Sample_List;

  begin
    DoblDobl_Sampling_Machine.Initialize(p);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0); -- 2 is too restrictive?
    DoblDobl_Sampling_Machine.Interactive_Tune_Sampler;
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    DoblDobl_Test_Samples(file,sps,sli,testsps);
    DoblDobl_Rectangular_Grid_Creator(file,sps,len,grid,eps,dst);
    DoblDobl_Trace_Interpolation1(file,grid,testsps,output);
    DoblDobl_Sampling_Machine.Clear;
  end Test_DoblDobl_Trace_Interpolation1;

  procedure Test_QuadDobl_Trace_Interpolation1
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test on the trace interpolation for curves with QuadDobl arithmetic.

    sli : constant QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Slices(p,dim);
    sps : QuadDobl_Sample_List := Create(sols,sli);
    len : constant natural32 := Length_Of(sps);
    grid : Array_of_QuadDobl_Sample_Lists(0..integer32(len));
    eps,dst : quad_double;
    testsps : QuadDobl_Sample_List;

  begin
    QuadDobl_Sampling_Machine.Initialize(p);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0); -- 2 too restrictive ?
    QuadDobl_Sampling_Machine.Interactive_Tune_Sampler;
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    new_line;
    put_line("See the output file for results...");
    new_line;
    QuadDobl_Test_Samples(file,sps,sli,testsps);
    QuadDobl_Rectangular_Grid_Creator(file,sps,len,grid,eps,dst);
    QuadDobl_Trace_Interpolation1(file,grid,testsps,output);
    QuadDobl_Sampling_Machine.Clear;
  end Test_QuadDobl_Trace_Interpolation1;

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

  procedure Test_DoblDobl_Trace_Interpolation
              ( file : in file_type;
                p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in DoblDobl_Complex_Solutions.Solution_List;
                dim : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test on the trace interpolation for general surfaces,
  --   with DoblDobl arithmetic.

    use DoblDobl_Stacked_Sample_Grids;

    deg : constant natural32 := DoblDobl_Complex_Solutions.Length_Of(sols);
    gsz : constant natural32 := Full_Grid_Size(dim+1,deg);
    sli : DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim)) := Slices(p,dim);
    sps : DoblDobl_Sample_List := Create(sols,sli);
    testsps : DoblDobl_Sample_List;
    grid : Stacked_Sample_Grid(integer32(dim),integer32(deg));

  begin
    put("Number of samples : "); put(gsz,1); new_line;
    new_line;
    DoblDobl_Sampling_Machine.Initialize(p);
    DoblDobl_Sampling_Machine.Default_Tune_Sampler(0); -- 2 too restrictive?
    DoblDobl_Sampling_Machine.Interactive_Tune_Sampler;
    DoblDobl_Sampling_Machine.Default_Tune_Refiner;
    new_line;
    new_line;
    put_line("See the output file for results...");
    new_line;
    DoblDobl_Test_Samples(file,sps,sli,testsps);
    DoblDobl_Stacked_Grid_Creator(file,sps,true,grid);
    DoblDobl_Trace_Interpolation(file,integer32(deg),grid,testsps,output);
    DoblDobl_Sampling_Machine.Clear;
  end Test_DoblDobl_Trace_Interpolation;

  procedure Test_QuadDobl_Trace_Interpolation
              ( file : in file_type;
                p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                sols : in QuadDobl_Complex_Solutions.Solution_List;
                dim : in natural32; output : in boolean ) is

  -- DESCRIPTION :
  --   Test on the trace interpolation for general surfaces,
  --   with QuadDobl arithmetic.

    use QuadDobl_Stacked_Sample_Grids;

    deg : constant natural32 := QuadDobl_Complex_Solutions.Length_Of(sols);
    gsz : constant natural32 := Full_Grid_Size(dim+1,deg);
    sli : QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim)) := Slices(p,dim);
    sps : QuadDobl_Sample_List := Create(sols,sli);
    testsps : QuadDobl_Sample_List;
    grid : Stacked_Sample_Grid(integer32(dim),integer32(deg));

  begin
    put("Number of samples : "); put(gsz,1); new_line;
    new_line;
    QuadDobl_Sampling_Machine.Initialize(p);
    QuadDobl_Sampling_Machine.Default_Tune_Sampler(0); -- 2 too restrictive ?
    QuadDobl_Sampling_Machine.Interactive_Tune_Sampler;
    QuadDobl_Sampling_Machine.Default_Tune_Refiner;
    new_line;
    new_line;
    put_line("See the output file for results...");
    new_line;
    QuadDobl_Test_Samples(file,sps,sli,testsps);
    QuadDobl_Stacked_Grid_Creator(file,sps,true,grid);
    QuadDobl_Trace_Interpolation(file,integer32(deg),grid,testsps,output);
    QuadDobl_Sampling_Machine.Clear;
  end Test_QuadDobl_Trace_Interpolation;

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
    st_lp : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
    dd_lp : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    qd_lp : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
    mp : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
    st_sols : Standard_Complex_Solutions.Solution_List;
    dd_sols : DoblDobl_Complex_Solutions.Solution_List;
    qd_sols : QuadDobl_Complex_Solutions.Solution_List;
    dim,deci,size : natural32 := 0;
    ans : character;
    output : boolean;

  begin
    new_line;
    put_line("Multivariate interpolation with trace formulas.");
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(outfile);
    new_line;
    put("Do you wish intermediate output during creation ? (y/n) ");
    Ask_Yes_or_No(ans);
    output := (ans = 'y');
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. standard double precision;");
    put_line("  1. double double precision;");
    put_line("  2. quad double precision;");
    put_line("  3. multiprecision;");
    put("Type 0, 1, 2, or 3 to select the precision : ");
    Ask_Alternative(ans,"0123");
    case ans is
      when '0' =>
        Standard_Read_Embedding(st_lp,st_sols,dim);
        if dim = 1 then
          Test_Standard_Trace_Interpolation1
            (outfile,st_lp.all,st_sols,dim,output);
        else
          Test_Standard_Trace_Interpolation
            (outfile,st_lp.all,st_sols,dim,output);
        end if;
      when '1' =>
        DoblDobl_Read_Embedding(dd_lp,dd_sols,dim);
        if dim = 1 then
          Test_DoblDobl_Trace_Interpolation1
            (outfile,dd_lp.all,dd_sols,dim,output);
        else
          Test_DoblDobl_Trace_Interpolation
            (outfile,dd_lp.all,dd_sols,dim,output);
        end if;
      when '2' =>
        QuadDobl_Read_Embedding(qd_lp,qd_sols,dim);
        if dim = 1 then
          Test_QuadDobl_Trace_Interpolation1
            (outfile,qd_lp.all,qd_sols,dim,output);
        else
          Test_QuadDobl_Trace_Interpolation
            (outfile,qd_lp.all,qd_sols,dim,output);
        end if;
      when '3' =>
        new_line;
        put("Give the number of decimal places : "); get(deci);
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        Get_Multprec_System(st_lp.all,mp,size,dim);
        if dim = 1 then
          Test_Multprec_Trace_Interpolation1
            (outfile,st_lp.all,mp.all,st_sols,dim,size,output);
        else
          Test_Multprec_Trace_Interpolation
            (outfile,st_lp.all,mp.all,st_sols,dim,size,output);
        end if;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_traces;
