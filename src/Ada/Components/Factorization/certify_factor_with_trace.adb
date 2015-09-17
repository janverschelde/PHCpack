with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Sample_Point_Lists;                 use Sample_Point_Lists;
with DoblDobl_Sample_Lists;              use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Lists;              use QuadDobl_Sample_Lists;
with Hypersurface_Sample_Grids;          use Hypersurface_Sample_Grids;
with DoblDobl_Gridded_Hypersurfaces;     use DoblDobl_Gridded_Hypersurfaces;
with QuadDobl_Gridded_Hypersurfaces;     use QuadDobl_Gridded_Hypersurfaces;
with Standard_Trace_Interpolators;       use Standard_Trace_Interpolators;
with DoblDobl_Trace_Interpolators;       use DoblDobl_Trace_Interpolators;
with QuadDobl_Trace_Interpolators;       use QuadDobl_Trace_Interpolators;

package body Certify_Factor_with_Trace is

  procedure Certify_Factor
               ( p : in Standard_Complex_Polynomials.Poly;
                 b,v,w : in Standard_Complex_Vectors.Vector;
                 testres : out double_float ) is

    grid : Array_of_Standard_Sample_Lists(0..2);
    use Standard_Complex_Numbers;

  begin
    Hypersurface_Sample_Grids.Initialize(p);
    grid := Parallel_Sample1(b,v,w,2);
    declare
      q : constant Standard_Complex_Vectors.Vector := Create(grid,1);
      val,eva : Complex_Number;
    begin
      Eval_Trace(q,1,grid(2),val,eva);
      testres := AbsVal(val-eva);
    end;
    Deep_Clear(grid);
    Hypersurface_Sample_Grids.Clear;
  end Certify_Factor;

  procedure Certify_Factor
               ( p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,w : in DoblDobl_Complex_Vectors.Vector;
                 testres : out double_float ) is

    grid : Array_of_DoblDobl_Sample_Lists(0..2);
    use DoblDobl_Complex_Numbers;

  begin
    DoblDobl_Gridded_Hypersurfaces.Initialize(p);
    grid := Parallel_Sample1(b,v,w,2);
    declare
      q : constant DoblDobl_Complex_Vectors.Vector := Create(grid,1);
      val,eva : Complex_Number;
      dd_testres : double_double;
    begin
      Eval_Trace(q,1,grid(2),val,eva);
      dd_testres := AbsVal(val-eva);
      testres := hi_part(dd_testres);
    end;
    Deep_Clear(grid);
    Hypersurface_Sample_Grids.Clear;
  end Certify_Factor;

  procedure Certify_Factor
               ( p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,w : in QuadDobl_Complex_Vectors.Vector;
                 testres : out double_float ) is

    grid : Array_of_QuadDobl_Sample_Lists(0..2);
    use QuadDobl_Complex_Numbers;

  begin
    QuadDobl_Gridded_Hypersurfaces.Initialize(p);
    grid := Parallel_Sample1(b,v,w,2);
    declare
      q : constant QuadDobl_Complex_Vectors.Vector := Create(grid,1);
      val,eva : Complex_Number;
      qd_testres : quad_double;
    begin
      Eval_Trace(q,1,grid(2),val,eva);
      qd_testres := AbsVal(val-eva);
      testres := hihi_part(qd_testres);
    end;
    Deep_Clear(grid);
    Hypersurface_Sample_Grids.Clear;
  end Certify_Factor;

  procedure Certify_Factor
               ( file : in file_type;
                 p : in Standard_Complex_Polynomials.Poly;
                 b,v,w : in Standard_Complex_Vectors.Vector;
                 testres : out double_float ) is

    grid : Array_of_Standard_Sample_Lists(0..2);
    use Standard_Complex_Numbers;

  begin
    Hypersurface_Sample_Grids.Initialize(p);
    grid := Parallel_Sample1(file,false,b,v,w,2);
    declare
      q : constant Standard_Complex_Vectors.Vector := Create(grid,1);
      val,eva : Complex_Number;
    begin
      Eval_Trace(q,1,grid(2),val,eva);
      put(file,"Value at trace : "); put(file,val); new_line(file);
      put(file,"Computed sum   : "); put(file,eva); new_line(file);
      testres := AbsVal(val-eva);
      put(file,"Absolute Value of Difference : "); put(file,testres);
      new_line(file);
    end;
    Deep_Clear(grid);
    Hypersurface_Sample_Grids.Clear;
  end Certify_Factor;

  procedure Certify_Factor
               ( file : in file_type;
                 p : in DoblDobl_Complex_Polynomials.Poly;
                 b,v,w : in DoblDobl_Complex_Vectors.Vector;
                 testres : out double_float ) is

    grid : Array_of_DoblDobl_Sample_Lists(0..2);
    use DoblDobl_Complex_Numbers;

  begin
    DoblDobl_Gridded_Hypersurfaces.Initialize(p);
    grid := Parallel_Sample1(file,false,b,v,w,2);
    declare
      q : constant DoblDobl_Complex_Vectors.Vector := Create(grid,1);
      val,eva : Complex_Number;
      dd_testres : double_double;
    begin
      Eval_Trace(q,1,grid(2),val,eva);
      put(file,"Value at trace : "); put(file,val); new_line(file);
      put(file,"Computed sum   : "); put(file,eva); new_line(file);
      dd_testres := AbsVal(val-eva);
      testres := hi_part(dd_testres);
      put(file,"Absolute Value of Difference : "); put(file,testres);
      new_line(file);
    end;
    Deep_Clear(grid);
    Hypersurface_Sample_Grids.Clear;
  end Certify_Factor;

  procedure Certify_Factor
               ( file : in file_type;
                 p : in QuadDobl_Complex_Polynomials.Poly;
                 b,v,w : in QuadDobl_Complex_Vectors.Vector;
                 testres : out double_float ) is

    grid : Array_of_QuadDobl_Sample_Lists(0..2);
    use QuadDobl_Complex_Numbers;

  begin
    QuadDobl_Gridded_Hypersurfaces.Initialize(p);
    grid := Parallel_Sample1(file,false,b,v,w,2);
    declare
      q : constant QuadDobl_Complex_Vectors.Vector := Create(grid,1);
      val,eva : Complex_Number;
      qd_testres : quad_double;
    begin
      Eval_Trace(q,1,grid(2),val,eva);
      put(file,"Value at trace : "); put(file,val); new_line(file);
      put(file,"Computed sum   : "); put(file,eva); new_line(file);
      qd_testres := AbsVal(val-eva);
      testres := hihi_part(qd_testres);
      put(file,"Absolute Value of Difference : "); put(file,testres);
      new_line(file);
    end;
    Deep_Clear(grid);
    Hypersurface_Sample_Grids.Clear;
  end Certify_Factor;

end Certify_Factor_with_Trace;
