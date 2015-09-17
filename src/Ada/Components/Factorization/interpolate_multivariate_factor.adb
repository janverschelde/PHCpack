with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
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
with Sample_Point_Lists_io;              use Sample_Point_Lists_io;
with DoblDobl_Sample_Lists;              use DoblDobl_Sample_Lists;
with QuadDobl_Sample_Lists;              use QuadDobl_Sample_Lists;
with Hypersurface_Sample_Grids;          use Hypersurface_Sample_Grids;
with DoblDobl_Gridded_Hypersurfaces;     use DoblDobl_Gridded_Hypersurfaces;
with QuadDobl_Gridded_Hypersurfaces;     use QuadDobl_Gridded_Hypersurfaces;
with Standard_Stacked_Sample_Grids;
with DoblDobl_Stacked_Sample_Grids; 
with QuadDobl_Stacked_Sample_Grids;
with Standard_Divided_Differences;
with DoblDobl_Divided_Differences;
with QuadDobl_Divided_Differences;
with Standard_Trace_Interpolators;
with DoblDobl_Trace_Interpolators;
with QuadDobl_Trace_Interpolators;

package body Interpolate_Multivariate_Factor is

-- AUXILIARY FUNCTIONS for Normalize :

  function Leading_Coefficient
              ( p : Standard_Complex_Polynomials.Poly;
                tol : double_float )
              return Standard_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the first coefficient in p larger than tol.

    use Standard_Complex_Numbers;
    res : Complex_Number;

    procedure Scan_Term
                ( t : in Standard_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      if AbsVal(t.cf) > tol then
        res := t.cf;
        continue := false;
      else
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is
      new Standard_Complex_Polynomials.Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Leading_Coefficient;

  function Leading_Coefficient
              ( p : DoblDobl_Complex_Polynomials.Poly;
                tol : double_float )
              return DoblDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the first coefficient in p larger than tol.

    use DoblDobl_Complex_Numbers;
    res : Complex_Number;

    procedure Scan_Term
                ( t : in DoblDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is

      val : constant double_double := AbsVal(t.cf);
 
    begin
      if val > tol then
        res := t.cf;
        continue := false;
      else
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is
      new DoblDobl_Complex_Polynomials.Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Leading_Coefficient;

  function Leading_Coefficient
              ( p : QuadDobl_Complex_Polynomials.Poly;
                tol : double_float )
              return QuadDobl_Complex_Numbers.Complex_Number is

  -- DESCRIPTION :
  --   Returns the first coefficient in p larger than tol.

    use QuadDobl_Complex_Numbers;
    res : Complex_Number;

    procedure Scan_Term
                ( t : in QuadDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is

      val : constant quad_double := AbsVal(t.cf);
 
    begin
      if val > tol then
        res := t.cf;
        continue := false;
      else
        continue := true;
      end if;
    end Scan_Term;
    procedure Scan_Terms is
      new QuadDobl_Complex_Polynomials.Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Leading_Coefficient;

-- TARGET ROUTINES :

  procedure Normalize ( p : in out Standard_Complex_Polynomials.Poly ) is

    use Standard_Complex_Numbers;

    tol : constant double_float := 1.0E-10;
    leadcff : constant Complex_Number := Leading_Coefficient(p,tol);

    procedure Normalize_Term
                ( t : in out Standard_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      t.cf := t.cf/leadcff;
      continue := true;
    end Normalize_Term;
    procedure Normalize_Terms is
      new Standard_Complex_Polynomials.Changing_Iterator(Normalize_Term);

  begin
    Normalize_Terms(p);
  end Normalize;

  procedure Normalize ( p : in out DoblDobl_Complex_Polynomials.Poly ) is

    use DoblDobl_Complex_Numbers;

    tol : constant double_float := 1.0E-10;
    leadcff : constant Complex_Number := Leading_Coefficient(p,tol);

    procedure Normalize_Term
                ( t : in out DoblDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      t.cf := t.cf/leadcff;
      continue := true;
    end Normalize_Term;
    procedure Normalize_Terms is
      new DoblDobl_Complex_Polynomials.Changing_Iterator(Normalize_Term);

  begin
    Normalize_Terms(p);
  end Normalize;

  procedure Normalize ( p : in out QuadDobl_Complex_Polynomials.Poly ) is

    use QuadDobl_Complex_Numbers;

    tol : constant double_float := 1.0E-10;
    leadcff : constant Complex_Number := Leading_Coefficient(p,tol);

    procedure Normalize_Term
                ( t : in out QuadDobl_Complex_Polynomials.Term;
                  continue : out boolean ) is
    begin
      t.cf := t.cf/leadcff;
      continue := true;
    end Normalize_Term;
    procedure Normalize_Terms is
      new QuadDobl_Complex_Polynomials.Changing_Iterator(Normalize_Term);

  begin
    Normalize_Terms(p);
  end Normalize;

  procedure Normalize ( p : in out Standard_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Normalize(p(i));
    end loop;
  end Normalize;

  procedure Normalize ( p : in out DoblDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Normalize(p(i));
    end loop;
  end Normalize;

  procedure Normalize ( p : in out QuadDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    for i in p'range loop
      Normalize(p(i));
    end loop;
  end Normalize;

  function Interpolate_Factor
              ( p : Standard_Complex_Polynomials.Poly;
                b,v,w : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Stacked_Sample_Grids;
    use Standard_Divided_Differences;
    use Standard_Trace_Interpolators;

    res : Standard_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;

  begin
    Hypersurface_Sample_Grids.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_Standard_Sample_Lists(0..d)
             := Parallel_Sample(b,v,w,d);
        q1 : Newton_Interpolator1 := Create(grid,Create(1.0));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      declare
        grid : constant Stacked_Sample_Grid(n-1,d) := Full_Sample(b,v,w);
        t : constant Trace_Interpolator := Create(grid,d);
      begin
        res := Expand(t);
      end;
    end if;
    Hypersurface_Sample_Grids.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

  function Interpolate_Factor
              ( p : DoblDobl_Complex_Polynomials.Poly;
                b,v,w : DoblDobl_Complex_Vectors.Vector )
              return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Stacked_Sample_Grids;
    use DoblDobl_Divided_Differences;
    use DoblDobl_Trace_Interpolators;

    res : DoblDobl_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;

  begin
    DoblDobl_Gridded_Hypersurfaces.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_DoblDobl_Sample_Lists(0..d)
             := Parallel_Sample(b,v,w,d);
        one : constant double_double := create(1.0);
        q1 : Newton_Interpolator1 := Create(grid,Create(one));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      declare
        grid : constant Stacked_Sample_Grid(n-1,d) := Full_Sample(b,v,w);
        t : constant Trace_Interpolator := Create(grid,d);
      begin
        res := Expand(t);
      end;
    end if;
    DoblDobl_Gridded_Hypersurfaces.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

  function Interpolate_Factor
              ( p : QuadDobl_Complex_Polynomials.Poly;
                b,v,w : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Stacked_Sample_Grids;
    use QuadDobl_Divided_Differences;
    use QuadDobl_Trace_Interpolators;

    res : QuadDobl_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;

  begin
    QuadDobl_Gridded_Hypersurfaces.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_QuadDobl_Sample_Lists(0..d)
             := Parallel_Sample(b,v,w,d);
        one : constant quad_double := create(1.0);
        q1 : Newton_Interpolator1 := Create(grid,Create(one));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      declare
        grid : constant Stacked_Sample_Grid(n-1,d) := Full_Sample(b,v,w);
        t : constant Trace_Interpolator := Create(grid,d);
      begin
        res := Expand(t);
      end;
    end if;
    QuadDobl_Gridded_Hypersurfaces.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

  function Interpolate_Factor
              ( file : file_type;
                p : Standard_Complex_Polynomials.Poly;
                b,v,w : Standard_Complex_Vectors.Vector )
              return Standard_Complex_Polynomials.Poly is

    use Standard_Complex_Numbers;
    use Standard_Divided_Differences;
    use Standard_Stacked_Sample_Grids;
    use Standard_Trace_Interpolators;

    res : Standard_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;
    max_err : double_float;

  begin
    Hypersurface_Sample_Grids.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_Standard_Sample_Lists(0..d)
             := Parallel_Sample(file,false,b,v,w,d);
        q1 : Newton_Interpolator1 := Create(grid,Create(1.0));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        max_err := Maximal_Error(q1,grid);
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      new_line(file);
      put(file,"Computing ");
      put(file,Full_Grid_Size(natural32(n),natural32(d)),1);
      put_line(file," samples...");
      flush(file);
      declare
        grid : constant Stacked_Sample_Grid(n-1,d)
             := Full_Sample(file,b,v,w);
        t : constant Trace_Interpolator := Create(file,grid,d);
      begin
        max_err := Maximal_Error(t,grid);
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(t);
      end;
    end if;
    Hypersurface_Sample_Grids.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

  function Interpolate_Factor
              ( file : file_type;
                p : DoblDobl_Complex_Polynomials.Poly;
                b,v,w : DoblDobl_Complex_Vectors.Vector )
              return DoblDobl_Complex_Polynomials.Poly is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Stacked_Sample_Grids;
    use DoblDobl_Divided_Differences;
    use DoblDobl_Trace_Interpolators;

    res : DoblDobl_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;
    max_err : double_float;

  begin
    DoblDobl_Gridded_Hypersurfaces.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_DoblDobl_Sample_Lists(0..d)
             := Parallel_Sample(file,false,b,v,w,d);
        one : constant double_double := create(1.0);
        q1 : Newton_Interpolator1 := Create(grid,Create(one));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        max_err := Maximal_Error(q1,grid);
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      new_line(file);
      put(file,"Computing ");
      put(file,Full_Grid_Size(natural32(n),natural32(d)),1);
      put_line(file," samples...");
      flush(file);
      declare
        grid : constant Stacked_Sample_Grid(n-1,d)
             := Full_Sample(file,b,v,w);
        t : constant Trace_Interpolator := Create(file,grid,d);
      begin
        max_err := hi_part(Maximal_Error(t,grid));
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(t);
      end;
    end if;
    DoblDobl_Gridded_Hypersurfaces.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

  function Interpolate_Factor
              ( file : file_type;
                p : QuadDobl_Complex_Polynomials.Poly;
                b,v,w : QuadDobl_Complex_Vectors.Vector )
              return QuadDobl_Complex_Polynomials.Poly is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Stacked_Sample_Grids;
    use QuadDobl_Divided_Differences;
    use QuadDobl_Trace_Interpolators;

    res : QuadDobl_Complex_Polynomials.Poly;
    n : constant integer32 := b'length;
    d : constant integer32 := w'length;
    max_err : double_float;

  begin
    QuadDobl_Gridded_Hypersurfaces.Initialize(p);
    if n = 2 then
      declare
        grid : Array_of_QuadDobl_Sample_Lists(0..d)
             := Parallel_Sample(file,false,b,v,w,d);
        one : constant quad_double := create(1.0);
        q1 : Newton_Interpolator1 := Create(grid,Create(one));
        eq : Newton_Form_Evaluator1 := Create(q1);
      begin
        max_err := Maximal_Error(q1,grid);
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(eq);
        Deep_Clear(grid);
        Clear(q1); Clear(eq);
      end;
    else
      new_line(file);
      put(file,"Computing ");
      put(file,Full_Grid_Size(natural32(n),natural32(d)),1);
      put_line(file," samples...");
      flush(file);
      declare
        grid : constant Stacked_Sample_Grid(n-1,d)
             := Full_Sample(file,b,v,w);
        t : constant Trace_Interpolator := Create(file,grid,d);
      begin
        max_err := hihi_part(Maximal_Error(t,grid));
        put(file,"The maximal error on the grid : ");
        put(file,max_err,3); new_line(file);
        res := Expand(t);
      end;
    end if;
    QuadDobl_Gridded_Hypersurfaces.Clear;
    Normalize(res);
    return res;
  end Interpolate_Factor;

end Interpolate_Multivariate_Factor;
