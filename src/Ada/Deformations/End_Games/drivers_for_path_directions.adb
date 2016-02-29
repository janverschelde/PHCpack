with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Floating_Vectors_io;       use Standard_Floating_Vectors_io;
with Double_Double_Vectors_io;           use Double_Double_Vectors_io;
with Quad_Double_Vectors_io;             use Quad_Double_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vector_Norms;      use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vector_Norms;      use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Matrices;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Standard_IncFix_Continuation;
with DoblDobl_IncFix_Continuation;
with QuadDobl_IncFix_Continuation;

package body Drivers_for_Path_Directions is

-- AUXILIARIES TO INSTANTIATE :

  function Toric_Evaluator ( x : Standard_Complex_Vectors.Vector;
                             t : Standard_Complex_Numbers.Complex_Number )
                           return Standard_Complex_Vectors.Vector is
  begin
    return Standard_Homotopy.Eval(x,t);
  end Toric_Evaluator;

  function Toric_Evaluator ( x : DoblDobl_Complex_Vectors.Vector;
                             t : DoblDobl_Complex_Numbers.Complex_Number )
                           return DoblDobl_Complex_Vectors.Vector is
  begin
    return DoblDobl_Homotopy.Eval(x,t);
  end Toric_Evaluator;

  function Toric_Evaluator ( x : QuadDobl_Complex_Vectors.Vector;
                             t : QuadDobl_Complex_Numbers.Complex_Number )
                           return QuadDobl_Complex_Vectors.Vector is
  begin
    return QuadDobl_Homotopy.Eval(x,t);
  end Toric_Evaluator;

  function Toric_Differentiator
               ( x : Standard_Complex_Vectors.Vector;
                 t : Standard_Complex_Numbers.Complex_Number )
               return Standard_Complex_Vectors.Vector is
  begin
    return Standard_Homotopy.Eval(x,t);
  end Toric_Differentiator;

  function Toric_Differentiator
               ( x : DoblDobl_Complex_Vectors.Vector;
                 t : DoblDobl_Complex_Numbers.Complex_Number )
               return DoblDobl_Complex_Vectors.Vector is
  begin
    return DoblDobl_Homotopy.Eval(x,t);
  end Toric_Differentiator;

  function Toric_Differentiator
               ( x : QuadDobl_Complex_Vectors.Vector;
                 t : QuadDobl_Complex_Numbers.Complex_Number )
               return QuadDobl_Complex_Vectors.Vector is
  begin
    return QuadDobl_Homotopy.Eval(x,t);
  end Toric_Differentiator;

  function Toric_Differentiator
               ( x : Standard_Complex_Vectors.Vector;
                 t : Standard_Complex_Numbers.Complex_Number )
               return Standard_Complex_Matrices.Matrix is
  begin
    return Standard_Homotopy.Diff(x,t);
  end Toric_Differentiator;

  function Toric_Differentiator
               ( x : DoblDobl_Complex_Vectors.Vector;
                 t : DoblDobl_Complex_Numbers.Complex_Number )
               return DoblDobl_Complex_Matrices.Matrix is
  begin
    return DoblDobl_Homotopy.Diff(x,t);
  end Toric_Differentiator;

  function Toric_Differentiator
               ( x : QuadDobl_Complex_Vectors.Vector;
                 t : QuadDobl_Complex_Numbers.Complex_Number )
               return QuadDobl_Complex_Matrices.Matrix is
  begin
    return QuadDobl_Homotopy.Diff(x,t);
  end Toric_Differentiator;

-- TARGET ROUTINES :

  procedure Init_Path_Directions
               ( n,nv : in natural32;
                 v : in out Standard_Floating_VecVecs.Link_to_VecVec;
                 errv : in out Standard_Floating_Vectors.Link_to_Vector ) is

  begin
    v := new Standard_Floating_VecVecs.VecVec(1..integer32(nv));
    for i in v'range loop
      v(i) := new Standard_Floating_Vectors.Vector'(1..integer32(n) => 0.0);
    end loop;
    errv := new Standard_Floating_Vectors.Vector'(1..integer32(nv) => 1.0);
  end Init_Path_Directions;

  procedure Init_Path_Directions
               ( n,nv : in natural32;
                 v : in out Double_Double_VecVecs.Link_to_VecVec;
                 errv : in out Double_Double_Vectors.Link_to_Vector ) is

    zero : constant double_double := create(integer32(0));
    one : constant double_double := create(integer32(1));

  begin
    v := new Double_Double_VecVecs.VecVec(1..integer32(nv));
    for i in v'range loop
      v(i) := new Double_Double_Vectors.Vector'(1..integer32(n) => zero);
    end loop;
    errv := new Double_Double_Vectors.Vector'(1..integer32(nv) => one);
  end Init_Path_Directions;

  procedure Init_Path_Directions
               ( n,nv : in natural32;
                 v : in out Quad_Double_VecVecs.Link_to_VecVec;
                 errv : in out Quad_Double_Vectors.Link_to_Vector ) is

    zero : constant quad_double := create(integer32(0));
    one : constant quad_double := create(integer32(1));

  begin
    v := new Quad_Double_VecVecs.VecVec(1..integer32(nv));
    for i in v'range loop
      v(i) := new Quad_Double_Vectors.Vector'(1..integer32(n) => zero);
    end loop;
    errv := new Quad_Double_Vectors.Vector'(1..integer32(nv) => one);
  end Init_Path_Directions;

  procedure Toric_Continue
               ( file : in file_type;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 proj,report : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Standard_Floating_VecVecs.VecVec;
                 errv : in out Standard_Floating_Vectors.Vector;
                 target : in Standard_Complex_Numbers.Complex_Number ) is

    use Standard_IncFix_Continuation;

    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Toric_Continue(Max_Norm,Toric_Evaluator,
                                Toric_Differentiator,Toric_Differentiator);
    procedure Rep_Cont is
      new Reporting_Toric_Continue(Max_Norm,Toric_Evaluator,
                                   Toric_Differentiator,Toric_Differentiator);
  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,w,v,errv,target=>target);
     else Sil_Cont(sols,proj,w,v,errv,target=>target);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"toric continuation"); new_line(file);
  end Toric_Continue;

  procedure Toric_Continue
               ( file : in file_type;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 proj,report : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Double_Double_VecVecs.VecVec;
                 errv : in out Double_Double_Vectors.Vector;
                 target : in DoblDobl_Complex_Numbers.Complex_Number ) is

    use DoblDobl_IncFix_Continuation;

    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Toric_Continue(Max_Norm,Toric_Evaluator,
                                Toric_Differentiator,Toric_Differentiator);
    procedure Rep_Cont is
      new Reporting_Toric_Continue(Max_Norm,Toric_Evaluator,
                                   Toric_Differentiator,Toric_Differentiator);
  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,w,v,errv,target=>target);
     else Sil_Cont(sols,proj,w,v,errv,target=>target);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"toric continuation"); new_line(file);
  end Toric_Continue;

  procedure Toric_Continue
               ( file : in file_type;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 proj,report : in boolean;
                 w : in out Standard_Integer_Vectors.Vector;
                 v : in out Quad_Double_VecVecs.VecVec;
                 errv : in out Quad_Double_Vectors.Vector;
                 target : in QuadDobl_Complex_Numbers.Complex_Number ) is

    use QuadDobl_IncFix_Continuation;

    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Toric_Continue(Max_Norm,Toric_Evaluator,
                                Toric_Differentiator,Toric_Differentiator);
    procedure Rep_Cont is
      new Reporting_Toric_Continue(Max_Norm,Toric_Evaluator,
                                   Toric_Differentiator,Toric_Differentiator);
  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,w,v,errv,target=>target);
     else Sil_Cont(sols,proj,w,v,errv,target=>target);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"toric continuation"); new_line(file);
  end Toric_Continue;

  procedure Write_Direction
               ( file : in file_type; w : in integer32;
                 v : in Standard_Floating_Vectors.Vector;
                 error : in double_float; i : in integer32 ) is
  begin
    put(file,"Computed direction of path ");
    put(file,i,1); put_line(file," :"); put_line(file,v);
    put(file,"with estimated winding number : ");
    put(file,w,1); new_line(file);
    put(file,"and error : "); put(file,error); new_line(file);
  end Write_Direction;

  procedure Write_Direction
               ( file : in file_type; w : in integer32;
                 v : in Double_Double_Vectors.Vector;
                 error : in double_double; i : in integer32 ) is
  begin
    put(file,"Computed direction of path ");
    put(file,i,1); put_line(file," :"); put_line(file,v);
    put(file,"with estimated winding number : ");
    put(file,w,1); new_line(file);
    put(file,"and error : "); put(file,error); new_line(file);
  end Write_Direction;

  procedure Write_Direction
               ( file : in file_type; w : in integer32;
                 v : in Quad_Double_Vectors.Vector;
                 error : in quad_double; i : in integer32 ) is
  begin
    put(file,"Computed direction of path ");
    put(file,i,1); put_line(file," :"); put_line(file,v);
    put(file,"with estimated winding number : ");
    put(file,w,1); new_line(file);
    put(file,"and error : "); put(file,error); new_line(file);
  end Write_Direction;

  procedure Write_Directions 
               ( file : in file_type;
                 w : in Standard_Integer_Vectors.Vector;
                 v : in Standard_Floating_VecVecs.VecVec;
                 errv : in Standard_Floating_Vectors.Vector ) is
  begin
    for i in v'range loop
      Write_Direction(file,w(i),v(i).all,errv(i),i);
    end loop;
  end Write_Directions;

  procedure Write_Directions 
               ( file : in file_type;
                 w : in Standard_Integer_Vectors.Vector;
                 v : in Double_Double_VecVecs.VecVec;
                 errv : in Double_Double_Vectors.Vector ) is
  begin
    for i in v'range loop
      Write_Direction(file,w(i),v(i).all,errv(i),i);
    end loop;
  end Write_Directions;

  procedure Write_Directions 
               ( file : in file_type;
                 w : in Standard_Integer_Vectors.Vector;
                 v : in Quad_Double_VecVecs.VecVec;
                 errv : in Quad_Double_Vectors.Vector ) is
  begin
    for i in v'range loop
      Write_Direction(file,w(i),v(i).all,errv(i),i);
    end loop;
  end Write_Directions;

end Drivers_for_Path_Directions;
