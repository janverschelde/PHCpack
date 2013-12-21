with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Norms_Equals;      use Standard_Complex_Norms_Equals;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;
with Standard_Homotopy;
with Standard_IncFix_Continuation;       use Standard_IncFix_Continuation;

package body Drivers_for_Path_Directions is

-- AUXILIARIES TO INSTANTIATE :

 -- function Scale ( v : Vector ) return Vector is
 --
  -- DESCRIPTION :
  --   Divides the vector by its largest component.
 --
 --   res : Vector(v'range);
 --   tol : constant double_float := 10.0**(-12);
 --   ind : integer := v'first;
 --   max : double_float := abs(v(ind));
--
 -- begin
 --   for i in v'range loop
 --     if abs(v(i)) > max
 --      then max := abs(v(i));
 --           ind := i;
 --     end if;
 --   end loop;
 --   if max > tol
 --    then for i in v'range loop
 --           res(i) := v(i)/max;
 --         end loop;
 --   end if;
 --   return res;
 -- end Scale;

  function Toric_Evaluator ( x : Standard_Complex_Vectors.Vector;
                             t : Complex_Number )
                           return Standard_Complex_Vectors.Vector is
  begin
    return Standard_Homotopy.Eval(x,t);
  end Toric_Evaluator;

  function Toric_Differentiator
               ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
               return Standard_Complex_Vectors.Vector is
  begin
    return Standard_Homotopy.Eval(x,t);
  end Toric_Differentiator;

  function Toric_Differentiator
               ( x : Standard_Complex_Vectors.Vector; t : Complex_Number )
               return Standard_Complex_Matrices.Matrix is
  begin
    return Standard_Homotopy.Diff(x,t);
  end Toric_Differentiator;

-- TARGET ROUTINES :

  procedure Init_Path_Directions
               ( n,nv : in natural32; v : in out Link_to_VecVec;
                 errv : in out Link_to_Vector ) is

  begin
    v := new VecVec(1..integer32(nv));
    for i in v'range loop
      v(i) := new Vector'(1..integer32(n) => 0.0);
    end loop;
    errv := new Vector'(1..integer32(nv) => 1.0);
  end Init_Path_Directions;

  procedure Toric_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj,report : in boolean; v : in out VecVec;
                 errv : in out Vector; target : in Complex_Number ) is

  -- DESCRIPTION :
  --   Performs the continuation with online toric compactifications.

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
     then Rep_Cont(file,sols,proj,v,errv,target);
     else Sil_Cont(sols,proj,v,errv,target);
    end if;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"toric continuation"); new_line(file);
  end Toric_Continue;

  procedure Write_Directions 
               ( file : in file_type; v : in VecVec; errv : in Vector ) is

    procedure Write ( file : in file_type; v : in Vector ) is
    begin
      for i in v'range loop
        put(file,v(i)); new_line(file);
      end loop;
    end Write;

    procedure Write_Direction
                 ( file : in file_type;
                   v : in Vector; error : in double_float; i : integer32 ) is
    begin
      put(file,"Computed direction of path ");
      put(file,i,1); put_line(file," :"); Write(file,v);
      put(file,"with error : "); put(file,error); new_line(file);
    end Write_Direction;

  begin
    for i in v'range loop
      Write_Direction(file,v(i).all,errv(i),i);
    end loop;
  end Write_Directions;

end Drivers_for_Path_Directions;
