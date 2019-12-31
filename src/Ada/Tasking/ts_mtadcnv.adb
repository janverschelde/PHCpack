with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_VecMats;
with Standard_Complex_Series_Vectors;
with Standard_Random_Series_Vectors;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;
with Standard_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Multitasking;

procedure ts_mtadcnv is

-- DESCRIPTION :
--   Development of algorithmic differentiation to evaluate and differentiate
--   polynomial systems at truncated power series with multitasking.

  function Allocate_Work_Space
             ( nbt,dim,deg : integer32 )
             return Standard_Speelpenning_Convolutions.VecVecVec is

  -- DESCRIPTION :
  --   Returns work space for nbt tasks to evaluate circuits of
  --   dimension dim at power series of degree deg.

    use Standard_Speelpenning_Convolutions;

    res : VecVecVec(1..nbt);

  begin
    for i in res'range loop
      declare
        cff : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
            := Allocate_Coefficients(dim+1,deg);
      begin
        res(i) := new Standard_Complex_VecVecs.VecVec'(cff);
      end;
    end loop;
    return res;
  end Allocate_Work_Space;

  procedure Standard_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Speelpenning_Convolutions.Convolution_Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                pwt : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat;
                output : in boolean := false ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the convolution circuits in c at x
  --   with multitasking.

  -- ON ENTRY :
  --   nbt      number of tasks;
  --   c        convolution circuits;
  --   x        coefficients of power series;
  --   pwt      the power table of x, as high as the degrees in c;
  --   vy       allocated space for evaluated c at x;
  --   vm       allocated space for differentiated c at x;
  --   output   true for the writing of intermediate output.

    use Standard_Speelpenning_Convolutions;

    dim : constant integer32 := c'last; -- assuming square circuits
    deg : constant integer32 := vm'last;
    yd : constant VecVecVec(1..nbt) := Allocate_Work_Space(nbt,dim,deg);
    done : Multitasking.boolean_array(1..nbt) := (1..nbt => false);

    procedure Silent_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   without intermediate output.
 
    begin
      for k in c'range loop
        if (k mod n) = i-1 then
          EvalDiff(c(k).all,x,pwt,yd(k).all);
        end if;
      end loop;
      done(i) := true;
    end Silent_Job;
    procedure silent_do_jobs is new Multitasking.Silent_Workers(Silent_Job);

    procedure Report_Job ( i,n : integer32 ) is

    -- DESCRIPTION :
    --   Task i out of n will evaluate and differentiate the circuit
    --   with intermediate output.

      idx : integer32 := i;
      ydi : Standard_Complex_VecVecs.Link_to_VecVec;
      vleft,vright : Standard_Complex_Vectors.Link_to_Vector;
      mleft : Standard_Complex_Matrices.Link_to_Matrix;
 
    begin
      put_line("number of circuits : " & Multitasking.to_string(c'last));
      put_line("number of tasks : " & Multitasking.to_string(n));
      while idx <= c'last loop
        put_line("considering circuit " & Multitasking.to_string(idx));
        put_line("task " & Multitasking.to_string(i)
                         & " computes circuit "
                         & Multitasking.to_string(idx));
        EvalDiff(c(idx).all,x,pwt,yd(i).all);
        put_line("task " & Multitasking.to_string(i)
                         & " done computing circuit "
                         & Multitasking.to_string(idx));
        vleft := vy(idx);
        ydi := yd(i);
        vright := ydi(x'last+1);
        for j in vleft'range loop
          vleft(j) := vright(j);
          vright(j) := Standard_Complex_Numbers.Create(0.0);
        end loop;
        for j in 1..x'last loop
          vright := ydi(j);
          for k in vm'range loop       -- k-th coefficient in matrix vm(k)
            mleft := vm(k);            -- row idx in vm(k) is the equation
            mleft(idx,j) := vright(k); -- column j in vm(k) is the variable
            vright(k) := Standard_Complex_Numbers.Create(0.0);
          end loop;
        end loop;
        put_line("idx before increment : " & Multitasking.to_string(idx));
        idx := idx + n;
        put_line("idx after increment : " & Multitasking.to_string(idx));
      end loop;
      put_line("task " & Multitasking.to_string(i) & " is done.");
      done(i) := true;
    end Report_Job;
    procedure report_do_jobs is new Multitasking.Silent_Workers(Report_Job);

  begin
    put("dim : "); put(dim,1); new_line;
    if output
     then report_do_jobs(nbt);
     else silent_do_jobs(nbt);
    end if;
   -- make sure main task does not terminate before all worker tasks
    while not Multitasking.all_true(nbt,done) loop
      delay 0.1;
    end loop;
  end Standard_Multitasked_EvalDiff;

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr,nbt : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks.

    use Standard_Speelpenning_Convolutions;

    c : constant Convolution_Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Create(xcff,mxe);
    yd : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Allocate_Coefficients(dim,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(1..dim)
        := Allocate_Coefficients(dim,deg);
    vm1 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    timer : timing_widget;
    err : double_float;

  begin
   -- if nbt = 1 then
      tstart(timer);
      EvalDiff(c,xcff,pwt,yd,vy1,vm1);
      tstop(timer);
   -- else
      tstart(timer);
      Standard_Multitasked_EvalDiff(nbt,c,xcff,pwt,vy2,vm2,true);
      tstop(timer);
   -- end if;
    err := Difference(vy1,vy2);
    put("the error of evaluation : "); put(err,3); new_line;
    err := Difference(vm1,vm2);
    put("the error of differentiation : "); put(err,3); new_line;
   -- new_line;
   -- print_times(standard_output,timer,"evaluation and differentiation");
  end Standard_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of tasks,
  --   for the dimensions of the problem,
  --   generates then the data for the problem, and
  --   launches the tasks.

    dim,deg,nbr,pwr,nbt : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the series : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    put("Give the highest power of each variable : "); get(pwr);
    put("Give the number of tasks : "); get(nbt);
    new_line;
    Standard_Random_Test(dim,deg,nbr,pwr,nbt);
  end Main;

begin
  Main;
end ts_mtadcnv;
