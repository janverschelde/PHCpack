with text_io;                            use text_io;
with duration_io;
with Ada.Calendar;
with Communications_with_User;           use Communications_with_User;
with Time_Stamps;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_VecMats;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_VecMats;
with Standard_Complex_Series_Vectors;
with Standard_Random_Series_Vectors;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Random_Series_Vectors;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Random_Series_Vectors;
with Series_Coefficient_Vectors;         use Series_Coefficient_Vectors;
with Standard_Speelpenning_Convolutions;
with DoblDobl_Speelpenning_Convolutions;
with QuadDobl_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;
with Evaluation_Differentiation_Errors;  use Evaluation_Differentiation_Errors;
with Multitasked_AlgoDiff_Convolutions;  use Multitasked_AlgoDiff_Convolutions;

procedure ts_mtadcnv is

-- DESCRIPTION :
--   Development of algorithmic differentiation to evaluate and differentiate
--   polynomial systems at truncated power series with multitasking.

  procedure Standard_Random_Test
              ( dim,deg,nbr,pwr,nbt : in integer32;
                output : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use Standard_Speelpenning_Convolutions;

    c : constant Convolution_Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    x : constant Standard_Complex_Series_Vectors.Vector(1..dim)
      := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant Standard_Complex_VecVecs.VecVec(1..dim)
         := Standard_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant Standard_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant Standard_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant Standard_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : double_float;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put("running with "); put(nbt,1); put_line(" tasks ...");
    multstart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    Standard_Multitasked_EvalDiff(nbt,c,xcff,pwt,vy2,vm2,output);
    multstop := Ada.Calendar.Clock;
    mult_elapsed := multstop - multstart;
    err := Difference(vy1,vy2);
    put("the error of evaluation : "); put(err,3); new_line;
    err := Difference(vm1,vm2);
    put("the error of differentiation : "); put(err,3); new_line;
    new_line;
    put_line("Elapsed time on multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    new_line;
    put_line("Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    new_line;
    if seri_elapsed + 1.0 /= 1.0 then
      speedup := seri_elapsed/mult_elapsed;
      put("The speedup : ");
      duration_io.put(speedup,1,3); new_line;
    end if;
  end Standard_Random_Test;

  procedure DoblDobl_Random_Test
              ( dim,deg,nbr,pwr,nbt : in integer32;
                output : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use DoblDobl_Speelpenning_Convolutions;

    c : constant Convolution_Circuits
      := DoblDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    x : constant DoblDobl_Complex_Series_Vectors.Vector(1..dim)
      := DoblDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant DoblDobl_Complex_VecVecs.VecVec(1..dim)
         := DoblDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant DoblDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant DoblDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant DoblDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : double_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put("running with "); put(nbt,1); put_line(" tasks ...");
    multstart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    DoblDobl_Multitasked_EvalDiff(nbt,c,xcff,pwt,vy2,vm2,output);
    multstop := Ada.Calendar.Clock;
    mult_elapsed := multstop - multstart;
    err := Difference(vy1,vy2);
    put("the error of evaluation : "); put(err,3); new_line;
    err := Difference(vm1,vm2);
    put("the error of differentiation : "); put(err,3); new_line;
    new_line;
    put_line("Elapsed time on multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    new_line;
    put_line("Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    new_line;
    if seri_elapsed + 1.0 /= 1.0 then
      speedup := seri_elapsed/mult_elapsed;
      put("The speedup : ");
      duration_io.put(speedup,1,3); new_line;
    end if;
  end DoblDobl_Random_Test;

  procedure QuadDobl_Random_Test
              ( dim,deg,nbr,pwr,nbt : in integer32;
                output : in boolean := true ) is

  -- DESCRIPTION :
  --   Generates a random convolution circuit and applies multitasking
  --   to evaluate and differentiate with multitasking,
  --   in double double precision.  Compares with the one task run.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   nbt      the number of tasks;
  --   output   flag for output during multitasking.

    use QuadDobl_Speelpenning_Convolutions;

    c : constant Convolution_Circuits
      := QuadDobl_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    x : constant QuadDobl_Complex_Series_Vectors.Vector(1..dim)
      := QuadDobl_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
    xcff : constant QuadDobl_Complex_VecVecs.VecVec(1..dim)
         := QuadDobl_Series_Coefficients(x);
    mxe : constant Standard_Integer_Vectors.Vector(1..dim)
        := Exponent_Maxima(c,dim);
    pwt : constant Link_to_VecVecVec := Allocate(mxe,deg);
    yd : constant QuadDobl_Complex_VecVecs.VecVec(1..dim+1)
       := Allocate_Coefficients(dim+1,deg);
    vy1 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vy2 : constant QuadDobl_Complex_VecVecs.VecVec(0..deg)
        := Linearized_Allocation(dim,deg);
    vm1 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    vm2 : constant QuadDobl_Complex_VecMats.VecMat(0..deg)
        := Allocate_Coefficients(dim,dim,deg);
    err : quad_double;
    multstart,multstop,seristart,seristop : Ada.Calendar.Time;
    seri_elapsed,mult_elapsed,speedup : Duration;

    use Ada.Calendar; -- for the difference operator on Duration

  begin
    put_line("running on one task ...");
    seristart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    EvalDiff(c,xcff,pwt,yd,vy1,vm1);
    seristop := Ada.Calendar.Clock;
    seri_elapsed := seristop - seristart;
    put("running with "); put(nbt,1); put_line(" tasks ...");
    multstart := Ada.Calendar.Clock;
    Compute(pwt,mxe,xcff);
    QuadDobl_Multitasked_EvalDiff(nbt,c,xcff,pwt,vy2,vm2,output);
    multstop := Ada.Calendar.Clock;
    mult_elapsed := multstop - multstart;
    err := Difference(vy1,vy2);
    put("the error of evaluation : "); put(err,3); new_line;
    err := Difference(vm1,vm2);
    put("the error of differentiation : "); put(err,3); new_line;
    new_line;
    put_line("Elapsed time on multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,multstart,multstop);
    new_line;
    put_line("Elapsed time without multitasking : ");
    Time_Stamps.Write_Elapsed_Time(standard_output,seristart,seristop);
    new_line;
    if seri_elapsed + 1.0 /= 1.0 then
      speedup := seri_elapsed/mult_elapsed;
      put("The speedup : ");
      duration_io.put(speedup,1,3); new_line;
    end if;
  end QuadDobl_Random_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of tasks,
  --   for the dimensions of the problem,
  --   generates then the data for the problem, and
  --   launches the tasks.

    dim,deg,nbr,pwr,nbt : integer32 := 0;
    answer,precision : character;

  begin
    new_line;
    put_line("MENU for the working precision :");
    put_line("  0. standard double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(precision,"012");
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the series : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    put("Give the highest power of each variable : "); get(pwr);
    put("Give the number of tasks : "); get(nbt);
    put("Output during multitasking ? (y/n) "); Ask_Yes_or_No(answer);
    new_line;
    case precision is
      when '0' => Standard_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
      when '1' => DoblDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
      when '2' => QuadDobl_Random_Test(dim,deg,nbr,pwr,nbt,answer = 'y');
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_mtadcnv;
