with text_io;                            use text_io;
with Timing_Package;                    use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_VecMats;
with Standard_Complex_Series_Vectors;
with Standard_Random_Series_Vectors;
with Series_Polynomial_Gradients;        use Series_Polynomial_Gradients;
with Standard_Speelpenning_Convolutions;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;

procedure ts_mtadcnv is

-- DESCRIPTION :
--   Development of algorithmic differentiation to evaluate and differentiate
--   polynomial systems at truncated power series with multitasking.

  procedure Standard_Multitasked_EvalDiff
              ( nbt : in integer32;
                c : in Standard_Speelpenning_Convolutions.Convolution_Circuits;
                x : in Standard_Complex_VecVecs.VecVec;
                pwt : in Standard_Speelpenning_Convolutions.Link_to_VecVecVec;
                vy : in Standard_Complex_VecVecs.VecVec;
                vm : in Standard_Complex_VecMats.VecMat ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the convolution circuits in c at x
  --   with multitasking.

  -- ON ENTRY :
  --   nbt      number of tasks;
  --   c        convolution circuits;
  --   x        coefficients of power series;
  --   pwt      the power table of x, as high as the degrees in c;
  --   vy       allocated space for evaluated c at x;
  --   vm       allocated space for differentiated c at x.

  begin
    null;
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
    vy : constant Standard_Complex_VecVecs.VecVec(1..dim)
       := Allocate_Coefficients(dim,deg);
    vm : constant Standard_Complex_VecMats.VecMat(0..deg)
       := Allocate_Coefficients(dim,dim,deg);
    timer : timing_widget;

  begin
    if nbt = 1 then
      tstart(timer);
      EvalDiff(c,xcff,pwt,yd,vy,vm);
      tstop(timer);
      print_times(standard_output,timer,"evaluation and differentiation");
    else
      tstart(timer);
      Standard_Multitasked_EvalDiff(nbt,c,xcff,pwt,vy,vm);
      tstop(timer);
      print_times(standard_output,timer,"evaluation and differentiation");
    end if;
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
