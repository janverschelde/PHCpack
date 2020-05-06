with text_io;                            use text_io;
with Timing_Package;                     use Timing_Package;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_VecVecs;
with Standard_Complex_VecVecs;
with Standard_Vector_Splitters;
with Standard_Complex_Series_Vectors;
with Standard_Random_Series_Vectors;
with Series_Coefficient_Vectors;
with Standard_Speelpenning_Convolutions;
with Standard_Coefficient_Convolutions;
with Standard_Convolution_Splitters;
with Random_Convolution_Circuits;        use Random_Convolution_Circuits;

procedure ts_perfcirc is

-- DESCRIPTION :
--   Test on the performance of the evaluation and differentiation of
--   convolution circuits, once with complex coefficients, and once
--   with coefficients with separate real and imaginary parts.

  procedure Test ( dim,deg,nbr,pwr,frq : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random system of exponents and coefficients.
  --   Run tests in standard double precision.

  -- ON ENTRY :
  --   dim      dimension of the exponent vectors;
  --   deg      degree of the power series;
  --   nbr      number of products;
  --   pwr      largest power of the variables;
  --   frq      frequency of the tests.

    c : constant Standard_Speelpenning_Convolutions.Circuits
      := Standard_Random_Convolution_Circuits(dim,deg,nbr,pwr);
    s : constant Standard_Speelpenning_Convolutions.Link_to_System
      := Standard_Speelpenning_Convolutions.Create(c,dim,deg);
    cs : constant Standard_Coefficient_Convolutions.Link_to_System
       := Standard_Convolution_Splitters.Split(s);
    x : Standard_Complex_Series_Vectors.Vector(1..dim);
    xcff : Standard_Complex_VecVecs.VecVec(1..dim);
    rx : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    ix : constant Standard_Floating_VecVecs.Link_to_VecVec
       := Standard_Vector_Splitters.Allocate_Floating_Coefficients(dim,deg);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for k in 1..frq loop
      x := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
      xcff := Series_Coefficient_Vectors.Standard_Series_Coefficients(x);
      Standard_Speelpenning_Convolutions.Compute(s.pwt,s.mxe,xcff);
      Standard_Speelpenning_Convolutions.EvalDiff(s,xcff);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"complex circuits");
    tstart(timer);
    for k in 1..frq loop
      x := Standard_Random_Series_Vectors.Random_Series_Vector(1,dim,deg);
      xcff := Series_Coefficient_Vectors.Standard_Series_Coefficients(x);
      Standard_Vector_Splitters.Complex_Parts(xcff,rx,ix);
      Standard_Coefficient_Convolutions.Compute(cs.rpwt,cs.ipwt,cs.mxe,rx,ix);
      Standard_Coefficient_Convolutions.EvalDiff(cs,rx.all,ix.all);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"splitted circuits");
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of a system of random 
  --   convolution circuits and then launches the tests.

    dim,deg,nbr,pwr,frq : integer32 := 0;

  begin
    new_line;
    put("Give the dimension : "); get(dim);
    put("Give the degree of the series : "); get(deg);
    put("Give the number of monomials : "); get(nbr);
    put("Give the highest power of each variable : "); get(pwr);
    put("Give the frequency of the tests : "); get(frq);
    Test(dim,deg,nbr,pwr,frq);
  end Main;

begin
  Main;
end ts_perfcirc;
