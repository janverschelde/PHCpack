with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Multprec_Complex_Vectors;          use Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;          use Multprec_Complex_VecVecs;
with Multprec_Random_VecVecs;           use Multprec_Random_VecVecs;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;   use Multprec_Complex_Poly_Functions;
with Multprec_Polynomial_Interpolators; use Multprec_Polynomial_Interpolators;

procedure ts_mupolin is

-- DESCRIPTION :
--   Test on the multi-precision polynomial interpolators.

  procedure Multiply ( point : in out Vector; factor : in Floating_Number ) is

  -- DESCRIPTION :
  --   Multiplies the point with a factor.

  begin
    for i in point'range loop
      Mul(point(i),factor);
    end loop;
  end Multiply;

  procedure Multiply ( points : in out VecVec; factor : in Floating_Number ) is

  -- DESCRIPTION :
  --   Multiplies every point with a factor.

  begin
    for i in points'range loop
      if points(i) /= null
       then Multiply(points(i).all,factor);
      end if;
    end loop;
  end Multiply;

  procedure Interpolation_Data
              ( n,m,sz : in natural32; p : in Poly;
                points : in out VecVec ) is

  -- DESCRIPTION :
  --   Generates or reads in the interpolation points.

  -- WARNING : m = #monomials = #points+1.

    ans : character;
    fac : Floating_Number;

  begin
    new_line;
    put_line("MENU for entering the interpolation points :");
    put_line("  0. Generate random choices for the points;");
    put_line("  1. Sample points from the polynomial with one coefficients;");
    put("Type your answer (0 or 1) : "); Ask_Alternative(ans,"01");
    case ans is
      when '0' => points := Random_VecVec(n,m-1,sz);
                  put("Do you want to scale the data with a factor? (y/n) ");
                  Ask_Yes_or_No(ans);
                  if ans = 'y'
                   then put("Give the multiplier : "); get(fac);
                        Multiply(points,fac);
                  end if;
      when '1' => points := Sample(p,m-1,sz);
      when others => null;
    end case;
  end Interpolation_Data;

  procedure Interpolate ( file : in file_type; d,n,m,sz : in natural32;
                          p : in Poly; points : VecVec ) is

  -- DESCRIPTION :
  --   Constructs an interpolating polynomial through randomly generated
  --   points of the given size.

  -- ON ENTRY :
  --   file      to write the results on;
  --   d         degree of the interpolating polynomial;
  --   n         number of variables;
  --   m         number of monomials;
  --   sz        size of the numbers;
  --   p         polynomial of degree d in n variables, coefficients all one;
  --   points    interpolation points;

    ip : Poly;
    eva : Complex_Number;
    rcond,dist : Floating_Number;
    timer : Timing_Widget;

  begin
    tstart(timer);
    Interpolate1(p,points,ip,rcond);
    tstop(timer);
    put_line("The interpolating polynomial p :"); put_line(ip);
    put_line(file,"The interpolating polynomial p :");
    put_line(file,ip);
    put_line("Evaluating the interpolating polynomial at the points :");
    put_line(file,"Evaluating the interpolating polynomial at the points :");
    for i in points'range loop
      eva := Eval(ip,points(i).all);
      put("p("); put(i,1); put(") : "); put(eva,3); new_line;
      put(file,"p("); put(file,i,1); put(file,") : ");
      put(file,eva,3); new_line(file);
      Clear(eva);
    end loop;
    put("inverse of condition number : "); put(rcond,3); new_line;
    put(file,"inverse of condition number : "); put(file,rcond,3);
    new_line(file);
    dist := Distance(p,ip);
    put("distance between p and ip : "); put(dist,3); new_line;
    put(file,"distance between p and ip : "); put(file,dist,3);
    new_line(file);
    new_line(file);
    print_times(file,timer,"Constructing the Interpolator");
  end Interpolate;

  function Choose_Coefficients return natural32 is

    cff : natural32 := 0;

  begin
    new_line;
    put_line("MENU for choice of coefficients :");
    put_line("  0. random complex coefficients of modulus one;");
    put_line("  1. all coefficients are equal to one;");
    put_line("  2. random real coefficients will be generated.");
    put("Give a natural number (0,1, or 2) : "); get(cff);
    return cff;
  end Choose_Coefficients;

  procedure Experiment ( file : in file_type; n : in natural32 ) is

    d1,d2,m,deci,size,cff : natural32 := 0;

  begin
    new_line;
    put_line("Interpolating coefficient=one polynomial for various degrees.");
    new_line;
    put("Give the lower degree in range : "); get(d1);
    put("Give the upper degree in range : "); get(d2);
    put("Give the number of decimal places : "); get(deci);
    size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
    cff := Choose_Coefficients;
    new_line;
    put_line("See the output file for results...");
    new_line;
    for i in d1..d2 loop
      put(file,"d "); put(file,i,1);
      m := Number_of_Terms(i,n);
      put(file," m  "); put(file,m,3);
      declare
        p : Poly := Create(i,n,cff);
        points : VecVec(1..integer32(m)-1) := Sample(p,m-1,size);
        ip : Poly;
        tmp : Floating_Number;
        rcond,dist : double_float;
      begin
        Interpolate1(p,points,ip,tmp);
        rcond := Round(tmp); Clear(tmp);
        put(file," rco :"); put(file,rcond,3);
        put(file," |LOG10| : "); put(file,AbsVal(LOG10(rcond)),2,3,0);
        tmp := Distance(p,ip);
        dist := Round(tmp); Clear(tmp);
        put(file," dst :"); put(file,dist,3);
        put(file," |LOG10| : "); put(file,AbsVal(LOG10(dist)),2,3,0);
        new_line(file);
        Clear(p); Clear(ip);
      end;
    end loop;
  end Experiment;

  procedure Main is

    file : file_type;
    d,n,m,sz,cff : natural32 := 0;
    p : Poly;
    ans : character;

  begin
    new_line;
    put_line("Interpolating polynomials with multi-precision arithmetic");
    new_line;
    put_line("Reading the name of the output file.");
    Read_Name_and_Create_File(file);
    new_line;
    put("Give the number of variables : "); get(n);
    new_line;
    put_line("MENU for testing interpolation methods :");
    put_line("  1. interactive testing;");
    put_line("  2. perform experiment on conditioning.");
    put("Type 1 or 2 to select : "); Ask_Alternative(ans,"12");
    if ans = '2' then
      Experiment(file,n);
    else
      put("Give the degree of the polynomial : "); get(d);
      cff := Choose_Coefficients;
      m := Number_of_Terms(d,n);
      put("The number of terms : "); put(m,1); new_line;
      put(file,"#vars : "); put(file,n,1);
      put(file,"  degree : "); put(file,d,1);
      put(file,"  #terms : "); put(file,m,1); new_line(file);
      p := Create(d,n,cff); 
      put_line("The monomial structure of the interpolating polynomial :");
      put(p); new_line;
      put_line(file,
               "The monomial structure of the interpolating polynomial :");
      put(file,p); new_line(file);
      put("Give the size of the numbers : "); get(sz);
      declare
        points : VecVec(1..integer32(m)-1);
      begin
        Interpolation_Data(n,m,sz,p,points);
        Interpolate(file,d,n,m,sz,p,points);
      end;
    end if;
    Close(file);
  end Main;

begin
  Main;
end ts_mupolin;
