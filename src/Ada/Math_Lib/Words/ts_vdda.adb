with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;
with Double_Double_Basics;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Double_Double_Vectors;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with DoblDobl_Random_Vectors;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Numbers_io;             use Quad_Double_Numbers_io;

procedure ts_vdda is

-- DESCRIPTION :
--   Test on the vectorized double double arithmetic.

  procedure Split ( v : in DoblDobl_Complex_Vectors.Vector;
                    vre0 : out Standard_Floating_Vectors.Vector;
                    vre1 : out Standard_Floating_Vectors.Vector;
                    vre2 : out Standard_Floating_Vectors.Vector;
                    vre3 : out Standard_Floating_Vectors.Vector;
                    vim0 : out Standard_Floating_Vectors.Vector;
                    vim1 : out Standard_Floating_Vectors.Vector;
                    vim2 : out Standard_Floating_Vectors.Vector;
                    vim3 : out Standard_Floating_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Splits the doubles of the real and imaginary parts of the
  --   complex double doubles in the vector v.

     nbr : double_double;
     flt : double_float;

  begin
    for i in v'range loop
      nbr := DoblDobl_Complex_Numbers.REAL_PART(v(i));
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,vre0(i),vre1(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,vre2(i),vre3(i));
      nbr := DoblDobl_Complex_Numbers.IMAG_PART(v(i));
      flt := hi_part(nbr);
      Double_Double_Basics.Split(flt,vim0(i),vim1(i));
      flt := lo_part(nbr);
      Double_Double_Basics.Split(flt,vim2(i),vim3(i));
    end loop;
  end Split;

  function to_double_double
             ( x0,x1,x2,x3 : in double_float )
             return double_double is

  -- DESCRIPTION :
  --   Uses a quad double to make a double double from
  --   the four given doubles.

    res : double_double;
    precision : constant natural32 := 32;
    s : string(1..integer(precision)+10);
    ends,pos : integer;
    x : constant quad_double := create(x0,x1,x2,x3);
    fail : boolean;

  begin
    to_string(x,precision,0,false,false,true,' ',s,ends);
    pos := s'first;
    read(s,pos,res,fail);
    return res;
  end to_double_double;

  procedure Test_Sum ( v : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Sums the real and imaginary parts twice,
  --   once with double double arithmetic, and
  --   once on the splitted high and low parts.

    dim : constant integer32 := v'last;
    re_sum,im_sum : double_double := create(0.0);
    vre0 : Standard_Floating_Vectors.Vector(v'range);
    vre1 : Standard_Floating_Vectors.Vector(v'range);
    vre2 : Standard_Floating_Vectors.Vector(v'range);
    vre3 : Standard_Floating_Vectors.Vector(v'range);
    vim0 : Standard_Floating_Vectors.Vector(v'range);
    vim1 : Standard_Floating_Vectors.Vector(v'range);
    vim2 : Standard_Floating_Vectors.Vector(v'range);
    vim3 : Standard_Floating_Vectors.Vector(v'range);
    re_sum0,re_sum1,re_sum2,re_sum3 : double_float := 0.0;
    im_sum0,im_sum1,im_sum2,im_sum3 : double_float := 0.0;
    qdre_sum,qdim_sum,err : quad_double;
    ddre_sum,ddim_sum : double_double;

  begin
    for i in v'range loop
      re_sum := re_sum + DoblDobl_Complex_Numbers.REAL_PART(v(i));
      im_sum := im_sum + DoblDobl_Complex_Numbers.IMAG_PART(v(i));
    end loop;
    new_line;
    put("Performed "); put(40*dim); put_line(" operations.");
    Split(v,vre0,vre1,vre2,vre3,vim0,vim1,vim2,vim3);
    for i in v'range loop
      re_sum0 := re_sum0 + vre0(i);
      re_sum1 := re_sum1 + vre1(i);
      re_sum2 := re_sum2 + vre2(i);
      re_sum3 := re_sum3 + vre3(i);
      im_sum0 := im_sum0 + vim0(i);
      im_sum1 := im_sum1 + vim1(i);
      im_sum2 := im_sum2 + vim2(i);
      im_sum3 := im_sum3 + vim3(i);
    end loop;
    put("after split: "); put(8*dim); put_line(" additions.");
    qdre_sum := create(re_sum0,re_sum1,re_sum2,re_sum3);
    ddre_sum := to_double_double(re_sum0,re_sum1,re_sum2,re_sum3);
    qdim_sum := create(im_sum0,im_sum1,im_sum2,im_sum3);
    ddim_sum := to_double_double(im_sum0,im_sum1,im_sum2,im_sum3);
    put("re(sum) : "); put(re_sum); new_line;
    put("re(sum) : "); put(qdre_sum,32); new_line;
    put("re(sum) : "); put(ddre_sum); new_line;
    err := abs(qdre_sum - re_sum);
    put("error : "); put(err,3); new_line;
    put("im(sum) : "); put(im_sum); new_line;
    put("im(sum) : "); put(qdim_sum,32); new_line;
    put("im(sum) : "); put(ddim_sum); new_line;
    err := abs(qdim_sum - im_sum);
    put("error : "); put(err,3); new_line;
  end Test_Sum;

  procedure Test_Product ( v : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Separates the real and imaginary parts of v into x and y
  --   and makes the inner products of x and y twice,
  --   once with double double arithmetic, and
  --   once on the splitted high and low parts.

    x,y : Double_Double_Vectors.Vector(v'range);
    prd : double_double := create(0.0);
    x0 : Standard_Floating_Vectors.Vector(v'range);
    x1 : Standard_Floating_Vectors.Vector(v'range);
    x2 : Standard_Floating_Vectors.Vector(v'range);
    x3 : Standard_Floating_Vectors.Vector(v'range);
    y0 : Standard_Floating_Vectors.Vector(v'range);
    y1 : Standard_Floating_Vectors.Vector(v'range);
    y2 : Standard_Floating_Vectors.Vector(v'range);
    y3 : Standard_Floating_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      x(i) := DoblDobl_Complex_Numbers.REAL_PART(v(i));
      y(i) := DoblDobl_Complex_Numbers.IMAG_PART(v(i));
      prd := prd + x(i)*y(i);
    end loop;
    Split(v,x0,x1,x2,x3,y0,y1,y2,y3);
  end Test_Product;

  procedure Test ( dim : in integer32 ) is

  -- DESCRIPTION :
  --   Generates a random vector of double double numbers
  --   of dimension dim.

    ddv : constant DoblDobl_Complex_Vectors.Vector(1..dim)
        := DoblDobl_Random_Vectors.Random_Vector(1,dim);

  begin
    if dim < 20 then
      put_line("A random complex vector :");
      put(ddv);
    else
      put("Testing on a random vector of dimension "); put(dim,1);
      put_line(" ...");
    end if;
    Test_Sum(ddv);
    Test_Product(ddv);
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the dimension and runs the tests.

    dim : integer32 := 0;

  begin
    put("Give the dimension : "); get(dim);
    Test(dim);
  end Main;

begin
  Main;
end ts_vdda;
