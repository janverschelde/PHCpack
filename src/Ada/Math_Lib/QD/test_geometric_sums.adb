with text_io;                           use text_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Octo_Double_Numbers_io;            use Octo_Double_Numbers_io;
with Hexa_Double_Numbers_io;            use Hexa_Double_Numbers_io;
with Geometric_Sums;                    use Geometric_Sums;

package body Test_Geometric_Sums is

  procedure Test_Double_Sum
              ( dim : in integer32; ratio : in double_float ) is

    gsum : constant double_float := Double_Sum(dim,ratio);
    exact : constant double_float
          := (ratio**natural(dim+1) - 1.0)/(ratio - 1.0);
    err : double_float;

  begin
    put("the sum : "); put(gsum); new_line;
    err := abs(exact - gsum);
    put("  error : "); put(err,2); new_line;
  end Test_Double_Sum;

  procedure Test_Double_Double_Sum
              ( dim : in integer32; ratio : in double_double ) is

    gsum : constant double_double := Double_Double_Sum(dim,ratio);
    exact : constant double_double
          := (ratio**natural(dim+1) - 1.0)/(ratio - 1.0);
    err : double_double;

  begin
    put("the sum : "); put(gsum); new_line;
    err := abs(exact - gsum);
    put("  error : "); put(err,2); new_line;
  end Test_Double_Double_Sum;

  procedure Test_Quad_Double_Sum
              ( dim : in integer32; ratio : in quad_double ) is

    gsum : constant quad_double := Quad_Double_Sum(dim,ratio);
    exact : constant quad_double
          := (ratio**natural(dim+1) - 1.0)/(ratio - 1.0);
    err : quad_double;

  begin
    put("the sum : "); put(gsum); new_line;
    err := abs(exact - gsum);
    put("  error : "); put(err,2); new_line;
  end Test_Quad_Double_Sum;

  procedure Test_Octo_Double_Sum
              ( dim : in integer32; ratio : in octo_double ) is

    gsum : constant octo_double := Octo_Double_Sum(dim,ratio);
    exact : constant octo_double
          := (ratio**natural(dim+1) - 1.0)/(ratio - 1.0);
    err : octo_double;

  begin
    put("the sum : "); put(gsum); new_line;
    err := abs(exact - gsum);
    put("  error : "); put(err,2); new_line;
  end Test_Octo_Double_Sum;

  procedure Test_Hexa_Double_Sum
              ( dim : in integer32; ratio : in hexa_double ) is

    gsum : constant hexa_double := Hexa_Double_Sum(dim,ratio);
    exact : constant hexa_double
          := (ratio**natural(dim+1) - 1.0)/(ratio - 1.0);
    err : hexa_double;

  begin
    put("the sum : "); put(gsum); new_line;
    err := abs(exact - gsum);
    put("  error : "); put(err,2); new_line;
  end Test_Hexa_Double_Sum;

  procedure Test is

    dim : constant integer32 := 1_000_000;
    ratio : constant double_float := 1.0 - 1.0E-5;
    ddratio : double_double;
    qdratio : quad_double;
    odratio : octo_double;
    hdratio : hexa_double;

  begin
    put_line("Computing a geometric sum of size" & dim'Image);
    put("ratio :"); put(ratio); new_line;
    put_line("in double arithmetic ...");
    Test_Double_Sum(dim,ratio);
    ddratio := create(ratio);
    put_line("in double double arithmetic ...");
    Test_Double_Double_Sum(dim,ddratio);
    qdratio := create(ratio);
    put_line("in quad double arithmetic ...");
    Test_Quad_Double_Sum(dim,qdratio);
    odratio := create(ratio);
    put_line("in octo double arithmetic ...");
    Test_Octo_Double_Sum(dim,odratio);
    hdratio := create(ratio);
    put_line("in hexa double arithmetic ...");
    Test_Hexa_Double_Sum(dim,hdratio);
  end Test;

end Test_Geometric_Sums;
