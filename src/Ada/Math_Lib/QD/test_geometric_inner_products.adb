with Ada.Text_IO;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with Octo_Double_Numbers;               use Octo_Double_Numbers;
with Octo_Double_Numbers_io;            use Octo_Double_Numbers_io;
with Hexa_Double_Numbers;               use Hexa_Double_Numbers;
with Hexa_Double_Numbers_io;            use Hexa_Double_Numbers_io;
with Geometric_Inner_Products;          use Geometric_Inner_Products;

package body Test_Geometric_Inner_Products is

  procedure Test_Double_Float is

    d : constant integer32 := 1_000_000_000;
    r : constant double_float := 1.0 + 1.0E-10;
    s : constant double_float := 1.0 - 1.0E-10;
    p : double_float;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Ada.Text_IO.Put_Line("with double arithmetic ...");
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put_Line("The inner product :" & p'Image);
  end Test_Double_Float;

  procedure Test_Double_Double is

    d : constant integer32 := 1_000_000_000;
    r : constant double_double
      := Double_Double_Numbers.create(1.0)
       + Double_Double_Numbers.create(1.0E-10);
    s : constant double_double
      := Double_Double_Numbers.create(1.0)
       - Double_Double_Numbers.create(1.0E-10);
    p : double_double;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Ada.Text_IO.Put_Line("with double double arithmetic ...");
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put("The inner product : "); put(p); 
    Ada.Text_IO.New_Line;
  end Test_Double_Double;

  procedure Test_Quad_Double is

    d : constant integer32 := 1_000_000_000;
    r : constant quad_double
      := Quad_Double_Numbers.create(1.0) + Quad_Double_Numbers.create(1.0E-10);
    s : constant quad_double
      := Quad_Double_Numbers.create(1.0) - Quad_Double_Numbers.create(1.0E-10);
    p : quad_double;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Ada.Text_IO.Put_Line("with quad double arithmetic ...");
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put("The inner product : "); put(p); 
    Ada.Text_IO.New_Line;
  end Test_Quad_Double;

  procedure Test_Octo_Double is

    d : constant integer32 := 1_000_000_000;
    r : constant octo_double
      := Octo_Double_Numbers.create(1.0) + Octo_Double_Numbers.create(1.0E-10);
    s : constant octo_double
      := Octo_Double_Numbers.create(1.0) - Octo_Double_Numbers.create(1.0E-10);
    p : octo_double;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Ada.Text_IO.Put_Line("with octo double arithmetic ...");
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put("The inner product : "); put(p); 
    Ada.Text_IO.New_Line;
  end Test_Octo_Double;

  procedure Test_Hexa_Double is

    d : constant integer32 := 1_000_000_000;
    r : constant hexa_double
      := Hexa_Double_Numbers.create(1.0) + Hexa_Double_Numbers.create(1.0E-10);
    s : constant hexa_double
      := Hexa_Double_Numbers.create(1.0) - Hexa_Double_Numbers.create(1.0E-10);
    p : hexa_double;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Ada.Text_IO.Put_Line("with hexa double arithmetic ...");
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put("The inner product : "); put(p); 
    Ada.Text_IO.New_Line;
  end Test_hexa_Double;

  procedure Test is
  begin
    Test_Double_Float;
    Test_Double_Double;
    Test_Quad_Double;
    Test_Octo_Double;
    Test_Hexa_Double;
  end Test;

end Test_Geometric_Inner_Products;
