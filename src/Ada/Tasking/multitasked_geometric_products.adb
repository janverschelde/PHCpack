with Ada.Text_IO;
with Generic_Task_Array;
with Double_Double_Numbers_io;           use Double_Double_Numbers_io;
with Quad_Double_Numbers_io;           use Quad_Double_Numbers_io;

package body Multitasked_Geometric_Products is

  function Inner_Product ( dim : integer32; rtx,rty : double_float ) 
                         return double_float is

    result : double_float := 0.0;
    x : double_float := 1.0;
    y : double_float := 1.0;

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( dim : integer32; rtx,rty : double_double ) 
                         return double_double is

    result : double_double := create(0.0);
    x : double_double := create(1.0);
    y : double_double := create(1.0);

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( dim : integer32; rtx,rty : quad_double ) 
                         return quad_double is

    result : quad_double := create(0.0);
    x : quad_double := create(1.0);
    y : quad_double := create(1.0);

  begin
    for i in 0..(dim-1) loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( idx_start,idx_end : integer32;
                           rtx,rty : double_float ) 
                         return double_float is

    result : double_float := 0.0;
    x : double_float := rtx**natural(idx_start);
    y : double_float := rty**natural(idx_start);

  begin
    for i in idx_start..idx_end loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( idx_start,idx_end : integer32;
                           rtx,rty : double_double ) 
                         return double_double is

    result : double_double := create(0.0);
    x : double_double := rtx**natural(idx_start);
    y : double_double := rty**natural(idx_start);

  begin
    for i in idx_start..idx_end loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  function Inner_Product ( idx_start,idx_end : integer32;
                           rtx,rty : quad_double ) 
                         return quad_double is

    result : quad_double := create(0.0);
    x : quad_double := rtx**natural(idx_start);
    y : quad_double := rty**natural(idx_start);

  begin
    for i in idx_start..idx_end loop
      result := result + x*y;
      x := x*rtx;
      y := y*rty;
    end loop;
    return result;
  end Inner_Product;

  procedure Double_Sequential_Run is

    d : constant integer32 := 1_000_000_000;
    r : constant double_float := 1.0 + 1.0E-10;
    s : constant double_float := 1.0 - 1.0E-10;
    p : double_float;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put_Line("The inner product :" & p'Image);
  end Double_Sequential_Run;

  procedure Double_Double_Sequential_Run is

    d : constant integer32 := 1_000_000_000;
    one : constant double_double := create(1.0);
    r : constant double_double
      := one + Double_Double_Numbers.Create(1.0E-10);
    s : constant double_double
      := one - Double_Double_Numbers.Create(1.0E-10);
    p : double_double;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put("The inner product : "); put(p);
    Ada.Text_IO.New_Line;
  end Double_Double_Sequential_Run;

  procedure Quad_Double_Sequential_Run is

    d : constant integer32 := 1_000_000_000;
    one : constant quad_double := create(1.0);
    r : constant quad_double
      := one + Quad_Double_Numbers.Create(1.0E-10);
    s : constant quad_double
      := one - Quad_Double_Numbers.Create(1.0E-10);
    p : quad_double;
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    p := Inner_Product(d,r,s);
    Ada.Text_IO.Put("The inner product : "); put(p);
    Ada.Text_IO.New_Line;
  end Quad_Double_Sequential_Run;

  procedure Double_Parallel_Run ( p : in integer32 ) is

    d : constant integer32 := 1_000_000_000;
    r : constant double_float := 1.0 + 1.0E-10;
    s : constant double_float := 1.0 - 1.0E-10;
    q : double_float := 0.0;
    a : array(1..p) of double_float;
    m : constant integer32 := d/p;

    procedure Do_Inner_Product ( idn : in integer ) is

    -- DESCRIPTION :
    --   Runs an inner product for start and end range,
    --   defined by the index number idn.

      id : constant integer32 := integer32(idn);
      idx_start : constant integer32 := m*(id-1);
      idx_end : constant integer32 := m*id - 1;

    begin
      Ada.Text_IO.Put_Line("Task" & idn'Image & " is computing ...");
      a(id) := Inner_Product(idx_start,idx_end,r,s);
    end Do_Inner_Product;

    procedure Do_Inner_Products is
      new Generic_Task_Array(Do_Inner_Product);
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Do_Inner_Products(integer(p));
    for i in a'range loop
      q := q + a(i);
    end loop;
    Ada.Text_IO.Put_Line("The inner product :" & q'Image);
  end Double_Parallel_Run;

  procedure Double_Double_Parallel_Run ( p : in integer32 ) is

    d : constant integer32 := 1_000_000_000;
    one : constant double_double := create(1.0);
    r : constant double_double
      := one + Double_Double_Numbers.create(1.0E-10);
    s : constant double_double
      := one - Double_Double_Numbers.Create(1.0E-10);
    q : double_double := create(0.0);
    a : array(1..p) of double_double;
    m : constant integer32 := d/p;

    procedure Do_Inner_Product ( idn : in integer ) is

    -- DESCRIPTION :
    --   Runs an inner product for start and end range,
    --   defined by the index number idn.

      id : constant integer32 := integer32(idn);
      idx_start : constant integer32 := m*(id-1);
      idx_end : constant integer32 := m*id - 1;

    begin
      Ada.Text_IO.Put_Line("Task" & idn'Image & " is computing ...");
      a(id) := Inner_Product(idx_start,idx_end,r,s);
    end Do_Inner_Product;

    procedure Do_Inner_Products is
      new Generic_Task_Array(Do_Inner_Product);
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Do_Inner_Products(integer(p));
    for i in a'range loop
      q := q + a(i);
    end loop;
    Ada.Text_IO.Put_Line("The inner product : "); put(q);
    Ada.Text_IO.New_Line;
  end Double_Double_Parallel_Run;

  procedure Quad_Double_Parallel_Run ( p : in integer32 ) is

    d : constant integer32 := 1_000_000_000;
    one : constant quad_double := create(1.0);
    r : constant quad_double
      := one + Quad_Double_Numbers.create(1.0E-10);
    s : constant quad_double
      := one - Quad_Double_Numbers.Create(1.0E-10);
    q : quad_double := create(0.0);
    a : array(1..p) of quad_double;
    m : constant integer32 := d/p;

    procedure Do_Inner_Product ( idn : in integer ) is

    -- DESCRIPTION :
    --   Runs an inner product for start and end range,
    --   defined by the index number idn.

      id : constant integer32 := integer32(idn);
      idx_start : constant integer32 := m*(id-1);
      idx_end : constant integer32 := m*id - 1;

    begin
      Ada.Text_IO.Put_Line("Task" & idn'Image & " is computing ...");
      a(id) := Inner_Product(idx_start,idx_end,r,s);
    end Do_Inner_Product;

    procedure Do_Inner_Products is
      new Generic_Task_Array(Do_Inner_Product);
 
  begin
    Ada.Text_IO.Put_Line("Computing an inner product of size" & d'Image);
    Do_Inner_Products(integer(p));
    for i in a'range loop
      q := q + a(i);
    end loop;
    Ada.Text_IO.Put_Line("The inner product : "); put(q);
    Ada.Text_IO.New_Line;
  end Quad_Double_Parallel_Run;

end Multitasked_Geometric_Products;
