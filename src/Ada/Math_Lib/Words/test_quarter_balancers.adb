with Ada.Text_IO;                        use Ada.Text_IO;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;
with Bits_of_Doubles;
with Test_Bits_of_Doubles;
with Random_Balanced_Quarters;
with Quarter_Balancers;

package body Test_Quarter_Balancers is

  procedure Test_Random_Double is

    x,x0,x1,x2,x3 : double_float;
    isbal,b01,b12,b23 : boolean;
    seed : natural32 := 0;

  begin
    new_line;
    put("Give seed for random numbers (0 is default) : ");
    get(seed);
    if seed /= 0
     then Standard_Random_Numbers.set_seed(seed);
    end if;
    x := abs(Standard_Random_Numbers.Random);
    new_line;
    put("x : "); put(x); new_line;
    put("xb : "); Bits_of_Doubles.write_52bits_expo(x); new_line;
    Bits_of_Doubles.Split(x,x0,x1,x2,x3);
    put("x0 : "); Bits_of_Doubles.write_52bits_expo(x0); new_line;
    put("x1 : "); Bits_of_Doubles.write_52bits_expo(x1); new_line;
    put("x2 : "); Bits_of_Doubles.write_52bits_expo(x2); new_line;
    put("x3 : "); Bits_of_Doubles.write_52bits_expo(x3); new_line;
    isbal := Random_Balanced_Quarters.Is_Balanced(0,x0,x1,x2,x3);
    if not isbal
     then put_line("unbalanced");
     else put_line("balanced");
    end if;
    put_line("Testing x0 and x1 ...");
    b01 := Quarter_Balancers.Is_Quarter_Balanced(x0,x1,1);
    put_line("Testing x1 and x2 ...");
    b12 := Quarter_Balancers.Is_Quarter_Balanced(x1,x2,1);
    put_line("Testing x2 and x3 ...");
    b23 := Quarter_Balancers.Is_Quarter_Balanced(x2,x3,1);
    if b01 and b12 and b23 then
      put_line("Quarters are balanced.");
    else
      put_line("Quarters are not balanced.");
      if not b01 then
        put_line("-> balancing x0 and x1 ...");
        Quarter_Balancers.Quarter_Balance(x0,x1,1);
        b01 := Quarter_Balancers.Is_Quarter_Balanced(x0,x1,1);
        if not b01
         then put_line("x0 and x1 are not balanced!  Bug?!");
        end if;
      end if;
      if not b12 then
        put_line("-> balancing x1 and x2 ...");
        Quarter_Balancers.Quarter_Balance(x1,x2,1);
        b12 := Quarter_Balancers.Is_Quarter_Balanced(x1,x2,1);
        if not b12
         then put_line("x1 and x2 are not balanced!  Bug?!");
        end if;
      end if;
      if not b23 then
        put_line("-> balancing x2 and x3 ...");
        Quarter_Balancers.Quarter_Balance(x2,x3,1);
        b23 := Quarter_Balancers.Is_Quarter_Balanced(x2,x3,1);
        if not b23
         then put_line("x2 and x3 are not balanced!  Bug?!");
        end if;
      end if;
    end if;
    Test_Bits_of_Doubles.Test_Bit_Split(x,x0,x1,x2,x3);
    seed := natural32(Standard_Random_Numbers.get_seed);
    new_line;
    put("seed used : "); put(seed,1); new_line;
  end Test_Random_Double;

  procedure Test_Quarter_Balance is

    x,x0,x1,x2,x3 : double_float;
    seed : natural32 := 0;

  begin
    new_line;
    put("Give seed for random numbers (0 is default) : ");
    get(seed);
    if seed /= 0
     then Standard_Random_Numbers.set_seed(seed);
    end if;
    x := abs(Standard_Random_Numbers.Random);
    new_line;
    put("x : "); put(x); new_line;
    Bits_of_Doubles.Split(x,x0,x1,x2,x3);
    Quarter_Balancers.Quarter_Balance(x0,x1,x2,x3,2);
    Test_Bits_of_Doubles.Test_Bit_Split(x,x0,x1,x2,x3);
    seed := natural32(Standard_Random_Numbers.get_seed);
    new_line;
    put("seed used : "); put(seed,1); new_line;
  end Test_Quarter_Balance;

  procedure Main is
  begin
    Test_Random_Double;
    Test_Quarter_Balance;
  end Main;

end Test_Quarter_Balancers;
