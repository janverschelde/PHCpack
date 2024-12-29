with Standard_Random_Numbers;

package body Balanced_Quarter_Doubles is

  function Thirteen_Random_Bits return integer64 is

    res : integer64 := 1;
    rnd : integer64;

  begin
    for i in 1..12 loop
      rnd := Standard_Random_Numbers.Random(0,1);
      res := 2*res + rnd;
    end loop;
    return res;
  end Thirteen_Random_Bits;

  function Random_Quarter ( e : in integer32 ) return double_float is

    frc : constant integer64 := Thirteen_Random_Bits;
    mrs : constant double_float := double_float(frc);
    res : constant double_float := double_float'compose(mrs,e);

  begin
    return res;
  end Random_Quarter;

  procedure Random ( x0,x1,x2,x3 : out double_float ) is
  begin
    x0 := Random_Quarter(-1);
    x1 := Random_Quarter(-13);
    x2 := Random_Quarter(-25);
    x3 := Random_Quarter(-37);
  end Random;

  procedure Random ( dim : in integer32;
                     x0,x1,x2,x3 : out Standard_Floating_Vectors.Vector ) is
  begin
    for i in 1..dim loop
      Random(x0(i),x1(i),x2(i),x3(i));
    end loop;
  end Random;

end Balanced_Quarter_Doubles;
