with Machines;                           use Machines;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Standard_Random_Numbers is

  a : constant integer32 := 16807;
  m : constant integer32 := 2147483647;
  pid : constant integer := Process_Id;
  pid32 : constant integer32 := integer32(pid);
 -- pid : constant integer := 186;  -- way to fix the seed
  initial_seed : integer32 := pid32;
  seed : integer32 := initial_seed;

  procedure Set_Seed ( n : natural32 ) is
  begin
    initial_seed := integer32(n);
    seed := initial_seed;
  end Set_Seed;

  function Get_Seed return integer32 is
  begin
    return initial_seed;
  end Get_Seed;

  function Random ( lower,upper : integer32 ) return integer32 is

    f : double_float;

  begin
    f := Random; -- f in [-1,1]
    f := 0.5*(double_float(upper-lower)*f + double_float(lower+upper));
                                         -- f in [lower,upper]
    return integer32(f); -- rounding the result to integer number
  end Random;

  function Random ( lower,upper : integer64 ) return integer64 is

    f : double_float;

  begin
    f := Random; -- f in [-1,1]
    f := 0.5*(double_float(upper-lower)*f + double_float(lower+upper));
                                         -- f in [lower,upper]
    return integer64(f); -- rounding the result to integer number
  end Random;

  function Random return double_float is

    x : double_float;

  begin
    seed := a*seed mod m;
    x := double_float(seed)/double_float(m);
    x := 2.0 * x - 1.0;
    return x;
  end Random;
  
  function Random return Complex_Number is
  begin
    return Create(Random,Random);
  end Random;

  function Random ( modulus : double_float ) return Complex_Number is

    arg : double_float;

  begin
    arg := PI*Random;  -- in [-pi,+pi]
    return Create(modulus*COS(arg),modulus*SIN(arg));
  end Random;

  function Random1 return Complex_Number is

    arg : double_float; 

  begin
    arg := PI*Random;  -- in [-pi,+pi]
    return Create(COS(arg),SIN(arg));
  end Random1;

  function Random_Magnitude ( m : natural32 ) return double_float is

    r : constant double_float := Random;
    low : constant integer32 := -integer32(m);
    upp : constant integer32 := integer32(m);
    k : constant integer32 := Random(low,upp);
    f : double_float := abs(r);

  begin
    if k < 0 then
      for i in 1..((-k)-1) loop
        f := f/10.0;
      end loop;
    else
      for i in 1..k loop
        f := f*10.0;
      end loop;
    end if;
    return f;
  end Random_Magnitude;

  function Random_Magnitude ( m : natural32 ) return Complex_Number is

    r : constant double_float := Random_Magnitude(m);
    c : constant Complex_Number := Random(r);

  begin
    return c;
  end Random_Magnitude;

end Standard_Random_Numbers;
