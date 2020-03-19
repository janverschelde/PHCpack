with Machines;                           use Machines;
with Standard_Mathematical_Functions;    use Standard_Mathematical_Functions;

package body Standard_Random_Numbers is

  a : constant integer32 := 16807;
  m : constant integer32 := 2147483647;
  pid : constant integer := Process_Id;
  pid32 : constant integer32 := integer32(pid);
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

  function Random ( probability : double_float := 0.5 ) return boolean is

    f : constant double_float := (Random + 1.0)/2.0;

  begin
    return (f <= probability);
  end Random;

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

  procedure Random_Integer_Number
              ( seed : in out integer32; lower,upper : in integer32;
                i : out integer32 ) is

    f : double_float;

  begin
    Random_Double_Float(seed,f); -- f is in [-1,+1]
    f := 0.5*(double_float(upper-lower)*f + double_float(lower+upper));
                                         -- f in [lower,upper]
    i := integer32(f); -- rounding the result to integer number
  end Random_Integer_Number;

  procedure Random_Integer_Number
              ( seed : in out integer32; lower,upper : in integer64;
                i : out integer64 ) is

    f : double_float;

  begin
    Random_Double_Float(seed,f); -- f is in [-1,+1]
    f := 0.5*(double_float(upper-lower)*f + double_float(lower+upper));
                                         -- f in [lower,upper]
    i := integer64(f); -- rounding the result to integer number
  end Random_Integer_Number;

  function Random return double_float is

    x : double_float;
    p : constant integer64 := integer64(a)*integer64(seed);
    m64 : constant integer64 := integer64(m);
    pmodm : constant integer64 := p mod m64;
     
  begin
    seed := integer32(pmodm);
    x := double_float(seed)/double_float(m);
    x := 2.0 * x - 1.0;
    return x;
  end Random;

  procedure Random_Double_Float 
              ( seed : in out integer32; f : out double_float ) is

    p : constant integer64 := integer64(a)*integer64(seed);
    m64 : constant integer64 := integer64(m);
    pmodm : constant integer64 := p mod m64;

  begin
    seed := integer32(pmodm); -- a*seed mod m;
    f := double_float(seed)/double_float(m);
    f := 2.0 * f - 1.0;
  end Random_Double_Float;
  
  function Random return Complex_Number is
  begin
    return Create(Random,Random);
  end Random;

  procedure Random_Complex_Number
              ( seed : in out integer32; c : out Complex_Number ) is

    cre,cim : double_float;

  begin
    Random_Double_Float(seed,cre);
    Random_Double_Float(seed,cim);
    c := Create(cre,cim);
  end Random_Complex_Number;

  function Random ( modulus : double_float ) return Complex_Number is

    arg : constant double_float := PI*Random;  -- in [-pi,+pi]

  begin
    return Create(modulus*COS(arg),modulus*SIN(arg));
  end Random;

  procedure Random_Complex_Number
              ( seed : in out integer32; modulus : in double_float;
                c : out Complex_Number ) is

    arg : double_float;

  begin
    Random_Double_Float(seed,arg);
    arg := PI*arg; -- in [-pi, +pi]
    c := Create(modulus*COS(arg),modulus*SIN(arg));
  end Random_Complex_Number;

  function Random1 return Complex_Number is

    arg : double_float; 

  begin
    arg := PI*Random;  -- in [-pi,+pi]
    return Create(COS(arg),SIN(arg));
  end Random1;

  procedure Random1_Complex_Number
              ( seed : in out integer32; c : out Complex_Number ) is

    arg : double_float; 

  begin
    Random_Double_Float(seed,arg);
    arg := PI*arg;  -- in [-pi,+pi]
    c := Create(COS(arg),SIN(arg));
  end Random1_Complex_Number;

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
