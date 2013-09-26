with Communications_with_User;           use Communications_with_User;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
--with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Numbers_io;                         use Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Interpolating_Homotopies;           use Interpolating_Homotopies;

procedure Interpolating_Homotopies_Driver
             ( file : in file_type; p : in Poly_Sys; z : in Partition;
               b : in out natural32; q : out Poly_Sys; 
               qsols : in out Solution_List ) is

  n : constant natural32 := natural32(p'last);
  interpols : Solution_List;
  ib,scalind : natural32;
  ans : character;

  function Random_Interpolating ( n,m : natural32 ) return Solution_List is

  -- DESCRIPTION :
  --   A list of m random n-dimensional vectors will be returned.
  --   The complex numbers will all have modulus one.

    res,res_last : Solution_List;
    s : Solution(integer32(n));

  begin
    s.t := Create(0.0);
    s.m := 1;
    s.err := 0.0; s.rco := 1.0; s.res := 0.0;
    for i in 1..m loop
      for j in 1..integer32(n) loop
        s.v(j) := random1;
      end loop;
      Append(res,res_last,s);
    end loop;
    return res;
  end Random_Interpolating;

  procedure Random_Linear_Scaler
              ( n : in natural32; p : in Poly;
                v : out vector; l : out integer32 ) is

  -- DESCRIPTION :
  --   Returns a random vector of dimension n+1, with range 0..n.
  --   There will be a nonzero entry only for those unknowns that occur in p.

    res : Vector(0..integer32(n));
    last : integer32 := 0;

  begin
    for i in res'range loop
      if Degree(p,i) > 0 then
        res(i) := random1;
        last := last + 1;
      else
        res(i) := Create(0.0);
      end if;
    end loop;
    v := res; l := last;
  end Random_Linear_Scaler;

  function Scale ( sc,v : Vector; last : integer32 ) return Complex_Number is

  -- DESCRIPTION :
  --   Returns the last component of the vector v, that is v(last),
  --   such that sc(0) + sum sc(i)*v(i), i in v'range, holds.

    res : Complex_Number := sc(0);

  begin
    for i in v'first..last-1 loop
      res := res + sc(i)*v(i);
    end loop;
    res := -res/sc(last);
    return res;
  end Scale;

  function Random_Interpolating
             ( n,m : natural32; scaler : Vector; scallast : integer32 ) 
             return Solution_List is

  -- DESCRIPTION :
  --   A list of m random n-dimensional vectors will be returned.
  --   The complex numbers will all have modulus one, except the last one,
  --   indicated by scallast, that has been chosen to satisfy the scaler 
  --   equation, defined by sum of scaler(i)*x_i = 0, with x_0 = 1.

    res,res_last : Solution_List;
    s : Solution(integer32(n));

  begin
    s.t := Create(0.0);
    s.m := 1;
    for i in 1..integer32(m) loop
      for j in 1..(integer32(n)-1) loop
        s.v(j) := random1;
      end loop;
      s.v(integer32(n)) := Scale(scaler,s.v,scallast);
      Append(res,res_last,s);
    end loop;
    return res;
  end Random_Interpolating;

  function Create ( v : Vector ) return Poly is

    res : Poly := Null_Poly;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(v'first+1..v'last => 0);
    for i in v'range loop
      t.cf := v(i);
      if i > v'first
       then t.dg(i) := 1;
      end if;
      Add(res,t);
      if i > v'first
       then t.dg(i) := 0;
      end if;
    end loop;
    Clear(t);
    return res;
  end Create;

  function Interpolating_by_User ( n,m : natural32 ) return Solution_List is

  -- DESCRIPTION :
  --   A list of m n-dimensional vectors will be read from standard input.

    res,res_last : Solution_List;
    s : Solution(integer32(n));
   -- f1,f2 : double_float;

  begin
    put("Reading "); put(m,1); put(" "); put(n,1);
    put_line("-dimensional complex vectors.");
    for i in 1..m loop
      s.t := Create(0.0);
      s.m := 1;
      s.err := 0.0; s.rco := 1.0; s.res := 0.0;
      put("Give the components of vector "); put(i,1); 
      put_line(" :");
      for j in 1..integer32(n) loop
       -- Read_Double_Float(f1);
       -- Read_Double_Float(f2);
       -- s.v(j) := Create(f1,f2);
        get(s.v(j));
      end loop;
      Append(res,res_last,s);
    end loop;
    return res;
  end Interpolating_by_User;

  procedure Driver_for_Interpolation is

  -- DESCRIPTION : interpolation without a scaling equation

    dp,ip : Poly_Sys(p'range);

  begin
    dp := Dense_Representation(p,z);
    ip := Independent_Representation(dp);
    ib := Independent_Roots(ip);
    if ib > b
     then ib := b;
    end if;
    put("The number of independent roots : "); put(ib,1); new_line;
    put("Do you want to give interpolation vectors by yourself ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then interpols := Interpolating_by_User(n,ib);
     else put_line("Random interpolating vectors will be generated.");
          interpols := Random_Interpolating(n,ib);
    end if;
    q := Interpolate(ip,ib,interpols);
    qsols := interpols; b := ib;
    Clear(dp); Clear(ip);
  end Driver_for_Interpolation;

  procedure Driver_for_Scaled_Interpolation is

  -- DESCRIPTION : interpolation with a scaling equation, p(scalind).

    dp,ip,pp,qq : Poly_Sys(p'range);
    scalvec : Vector(0..integer32(n));
    scalveclast : integer32;

  begin
    Random_Linear_Scaler(n,p(integer32(scalind)),scalvec,scalveclast);
    for i in p'range loop
      if i = integer32(scalind)
       then pp(i) := Null_Poly;
       else pp(i) := p(i);
      end if;
    end loop;
    dp := Dense_Representation(pp,z);
    ip := Independent_Representation(dp);
    ib := Independent_Roots(ip,natural32(scalveclast));
    if ib > b
     then ib := b;
    end if;
    put("The number of independent roots : "); put(ib,1); new_line;
    put("Do you want to give interpolation vectors by yourself ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then interpols := Interpolating_by_User(n,ib);
     else put_line("Random interpolating vectors will be generated.");
          interpols := Random_Interpolating(n,ib,scalvec,scalveclast);
    end if;
    qq := Interpolate(ip,natural32(scalveclast),ib,interpols);
    qq(integer32(scalind)) := Create(scalvec);
    qsols := interpols; b := ib; q := qq;
    Clear(dp); Clear(ip);
  end Driver_for_Scaled_Interpolation;

begin
  new_line;
  put("Give the number of the linear scaling equation (0 if none) : ");
  Read_Natural(scalind);
  if scalind = 0
   then Driver_for_Interpolation;
   else Driver_for_Scaled_Interpolation;
  end if;
end Interpolating_Homotopies_Driver;
