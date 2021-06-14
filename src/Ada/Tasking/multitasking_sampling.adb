with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Random_Numbers;             use Standard_Random_Numbers;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;         use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Standard_Complex_Polynomials;
with Standard_System_and_Solutions_io;    use Standard_System_and_Solutions_io;
with Standard_Homotopy;
with Planes_and_Polynomials;
with Continuation_Parameters;
with Witness_Sets,Monodromy_Partitions;
with Multitasking_Continuation;

package body Multitasking_Sampling is

-- NOTE : the samplers are actually wrappers around the path trackers
--   of the package Multitasking_Continuation.

  function Change_Slices
               ( p : in Poly_Sys;
                 h : in Standard_Complex_VecVecs.VecVec )
               return Poly_Sys is

    d : constant integer32 := h'last;
    res : Poly_Sys(p'range);

  begin
    for i in p'first..p'last-d loop
      Standard_Complex_Polynomials.Copy(p(i),res(i));
    end loop;
    for i in h'range loop
      res(p'last-d+i) := Planes_and_Polynomials.Hyperplane(h(i).all);
    end loop;
    return res;
  end Change_Slices;

  procedure Silent_Sampler
               ( s : in out Solution_List; n,d : in integer32;
                 h : in Standard_Complex_VecVecs.VecVec;
                 gamma : in Complex_Number;
                 q : in Poly_Sys; p : out Poly_Sys ) is
  begin
    p := Change_Slices(q,h);
    Standard_Homotopy.Create(p,q,2,gamma); 
    Multitasking_Continuation.Silent_Multitasking_Path_Tracker(s,n);
    Standard_Homotopy.Clear;
  end Silent_Sampler;

  procedure Reporting_Sampler
               ( s : in out Solution_List; n,d : in integer32;
                 h : in Standard_Complex_VecVecs.VecVec;
                 gamma : in Complex_Number;
                 q : in Poly_Sys; p : out Poly_Sys ) is
  begin
    p := Change_Slices(q,h);
    Standard_Homotopy.Create(p,q,2,gamma); 
    Multitasking_Continuation.Reporting_Multitasking_Path_Tracker(s,n);
    Standard_Homotopy.Clear;
  end Reporting_Sampler;

  procedure Driver_to_Sampler
               ( file : in file_type; n,d : in integer32;
                 ep : in Poly_Sys; esols : in out Solution_List ) is

    target : Poly_Sys(ep'range);
    h : constant Standard_Complex_VecVecs.VecVec(1..d)
      := Witness_Sets.Random_Hyperplanes(natural32(d),natural32(ep'last));
    gamma : constant Complex_Number := Random1;
    ans : character;
    timer : Timing_Widget;

  begin
    Continuation_Parameters.Tune(2);
    put("Do you want to monitor the progress on screen ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y' then
      tstart(timer);
      Reporting_Sampler(esols,n,d,h,gamma,ep,target);
      tstop(timer);
    else
      tstart(timer);
      Silent_Sampler(esols,n,d,h,gamma,ep,target);
      tstop(timer);
    end if;
    put(file,target,esols);
    new_line(file);
    print_times(file,timer,"multitasking continuation");
  end Driver_to_Sampler;

-- AUXILIARIES for monodromy breakup :

  function Keys ( k : Standard_Complex_Vectors.Vector; s : Solution_List )
                return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector of inner products of the vector in key 
  --   with the solutions in s.

    res : Standard_Complex_Vectors.Vector(1..integer32(Length_Of(s)));
    tmp : Solution_List := s;
    ls : Link_to_Solution;

    use Standard_Complex_Vectors;

  begin
    for i in res'range loop
      ls := Head_Of(tmp);
      res(i) := k*ls.v(k'range);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Keys;

  function Map ( k : Standard_Complex_Vectors.Vector;
                 s1,s2 : Solution_List; tol : double_float )
               return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the permutation connecting s1 with s2,
  --   using the random vector in k to get to the keys.

  -- REQUIRED : Length_Of(s1) = Length_Of(s2).

    dim : constant integer32 := integer32(Length_Of(s1));
    res : Standard_Natural_Vectors.Vector(1..dim);
    k1 : constant Standard_Complex_Vectors.Vector(1..dim) := Keys(k,s1);
    k2 : constant Standard_Complex_Vectors.Vector(1..dim) := Keys(k,s2);

  begin
    res := Monodromy_Partitions.Map(k1,k2,tol);
    return res;
  end Map;

  procedure Update_Decomposition
               ( deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                 nbfac : in out integer32;
                 k : in Standard_Complex_Vectors.Vector;
                 s1,s2 : in Solution_List; tol : double_float ) is

  -- DESCRIPTION :
  --   Updates the current decomposition with two solution lists
  --   obtained after one loop.

    p : Standard_Natural_Vectors.Vector(1..integer32(Length_Of(s1)))
      := Map(k,s1,s2,tol);

  begin
    put_line("the permutation : "); put(p); new_line;
    Monodromy_Partitions.Add_Map(deco,natural32(nbfac),p);
    put("number of irreducible factors : "); put(nbfac,1); put_line(" :");
    Monodromy_Partitions.Write_Factors(standard_output,deco.all);
  end Update_Decomposition;

  procedure Driver_to_Monodromy
               ( file : in file_type; n,d : in integer32;
                 ep : in Poly_Sys; esols : in out Solution_List ) is

    tol : constant double_float := 1.0E-8;
    target_alpha,target_beta : Poly_Sys(ep'range);
    s : Standard_Complex_VecVecs.VecVec(1..d)
      := Witness_Sets.Slices(ep,natural32(d));
    h : Standard_Complex_VecVecs.VecVec(1..d);
    dim : constant integer32 := ep'last-d;
    key : Standard_Complex_Vectors.Vector(1..dim)
        := Standard_Random_Vectors.Random_Vector(1,dim);
    alpha,beta : Complex_Number;
    timer : Timing_Widget;
    sols0 : Solution_List;
    degree : constant natural32 := Length_Of(esols);
    deco : Standard_Natural_VecVecs.Link_to_VecVec
         := Monodromy_Partitions.Init_Factors(degree);
    nbfac : integer32 := integer32(degree);

  begin
    Copy(esols,sols0);
    Continuation_Parameters.Tune(2);
    tstart(timer);
    for i in 1..10 loop
      -- put_line("going forward ...");
       h := Witness_Sets.Random_Hyperplanes(natural32(d),natural32(ep'last));
       alpha := Random1; beta := Random1;
       Silent_Sampler(esols,n,d,h,alpha,ep,target_alpha);
      -- put(file,target_alpha,esols);
      -- put_line("going back ...");
       Silent_Sampler(esols,n,d,s,beta,target_alpha,target_beta);
      -- put(file,target_beta,esols);
       Update_Decomposition(deco,nbfac,key,sols0,esols,tol);
       Copy(sols0,esols);
    end loop;
    tstop(timer);
    new_line(file);
    print_times(file,timer,"multitasking monodromy");
  end Driver_to_Monodromy;

end Multitasking_Sampling;
