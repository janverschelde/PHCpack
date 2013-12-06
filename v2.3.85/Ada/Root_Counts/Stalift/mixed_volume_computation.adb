with Standard_Integer_Vectors_io;        use Standard_Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;

with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Transforming_Integer32_Vector_Lists;
 use Transforming_Integer32_Vector_Lists;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Integer_Lifting_Functions;          use Integer_Lifting_Functions;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Integer_Pruning_Methods;            use Integer_Pruning_Methods;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Floating_Integer_Convertors;        use Floating_Integer_Convertors;

-- to fix an "Program received signal SIGFPE, Arithmetic exception."
with Standard_Integer64_Matrices;
with Standard_Integer64_Linear_Solvers;
-- to fix Program received signal SIGFPE, Arithmetic exception.
-- 0x000000010010cafb in standard_integer64_linear_solvers.det ()
with Multprec_Integer_Numbers;
with Multprec_Integer_Matrices;
with Multprec_Integer_Linear_Solvers;
-- but even that does not seem to help!
with Standard_Floating_Numbers;
with Standard_Floating_Matrices;
with Standard_Floating_Linear_Solvers;
-- for debugging:
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;

package body Mixed_Volume_Computation is

-- AUXILIAIRY OUTPUT ROUTINES :

  procedure put ( file : in file_type; points : in Array_of_Lists;
                  n : in natural32; mix : in Vector;
                  mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                  mv : out natural32;
                  multprec_hermite : in boolean := false ) is
  begin
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    new_line(file);
    put(file,points);
    new_line(file);
    put_line(file,"THE MIXED SUBDIVISION :");
    new_line(file);
    put(file,n,mix,mixsub,mv); --,multprec_hermite);
  end put;

  procedure Sort ( supports : in out Array_of_Lists; k,nb,n : in integer32;
                   mxt,perm : in out Vector ) is

  -- DESCRIPTION :
  --   Auxiliary operation for Compute_Mixture.
  --   Compares the kth support with the following supports.
  --   Already nb different supports have been found.

  begin
    for l in (k+1)..n loop
      if Equal(Supports(k),Supports(l)) then
        if l /= k + mxt(nb) then
          declare
            pos : constant integer32 := k + mxt(nb);
            tmpdl : constant List := supports(l);
            tmppos : integer32;
          begin
            supports(l) := supports(pos);
            supports(pos) := tmpdl;
            tmppos := perm(l);
            perm(l) := perm(pos);
            perm(pos) := tmppos;
          end;
        end if;
        mxt(nb) := mxt(nb) + 1;
      end if;
    end loop;
  end Sort;

-- COMPUTING MIXED SUBDIVISION :

  procedure Mixed_Coherent_Subdivision
               ( n : in integer32; mix : in Vector; points : in Array_of_Lists;
                 linear : in boolean; low,upp : in Vector;
                 lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Standard_Floating_Vectors.Vector;
                 mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision ) is

    fa : Array_of_Faces(mix'range);
    index : integer32 := points'first;

  begin
    for k in lifted'range loop                        -- compute lower faces
      if linear
       then lifted(k) := Random_Linear_Lift(low(k),upp(k),points(index));
       else lifted(k) := Random_Lift(low(k),upp(k),points(index));
      end if;
      fa(k) := Create_Lower(mix(k),integer32(n)+1,lifted(k));
      index := index + mix(k);
    end loop;
    Create_CS(n,mix,fa,lifted,nbsucc,nbfail,mixsub); -- prune for mixed cells
    Shallow_Clear(fa);
  end Mixed_Coherent_Subdivision;

-- TARGET ROUTINES :

  procedure Compute_Mixture ( supports : in out Array_of_Lists;
                              mix,perms : out Link_to_Vector ) is

    n : constant integer32 := supports'last;
    cnt : integer32 := 0;           -- counts the number of different supports
    mxt : Vector(supports'range)    -- counts the number of occurrencies
        := (supports'range => 1);
    perm : constant Link_to_Vector  -- keeps track of the permutations
         := new Standard_Integer_Vectors.Vector(supports'range);
    index : integer32 := supports'first;

  begin
    for k in perm'range loop
      perm(k) := k;
    end loop;
    while index <= supports'last loop
      cnt := cnt + 1;
      Sort(supports,index,cnt,n,mxt,perm.all);
      index := index + mxt(cnt);
    end loop;
    mix := new Standard_Integer_Vectors.Vector'(mxt(mxt'first..cnt));
    perms := perm;
  end Compute_Mixture;

  function Compute_Index ( k : integer32; mix : Vector ) return integer32 is

  -- DESCRIPTION :
  --   Returns the index of k w.r.t. to the type of mixture.

    index : integer32 := mix(mix'first);

  begin
    if k <= index then
      return mix'first;
    else
      for l in (mix'first+1)..mix'last loop
        index := index + mix(l);
        if k <= index
         then return l;
        end if;
      end loop;
      return mix'last;
    end if;
  end Compute_Index;

  function Compute_Permutation ( n : integer32; mix : Vector;
                                 supports : Array_of_Lists )
                               return Link_to_Vector is

    perms : constant Link_to_Vector := new Vector(1..n);

  begin
    for k in perms'range loop
      perms(k) := k;
    end loop;
    return perms;
  end Compute_Permutation;

  function Permute ( p : Poly_Sys; perm : Link_to_Vector ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := p(perm(k));
    end loop;
    return res;
  end Permute;

  function Permute ( p : Laur_Sys; perm : Link_to_Vector ) return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := p(perm(k));
    end loop;
    return res;
  end Permute;

  function Permute ( supports : Array_of_Lists; perm : Link_to_Vector )
                   return Array_of_Lists is

    res : Array_of_Lists(supports'range);

  begin
    for k in supports'range loop
      res(k) := supports(perm(k));
    end loop;
    return res;
  end Permute;

  function Typed_Lists ( mix : Vector; points : Array_of_Lists )
                       return Array_of_Lists is

    res : Array_of_Lists(mix'range);
    ind : integer32 := res'first;

  begin
    for i in mix'range loop
      res(i) := points(ind);
      ind := ind + mix(i);
    end loop;
    return res;
  end Typed_Lists;

-- MIXED VOLUME COMPUTATIONS BASED ON SUBDIVISIONS :

-- AUXILIARIES :

  function Is_Fine ( mix : Vector;
                     mic : Integer_Mixed_Subdivisions.Mixed_Cell )
                   return boolean is

  -- DESCRIPTION :
  --   Returns true if the mixed volume can be computed by a determinant.

    fine : boolean := true;

  begin
    for k in mic.pts'range loop
      fine := (Length_Of(mic.pts(k)) = natural32(mix(k)) + 1);
      exit when not fine;
    end loop;
    return fine;
  end Is_Fine;

  function Reduced_Supports
              ( n : integer32; mix : Vector;
                mic : Integer_Mixed_Subdivisions.Mixed_Cell )
              return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the supports of the cell without the lifting values.

    res : Array_of_Lists(1..n);
    cnt : integer32 := 1;

  begin
    for k in mic.pts'range loop
      res(cnt) := Reduce(mic.pts(k),n+1);
      for l in 1..mix(k)-1 loop
        Copy(res(cnt),res(cnt+l));
      end loop;
      cnt := cnt + mix(k);
    end loop;
    return res;
  end Reduced_Supports;

  function Determinant
              ( A : Standard_Integer_Matrices.Matrix ) return integer32 is

    use Standard_Floating_Numbers;
    B : Standard_Floating_Matrices.Matrix(A'range(1),A'range(2));
    det : double_float;
    res : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := double_float(A(i,j));
      end loop;
    end loop;
    Standard_Floating_Linear_Solvers.Triangulate(B,B'last(1),B'last(2));
    det := 1.0;
    for i in B'range(1) loop
      det := det*B(i,i);
    end loop;
    res := integer32(det);
    return res;
  end Determinant;

  function Det64 ( A : Standard_Integer_Matrices.Matrix ) return integer32 is

  -- DESCRIPTION :
  --   Computes the determinant using 64-bit integer arithmetic.

    B : Standard_Integer64_Matrices.Matrix(A'range(1),A'range(2));
    res : integer64;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := integer64(A(i,j));
      end loop;
    end loop;
    res := Standard_Integer64_Linear_Solvers.Det(B);
    return integer32(res);
  end Det64;

  function DetMP ( A : Standard_Integer_Matrices.Matrix ) return integer32 is

  -- DESCRIPTION :
  --   Computes the determinant using multiprecision integer arithmetic.

    B : Multprec_Integer_Matrices.Matrix(A'range(1),A'range(2));
    det : Multprec_Integer_Numbers.Integer_Number;
    res : integer32;

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        B(i,j) := Multprec_Integer_Numbers.Create(A(i,j));
      end loop;
    end loop;
    det := Multprec_Integer_Linear_Solvers.Det(B);
    res := Multprec_Integer_Numbers.create(det);
    Multprec_Integer_Matrices.Clear(B);
    Multprec_Integer_Numbers.Clear(det);
    return res;
  end DetMP;

  function Fine_Mixed_Volume
              ( n : integer32; mix : Vector;
                mic : Integer_Mixed_Subdivisions.Mixed_Cell;
                multprec_hermite : boolean := false )
              return natural32 is

  -- DESCRIPTION :
  --   Computes the mixed volume for a cell that is fine mixed.

  -- REQUIRED : Fine(mix,mic).

    res : natural32;
    mat : matrix(1..n,1..n);
    count,detmat : integer32;
    tmp : List;
    sh,pt : Link_to_Vector;

  begin
    count := 1;
    for k in mic.pts'range loop
      sh := Head_Of(mic.pts(k));
      tmp := Tail_Of(mic.pts(k));
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        for j in 1..n loop
          mat(count,j) := pt(j) - sh(j);
        end loop;
        tmp := Tail_Of(tmp);
        count := count + 1;
      end loop;
    end loop;
    if multprec_hermite then
     -- put_line("using multiprecision Hermite normal form ...");
      detmat := DetMP(mat); -- Det64(mat); -- Det(mat);
    else
     -- put_line("using standard arithmetic ...");
       declare
       begin
         detmat := Det(mat);
       exception
         when others
           => put_line("exception raised when computing determinant");
              detmat := Determinant(mat);
             put("floating point computation gives "); put(detmat,1); new_line;
             put_line("but multiprecision Hermite normal form will be sure...");
       end;
    end if;
   -- put("correct determinant gives "); put(detmat,1); new_line;
    if detmat >= 0
     then res := natural32(detmat);
     else res := natural32(-detmat);
    end if;
    return res;
  end Fine_Mixed_Volume;

  function Mixed_Volume 
              ( n : integer32; mix : Vector;
                mic : Integer_Mixed_Subdivisions.Mixed_Cell;
                multprec_hermite : in boolean := false )
              return natural32 is

  -- ALGORITHM :
  --   First check if the cell has a refinement, if so, then use it,
  --   if not, then check if the cell is fine mixed.
  --   If the cell is fine mixed, only a determinant needs to be computed,
  --   otherwise the cell will be refined.
   
    use Integer_Mixed_Subdivisions;
    res : natural32;

  begin
    if (mic.sub /= null) and then not Is_Null(mic.sub.all) then
      res := Mixed_Volume_Computation.Mixed_Volume
                (n,mix,mic.sub.all,multprec_hermite);
    elsif Is_Fine(mix,mic) then
      res := Fine_Mixed_Volume(n,mix,mic,multprec_hermite);
    else
      declare
        rcell : Array_of_Lists(1..n) := Reduced_Supports(n,mix,mic);
      begin
        res := Mixed_Volume_Computation.Mixed_Volume(n,rcell);
        Deep_Clear(rcell);
      end;
    end if;
    return res;
  end Mixed_Volume;

  function Mixed_Volume
              ( n : integer32; mix : Vector;
                mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
                multprec_hermite : boolean := false )
              return natural32 is

    use Integer_Mixed_Subdivisions;

    res : natural32 := 0;
    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      res := res + Mixed_Volume(n,mix,Head_Of(tmp),multprec_hermite);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Mixed_Volume;

  procedure Mixed_Volume
               ( n : in integer32; mix : in Vector;
                 mic : in out Integer_Mixed_Subdivisions.Mixed_Cell;
                 mv : out natural32;
                 multprec_hermite : in boolean := false ) is

    use Integer_Mixed_Subdivisions;

  begin
    if (mic.sub /= null) and then not Is_Null(mic.sub.all) then
      mv := Mixed_Volume_Computation.Mixed_Volume
              (n,mix,mic.sub.all,multprec_hermite);
    elsif Is_Fine(mix,mic) then
      mv := Fine_Mixed_Volume(n,mix,mic,multprec_hermite);
    else -- NOTE : keep the same type of mixture!
      declare
        rcell : Array_of_Lists(1..n) := Reduced_Supports(n,mix,mic);
        lifted : Array_of_Lists(mix'range);
        mixsub : Mixed_Subdivision;
      begin
        Mixed_Volume_Computation.Mixed_Volume
          (n,mix,rcell,lifted,mixsub,mv,multprec_hermite);
        mic.sub := new Mixed_Subdivision'(mixsub);
        Deep_Clear(rcell);  Deep_Clear(lifted);
      end;
    end if;
  end Mixed_Volume;

  procedure Mixed_Volume
               ( n : in integer32; mix : in Vector;
                 mic : in out Floating_Mixed_Subdivisions.Mixed_Cell;
                 mv : out natural32;
                 multprec_hermite : in boolean := false ) is

    intmic : Integer_Mixed_Subdivisions.Mixed_Cell;
    intsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    use Integer_Mixed_Subdivisions;
    use Floating_Mixed_Subdivisions;

  begin
    if mic.sub /= null then
      intsub := Convert(mic.sub.all);
      Mixed_Volume(n,mix,intsub,mv,multprec_hermite);
      Deep_Clear(intsub);
    else
      intmic := Convert(mic);
      Mixed_Volume(n,mix,intmic,mv,multprec_hermite);
      if intmic.sub /= null
       then mic.sub := new Floating_Mixed_Subdivisions.
                               Mixed_Subdivision'(Convert(intmic.sub.all));
      end if;
      Deep_Clear(intmic);
    end if;
  end Mixed_Volume;

  procedure Mixed_Volume
              ( n : in integer32; mix : in Vector;
                mixsub : in out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32;
                multprec_hermite : in boolean := false ) is

    use Integer_Mixed_Subdivisions;
    res : natural32 := 0;
    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        mmv : natural32;
      begin
        Mixed_Volume(n,mix,mic,mmv,multprec_hermite);
        Set_Head(tmp,mic);
        res := res + mmv;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    mv := res;
  end Mixed_Volume;

  procedure Mixed_Volume
              ( n : in integer32; mix : in Vector;
                mixsub : in out Floating_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32;
                multprec_hermite : in boolean := false ) is

    use Floating_Mixed_Subdivisions;
    res : natural32 := 0;
    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        mmv : natural32;
      begin
        Mixed_Volume(n,mix,mic,mmv,multprec_hermite);
        Set_Head(tmp,mic);
        res := res + mmv;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    mv := res;
  end Mixed_Volume;

-- MIXED VOLUME COMPUTATIONS BASED ON SUPPORTS :

  function Mixed_Volume ( n : integer32; supports : Array_of_Lists )
			return natural32 is

    mv : natural32;
    mix,perm : Link_to_Vector;
    permsupp : Array_of_Lists(supports'range);

  begin
    Copy(supports,permsupp);
    Compute_Mixture(permsupp,mix,perm);
    mv := Mixed_Volume(n,mix.all,permsupp);
    Clear(mix); Clear(perm); Deep_Clear(permsupp);
    return mv;
  end Mixed_Volume;

  function Mixed_Volume ( file : file_type; n : integer32;
                          supports : Array_of_Lists ) return natural32 is

    mv : natural32;
    mix,perm : Link_to_Vector;
    permsupp : Array_of_Lists(supports'range);

  begin
    Copy(supports,permsupp);
    Compute_Mixture(permsupp,mix,perm);
    mv := Mixed_Volume(file,n,mix.all,permsupp);
    Clear(mix); Clear(perm); Deep_Clear(permsupp);
    return mv;
  end Mixed_Volume;

  function Mixed_Volume ( n : integer32; mix : Vector;
                          supports : Array_of_Lists ) return natural32 is

    res : natural32;
    mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    lifted : Array_of_Lists(mix'range);

  begin
    Mixed_Volume_Computation.Mixed_Volume(n,mix,supports,lifted,mixsub,res);
    Deep_Clear(lifted);
    Integer_Mixed_Subdivisions.Shallow_Clear(mixsub);
    return res;
  end Mixed_Volume;

  procedure Mixed_Volume 
              ( n : in integer32; mix : in Vector;
                supports : in Array_of_Lists; lifted : out Array_of_Lists;
                mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32;
                multprec_hermite : in boolean := false ) is

    low : constant Vector := (mix'range => 0);
    upp : constant Vector := Adaptive_Lifting(supports);
    nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                  := (mix'range => 0.0);
    liftsupp : Array_of_Lists(mix'range);
    sub : Integer_Mixed_Subdivisions.Mixed_Subdivision;

  begin
    Mixed_Coherent_Subdivision
      (n,mix,supports,false,low,upp,liftsupp,nbsucc,nbfail,sub);
    Mixed_Volume(n,mix,sub,mv,multprec_hermite);
    lifted := liftsupp;
    mixsub := sub;
  end Mixed_Volume;

  function Mixed_Volume ( file : file_type; n : integer32; mix : Vector;
                          supports : Array_of_Lists ) return natural32 is
  
    res : natural32;
    mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    lifted : Array_of_Lists(mix'range);

  begin
    Mixed_Volume_Computation.Mixed_Volume
      (file,n,mix,supports,lifted,mixsub,res);
    Deep_Clear(lifted);
    Integer_Mixed_Subdivisions.Shallow_Clear(mixsub);
    return res;
  end Mixed_Volume;

  procedure Mixed_Volume
              ( file : in file_type; n : in integer32; mix : in Vector;
                supports : in Array_of_Lists; lifted : out Array_of_Lists;
                mixsub : out Integer_Mixed_Subdivisions.Mixed_Subdivision;
                mv : out natural32;
                multprec_hermite : in boolean := false ) is

    low : constant Vector := (mix'range => 0);
    upp : constant Vector := Adaptive_Lifting(supports);
    sub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    nbsucc,nbfail : Standard_Floating_Vectors.Vector(mix'range)
                  := (mix'range => 0.0);
    liftsupp : Array_of_Lists(mix'range);

  begin
    Mixed_Coherent_Subdivision
      (n,mix,supports,false,low,upp,liftsupp,nbsucc,nbfail,sub);
    lifted := liftsupp;
    mixsub := sub;
    put(file,liftsupp,natural32(n),mix,sub,mv,multprec_hermite);
  end Mixed_Volume;

end Mixed_Volume_Computation;
