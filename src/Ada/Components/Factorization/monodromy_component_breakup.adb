with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with DoblDobl_Random_Vectors;           use DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;           use QuadDobl_Random_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Witness_Sets;                      use Witness_Sets;
with Sampling_Machine; 
with Sampling_Laurent_Machine; 
with DoblDobl_Sampling_Machine; 
with DoblDobl_Sampling_Laurent_Machine; 
with QuadDobl_Sampling_Machine; 
with QuadDobl_Sampling_Laurent_Machine; 
with Sample_Points;                     use Sample_Points;
with Sample_Point_Grids;                use Sample_Point_Grids;
with DoblDobl_Sample_Points;            use DoblDobl_Sample_Points;
with DoblDobl_Sample_Grids;                use DoblDobl_Sample_Grids;
with QuadDobl_Sample_Points;            use QuadDobl_Sample_Points;
with QuadDobl_Sample_Grids;                use QuadDobl_Sample_Grids;
with Rectangular_Sample_Grids;          use Rectangular_Sample_Grids;
with DoblDobl_Rectangular_Sample_Grids; use DoblDobl_Rectangular_Sample_Grids;
with QuadDobl_Rectangular_Sample_Grids; use QuadDobl_Rectangular_Sample_Grids;
with Drivers_to_Grid_Creators;          use Drivers_to_Grid_Creators;
with Monodromy_Partitions;              use Monodromy_Partitions;
with Standard_Trace_Interpolators;      use Standard_Trace_Interpolators;
with DoblDobl_Trace_Interpolators;      use DoblDobl_Trace_Interpolators;
with QuadDobl_Trace_Interpolators;      use QuadDobl_Trace_Interpolators;

package body Monodromy_Component_Breakup is

-- INTERNAL STATE :

  use_laurent : boolean := false;

-- AUXILIARY ROUTINES FOR MONODROMY :

  function Keys ( k : Standard_Complex_Vectors.Vector;
                  len : in integer32; sps : Standard_Sample_List )
                return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   A list of keys is obtained from the inner product of the vector v
  --   with all solutions (length = len) in the list of samples.

    res : Standard_Complex_Vectors.Vector(1..len);
    tmp : Standard_Sample_List := sps;
    use Standard_Complex_Vectors;

  begin
    for i in res'range loop
      declare
        s : constant Standard_Sample := Head_Of(tmp);
      begin
        res(i) := k*Sample_Point(s).v(k'range); 
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Keys;

  function Keys ( k : DoblDobl_Complex_Vectors.Vector;
                  len : in integer32; sps : DoblDobl_Sample_List )
                return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   A list of keys is obtained from the inner product of the vector v
  --   with all solutions (length = len) in the list of samples.

    res : DoblDobl_Complex_Vectors.Vector(1..len);
    tmp : DoblDobl_Sample_List := sps;
    use DoblDobl_Complex_Vectors;

  begin
    for i in res'range loop
      declare
        s : constant DoblDobl_Sample := Head_Of(tmp);
      begin
        res(i) := k*Sample_Point(s).v(k'range); 
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Keys;

  function Keys ( k : QuadDobl_Complex_Vectors.Vector;
                  len : in integer32; sps : QuadDobl_Sample_List )
                return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   A list of keys is obtained from the inner product of the vector v
  --   with all solutions (length = len) in the list of samples.

    res : QuadDobl_Complex_Vectors.Vector(1..len);
    tmp : QuadDobl_Sample_List := sps;
    use QuadDobl_Complex_Vectors;

  begin
    for i in res'range loop
      declare
        s : constant QuadDobl_Sample := Head_Of(tmp);
      begin
        res(i) := k*Sample_Point(s).v(k'range); 
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Keys;

  procedure Update ( f : in Standard_Natural_VecVecs.Link_to_VecVec;
                     nf,cnt,nit : in out natural32; tol : in double_float;
                     k0,rvk : in Standard_Complex_Vectors.Vector;
                     len : in integer32; sps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the permutation map, updates the factorization and
  --   maintains the counters.

  -- ON ENTRY :
  --   f             current factorization;
  --   nf            number of factors in f;
  --   cnt           counts the number of times f did not change;
  --   nit           counts the number of monodromy loops;
  --   tol           tolerance to decide equality of floats;
  --   k0            keys for the samples at the start;
  --   rvk           random vector used to make keys;
  --   len           number of samples in the list sps;
  --   sps           new list of samples at the arrival of one loop.

  -- ON RETURN :
  --   f             current factorization updated using the loop;
  --   nf            number of factors in f;
  --   cnt           updated counter of number of times f did not change;
  --   nit           augmented with one.

    k1 : constant Standard_Complex_Vectors.Vector(1..len) := Keys(rvk,len,sps);
    mp : constant Standard_Natural_Vectors.Vector(1..len) := Map(k0,k1,tol);
    n0 : constant natural32 := nf;

  begin
    nit := nit + 1;
    Add_Map(f,nf,mp);
    if n0 = nf
     then cnt := cnt + 1;
     else cnt := 0;
    end if;
  end Update;

  procedure Update ( f : in Standard_Natural_VecVecs.Link_to_VecVec;
                     nf,cnt,nit : in out natural32; tol : in double_float;
                     k0,rvk : in DoblDobl_Complex_Vectors.Vector;
                     len : in integer32; sps : in DoblDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the permutation map, updates the factorization and
  --   maintains the counters.

  -- ON ENTRY :
  --   f             current factorization;
  --   nf            number of factors in f;
  --   cnt           counts the number of times f did not change;
  --   nit           counts the number of monodromy loops;
  --   tol           tolerance to decide equality of floats;
  --   k0            keys for the samples at the start;
  --   rvk           random vector used to make keys;
  --   len           number of samples in the list sps;
  --   sps           new list of samples at the arrival of one loop.

  -- ON RETURN :
  --   f             current factorization updated using the loop;
  --   nf            number of factors in f;
  --   cnt           updated counter of number of times f did not change;
  --   nit           augmented with one.

    k1 : constant DoblDobl_Complex_Vectors.Vector(1..len) := Keys(rvk,len,sps);
    mp : constant Standard_Natural_Vectors.Vector(1..len) := Map(k0,k1,tol);
    n0 : constant natural32 := nf;

  begin
    nit := nit + 1;
    Add_Map(f,nf,mp);
    if n0 = nf
     then cnt := cnt + 1;
     else cnt := 0;
    end if;
  end Update;

  procedure Update ( f : in Standard_Natural_VecVecs.Link_to_VecVec;
                     nf,cnt,nit : in out natural32; tol : in double_float;
                     k0,rvk : in QuadDobl_Complex_Vectors.Vector;
                     len : in integer32; sps : in QuadDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the permutation map, updates the factorization and
  --   maintains the counters.

  -- ON ENTRY :
  --   f             current factorization;
  --   nf            number of factors in f;
  --   cnt           counts the number of times f did not change;
  --   nit           counts the number of monodromy loops;
  --   tol           tolerance to decide equality of floats;
  --   k0            keys for the samples at the start;
  --   rvk           random vector used to make keys;
  --   len           number of samples in the list sps;
  --   sps           new list of samples at the arrival of one loop.

  -- ON RETURN :
  --   f             current factorization updated using the loop;
  --   nf            number of factors in f;
  --   cnt           updated counter of number of times f did not change;
  --   nit           augmented with one.

    k1 : constant QuadDobl_Complex_Vectors.Vector(1..len) := Keys(rvk,len,sps);
    mp : constant Standard_Natural_Vectors.Vector(1..len) := Map(k0,k1,tol);
    n0 : constant natural32 := nf;

  begin
    nit := nit + 1;
    Add_Map(f,nf,mp);
    if n0 = nf
     then cnt := cnt + 1;
     else cnt := 0;
    end if;
  end Update;

  procedure Update ( file : in file_type;
                     f : in Standard_Natural_VecVecs.Link_to_VecVec;
                     nf,cnt,nit : in out natural32; tol : in double_float;
                     k0,rvk : in Standard_Complex_Vectors.Vector;
                     len : in integer32; sps : in Standard_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the permutation map, updates the factorization and
  --   maintains the counters.

  -- ON ENTRY :
  --   file          for intermediate output and diagnostics;
  --   f             current factorization;
  --   nf            number of factors in f;
  --   cnt           counts the number of times f did not change;
  --   nit           counts the number of monodromy loops;
  --   tol           tolerance to decide equality of floats;
  --   k0            keys for the samples at the start;
  --   rvk           random vector used to make keys;
  --   len           number of samples in the list sps;
  --   sps           new list of samples at the arrival of one loop.

  -- ON RETURN :
  --   f             current factorization updated using the loop;
  --   nf            number of factors in f;
  --   cnt           updated counter of number of times f did not change;
  --   nit           augmented with one.

    k1 : constant Standard_Complex_Vectors.Vector(1..len) := Keys(rvk,len,sps);
    mp : constant Standard_Natural_Vectors.Vector(1..len) := Map(k0,k1,tol);
    n0 : constant natural32 := nf;

  begin
    nit := nit + 1;
    put(file,"Monodromy Permutation #"); put(file,nit,1);
    put(file," : "); put(file,mp); new_line(file);
    Add_Map(f,nf,mp);
    put_line(file,"Current factorization : ");
    Write_Factors(file,f.all);
    put(file,"Drop in #factors : ");
    put(file,n0,1); put(file," -> "); put(file,nf,1);
    if n0 = nf
     then cnt := cnt + 1;
     else cnt := 0;
    end if;
    put(file,"  count : "); put(file,cnt,1); new_line(file);
  end Update;

  procedure Update ( file : in file_type;
                     f : in Standard_Natural_VecVecs.Link_to_VecVec;
                     nf,cnt,nit : in out natural32; tol : in double_float;
                     k0,rvk : in DoblDobl_Complex_Vectors.Vector;
                     len : in integer32; sps : in DoblDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the permutation map, updates the factorization and
  --   maintains the counters.

  -- ON ENTRY :
  --   file          for intermediate output and diagnostics;
  --   f             current factorization;
  --   nf            number of factors in f;
  --   cnt           counts the number of times f did not change;
  --   nit           counts the number of monodromy loops;
  --   tol           tolerance to decide equality of floats;
  --   k0            keys for the samples at the start;
  --   rvk           random vector used to make keys;
  --   len           number of samples in the list sps;
  --   sps           new list of samples at the arrival of one loop.

  -- ON RETURN :
  --   f             current factorization updated using the loop;
  --   nf            number of factors in f;
  --   cnt           updated counter of number of times f did not change;
  --   nit           augmented with one.

    k1 : constant DoblDobl_Complex_Vectors.Vector(1..len) := Keys(rvk,len,sps);
    mp : constant Standard_Natural_Vectors.Vector(1..len) := Map(k0,k1,tol);
    n0 : constant natural32 := nf;

  begin
    nit := nit + 1;
    put(file,"Monodromy Permutation #"); put(file,nit,1);
    put(file," : "); put(file,mp); new_line(file);
    Add_Map(f,nf,mp);
    put_line(file,"Current factorization : ");
    Write_Factors(file,f.all);
    put(file,"Drop in #factors : ");
    put(file,n0,1); put(file," -> "); put(file,nf,1);
    if n0 = nf
     then cnt := cnt + 1;
     else cnt := 0;
    end if;
    put(file,"  count : "); put(file,cnt,1); new_line(file);
  end Update;

  procedure Update ( file : in file_type;
                     f : in Standard_Natural_VecVecs.Link_to_VecVec;
                     nf,cnt,nit : in out natural32; tol : in double_float;
                     k0,rvk : in QuadDobl_Complex_Vectors.Vector;
                     len : in integer32; sps : in QuadDobl_Sample_List ) is

  -- DESCRIPTION :
  --   Computes the permutation map, updates the factorization and
  --   maintains the counters.

  -- ON ENTRY :
  --   file          for intermediate output and diagnostics;
  --   f             current factorization;
  --   nf            number of factors in f;
  --   cnt           counts the number of times f did not change;
  --   nit           counts the number of monodromy loops;
  --   tol           tolerance to decide equality of floats;
  --   k0            keys for the samples at the start;
  --   rvk           random vector used to make keys;
  --   len           number of samples in the list sps;
  --   sps           new list of samples at the arrival of one loop.

  -- ON RETURN :
  --   f             current factorization updated using the loop;
  --   nf            number of factors in f;
  --   cnt           updated counter of number of times f did not change;
  --   nit           augmented with one.

    k1 : constant QuadDobl_Complex_Vectors.Vector(1..len) := Keys(rvk,len,sps);
    mp : constant Standard_Natural_Vectors.Vector(1..len) := Map(k0,k1,tol);
    n0 : constant natural32 := nf;

  begin
    nit := nit + 1;
    put(file,"Monodromy Permutation #"); put(file,nit,1);
    put(file," : "); put(file,mp); new_line(file);
    Add_Map(f,nf,mp);
    put_line(file,"Current factorization : ");
    Write_Factors(file,f.all);
    put(file,"Drop in #factors : ");
    put(file,n0,1); put(file," -> "); put(file,nf,1);
    if n0 = nf
     then cnt := cnt + 1;
     else cnt := 0;
    end if;
    put(file,"  count : "); put(file,cnt,1); new_line(file);
  end Update;

-- TARGET ROUTINES :

-- AUXILIARIES FOR LINEAR TRACE CERTIFICATES :

  function Create ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..2);
    sli : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);

  begin
    Rectangular_Sample_Grids.Set_Polynomial_Type(false);
    res := Create1(sps,2);
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..2);
    sli : constant DoblDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant DoblDobl_Sample_List := Create(sols,sli);

  begin
    DoblDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(false);
    res := Create1(sps,2);
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..2);
    sli : constant QuadDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant QuadDobl_Sample_List := Create(sols,sli);

  begin
    QuadDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(false);
    res := Create1(sps,2);
    return res;
  end Create;

  function Create ( file : file_type;
                    p : Standard_Complex_Poly_Systems.Poly_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..2);
    sli : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    eps,dist : double_float;

  begin
    new_line(file);
    put_line(file,"GRID CREATION FOR LINEAR TRACES VALIDATION");
    Rectangular_Sample_Grids.Set_Polynomial_Type(false);
    Standard_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    return res;
  end Create;

  function Create ( file : file_type;
                    p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..2);
    sli : constant DoblDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant DoblDobl_Sample_List := Create(sols,sli);
    eps,dist : double_double;

  begin
    new_line(file);
    put_line(file,"GRID CREATION FOR LINEAR TRACES VALIDATION");
    DoblDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(false);
    DoblDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    return res;
  end Create;

  function Create ( file : file_type;
                    p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..2);
    sli : constant QuadDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant QuadDobl_Sample_List := Create(sols,sli);
    eps,dist : quad_double;

  begin
    new_line(file);
    put_line(file,"GRID CREATION FOR LINEAR TRACES VALIDATION");
    QuadDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(false);
    QuadDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    return res;
  end Create;

  function Create ( p : Standard_Complex_Laur_Systems.Laur_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..2);
    sli : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);

  begin
    Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    res := Create1(sps,2);
    return res;
  end Create;

  function Create ( p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..2);
    sli : constant DoblDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant DoblDobl_Sample_List := Create(sols,sli);

  begin
    DoblDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    res := Create1(sps,2);
    return res;
  end Create;

  function Create ( p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..2);
    sli : constant QuadDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant QuadDobl_Sample_List := Create(sols,sli);

  begin
    QuadDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    res := Create1(sps,2);
    return res;
  end Create;

  function Create ( file : file_type;
                    p : Standard_Complex_Laur_Systems.Laur_Sys;
                    sols : Standard_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_Standard_Sample_Lists is

    res : Array_of_Standard_Sample_Lists(0..2);
    sli : constant Standard_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant Standard_Sample_List := Create(sols,sli);
    eps,dist : double_float;

  begin
    new_line(file);
    put_line(file,"GRID CREATION FOR LINEAR TRACES VALIDATION");
    Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    Standard_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    return res;
  end Create;

  function Create ( file : file_type;
                    p : DoblDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : DoblDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_DoblDobl_Sample_Lists is

    res : Array_of_DoblDobl_Sample_Lists(0..2);
    sli : constant DoblDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant DoblDobl_Sample_List := Create(sols,sli);
    eps,dist : double_double;

  begin
    new_line(file);
    put_line(file,"GRID CREATION FOR LINEAR TRACES VALIDATION");
    DoblDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    DoblDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    return res;
  end Create;

  function Create ( file : file_type;
                    p : QuadDobl_Complex_Laur_Systems.Laur_Sys;
                    sols : QuadDobl_Complex_Solutions.Solution_List;
                    dim : natural32 )
                  return Array_of_QuadDobl_Sample_Lists is

    res : Array_of_QuadDobl_Sample_Lists(0..2);
    sli : constant QuadDobl_Complex_VecVecs.VecVec := Slices(p,dim);
    sps : constant QuadDobl_Sample_List := Create(sols,sli);
    eps,dist : quad_double;

  begin
    new_line(file);
    put_line(file,"GRID CREATION FOR LINEAR TRACES VALIDATION");
    QuadDobl_Rectangular_Sample_Grids.Set_Polynomial_Type(true);
    QuadDobl_Rectangular_Grid_Creator(file,sps,2,res,eps,dist);
    return res;
  end Create;

-- VALIDATION OF BREAKUP WITH LINEAR TRACES :

  function Trace_Sum_Difference
             ( f : Standard_Natural_Vectors.Vector;
               grid : Array_of_Standard_Sample_Lists ) return double_float is

    use Standard_Complex_Numbers;

    sub_grid : constant Array_of_Standard_Sample_Lists(grid'range)
             := Select_Sub_Grid(grid,f);
    q : constant Standard_Complex_Vectors.Vector := Create(sub_grid,1);
    val,eva : Complex_Number;

  begin
    Eval_Trace(q,1,sub_grid(2),val,eva);
    return AbsVal(val-eva);
  end Trace_Sum_Difference;

  function Trace_Sum_Difference
             ( file : file_type;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_Standard_Sample_Lists ) return double_float is

    use Standard_Complex_Numbers;

    sub_grid : constant Array_of_Standard_Sample_Lists(grid'range)
             := Select_Sub_Grid(grid,f);
    q : constant Standard_Complex_Vectors.Vector := Create(sub_grid,1);
    val,eva : Complex_Number;

  begin
    Eval_Trace(q,1,sub_grid(2),val,eva);
    put(file,"calculated sum at samples : ");
    put(file,val); new_line(file);
    put(file,"value at the linear trace : ");
    put(file,eva); new_line(file);
    return AbsVal(val-eva);
  end Trace_Sum_Difference;

  function Trace_Sum_Difference
             ( f : Standard_Natural_Vectors.Vector;
               grid : Array_of_DoblDobl_Sample_Lists ) return double_double is

    use DoblDobl_Complex_Numbers;

    sub_grid : constant Array_of_DoblDobl_Sample_Lists(grid'range)
             := Select_Sub_Grid(grid,f);
    q : constant DoblDobl_Complex_Vectors.Vector := Create(sub_grid,1);
    val,eva : Complex_Number;

  begin
    Eval_Trace(q,1,sub_grid(2),val,eva);
    return AbsVal(val-eva);
  end Trace_Sum_Difference;

  function Trace_Sum_Difference
             ( file : file_type;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_DoblDobl_Sample_Lists ) return double_double is

    use DoblDobl_Complex_Numbers;

    sub_grid : constant Array_of_DoblDobl_Sample_Lists(grid'range)
             := Select_Sub_Grid(grid,f);
    q : constant DoblDobl_Complex_Vectors.Vector := Create(sub_grid,1);
    val,eva : Complex_Number;

  begin
    Eval_Trace(q,1,sub_grid(2),val,eva);
    put(file,"calculated sum at samples : ");
    put(file,val); new_line(file);
    put(file,"value at the linear trace : ");
    put(file,eva); new_line(file);
    return AbsVal(val-eva);
  end Trace_Sum_Difference;

  function Trace_Sum_Difference
             ( f : Standard_Natural_Vectors.Vector;
               grid : Array_of_QuadDobl_Sample_Lists ) return quad_double is

    use QuadDobl_Complex_Numbers;

    sub_grid : constant Array_of_QuadDobl_Sample_Lists(grid'range)
             := Select_Sub_Grid(grid,f);
    q : constant QuadDobl_Complex_Vectors.Vector := Create(sub_grid,1);
    val,eva : Complex_Number;

  begin
    Eval_Trace(q,1,sub_grid(2),val,eva);
    return AbsVal(val-eva);
  end Trace_Sum_Difference;

  function Trace_Sum_Difference
             ( file : file_type;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_QuadDobl_Sample_Lists ) return quad_double is

    use QuadDobl_Complex_Numbers;

    sub_grid : constant Array_of_QuadDobl_Sample_Lists(grid'range)
             := Select_Sub_Grid(grid,f);
    q : constant QuadDobl_Complex_Vectors.Vector := Create(sub_grid,1);
    val,eva : Complex_Number;

  begin
    Eval_Trace(q,1,sub_grid(2),val,eva);
    put(file,"calculated sum at samples : ");
    put(file,val); new_line(file);
    put(file,"value at the linear trace : ");
    put(file,eva); new_line(file);
    return AbsVal(val-eva);
  end Trace_Sum_Difference;

  function Certify_Factor
             ( tol : double_float;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_Standard_Sample_Lists ) return boolean is
  begin
    return (Trace_Sum_Difference(f,grid) < tol);
  end Certify_Factor;

  function Certify_Factor
             ( file : file_type; tol : double_float;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_Standard_Sample_Lists ) return boolean is

    d : constant double_float := Trace_Sum_Difference(file,f,grid);
 
  begin
    put(file,"The witness points"); put(file,f);
    if d < tol then
      put_line(file," define a factor.");
      return true;
    else
      put_line(file," do not define a factor.");
      return false;
    end if;
  end Certify_Factor;

  function Certify_Factor
             ( tol : double_float;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_DoblDobl_Sample_Lists ) return boolean is
  begin
    return (Trace_Sum_Difference(f,grid) < tol);
  end Certify_Factor;

  function Certify_Factor
             ( file : file_type; tol : double_float;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_DoblDobl_Sample_Lists ) return boolean is

    d : constant double_double := Trace_Sum_Difference(file,f,grid);
 
  begin
    put(file,"The witness points"); put(file,f);
    if d < tol then
      put_line(file," define a factor.");
      return true;
    else
      put_line(file," do not define a factor.");
      return false;
    end if;
  end Certify_Factor;

  function Certify_Factor
             ( tol : double_float;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_QuadDobl_Sample_Lists ) return boolean is
  begin
    return (Trace_Sum_Difference(f,grid) < tol);
  end Certify_Factor;

  function Certify_Factor
             ( file : file_type; tol : double_float;
               f : Standard_Natural_Vectors.Vector;
               grid : Array_of_QuadDobl_Sample_Lists ) return boolean is

    d : constant quad_double := Trace_Sum_Difference(file,f,grid);
 
  begin
    put(file,"The witness points"); put(file,f);
    if d < tol then
      put_line(file," define a factor.");
      return true;
    else
      put_line(file," do not define a factor.");
      return false;
    end if;
  end Certify_Factor;

  function Is_Factorization
                ( tol : double_float;
                  f : Standard_Natural_VecVecs.VecVec;
                  grid : Array_of_Standard_Sample_Lists ) return boolean is

    use Standard_Natural_Vectors;

  begin
    for i in f'range loop
      if f(i) /= null then
        if not Certify_Factor(tol,f(i).all,grid)
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Is_Factorization;

  function Is_Factorization
                ( file : file_type; tol : double_float;
                  f : Standard_Natural_VecVecs.VecVec;
                  grid : Array_of_Standard_Sample_Lists ) return boolean is

    use Standard_Natural_Vectors;

  begin
    for i in f'range loop
      if f(i) /= null then
        if not Certify_Factor(file,tol,f(i).all,grid) then
          put_line(file,"The factorization cannot be certified.");
          return false;
        end if;
      end if;
    end loop;
    put_line(file,"The factorization is certified.");
    return true;
  end Is_Factorization;

  function Is_Factorization
                ( tol : double_float;
                  f : Standard_Natural_VecVecs.VecVec;
                  grid : Array_of_DoblDobl_Sample_Lists ) return boolean is

    use Standard_Natural_Vectors;

  begin
    for i in f'range loop
      if f(i) /= null then
        if not Certify_Factor(tol,f(i).all,grid)
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Is_Factorization;

  function Is_Factorization
                ( file : file_type; tol : double_float;
                  f : Standard_Natural_VecVecs.VecVec;
                  grid : Array_of_DoblDobl_Sample_Lists ) return boolean is

    use Standard_Natural_Vectors;

  begin
    for i in f'range loop
      if f(i) /= null then
        if not Certify_Factor(file,tol,f(i).all,grid) then
          put_line(file,"The factorization cannot be certified.");
          return false;
        end if;
      end if;
    end loop;
    put_line(file,"The factorization is certified.");
    return true;
  end Is_Factorization;

  function Is_Factorization
                ( tol : double_float;
                  f : Standard_Natural_VecVecs.VecVec;
                  grid : Array_of_QuadDobl_Sample_Lists ) return boolean is

    use Standard_Natural_Vectors;

  begin
    for i in f'range loop
      if f(i) /= null then
        if not Certify_Factor(tol,f(i).all,grid)
         then return false;
        end if;
      end if;
    end loop;
    return true;
  end Is_Factorization;

  function Is_Factorization
                ( file : file_type; tol : double_float;
                  f : Standard_Natural_VecVecs.VecVec;
                  grid : Array_of_QuadDobl_Sample_Lists ) return boolean is

    use Standard_Natural_Vectors;

  begin
    for i in f'range loop
      if f(i) /= null then
        if not Certify_Factor(file,tol,f(i).all,grid) then
          put_line(file,"The factorization cannot be certified.");
          return false;
        end if;
      end if;
    end loop;
    put_line(file,"The factorization is certified.");
    return true;
  end Is_Factorization;

-- APPLICATION OF MONODROMY FOLLOWED BY LINEAR TRACES :

  procedure Monodromy_Breakup
                ( grid : in Array_of_Standard_Sample_Lists;
                  dim,threshold : in natural32; tol : in double_float;
                  f : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    sps : constant Standard_Sample_List := grid(0);
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    len : constant integer32 := integer32(Length_Of(sps));
    rvk : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    k0 : Standard_Complex_Vectors.Vector(1..len);
    nf : natural32 := natural32(len);
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : Standard_Complex_VecVecs.VecVec(1..integer32(dim));
    mapsps,mapsps_last,prevsps : Standard_Sample_List;
    sl,sl_last,tmp : Standard_Sample_Grid;
    cnt,nit : natural32 := 0;
    done : boolean;

  begin
    loop
      k0 := Keys(rvk,len,sps);
      newhyp := Random_Hyperplanes(dim,natural32(n)); 
      declare -- create new sample list
        newsps,newsps_last : Standard_Sample_List;
      begin
        Sample(sps,newhyp,newsps,newsps_last);
        if use_laurent
         then Sampling_Laurent_Machine.Change_Slices(newhyp);
         else Sampling_Machine.Change_Slices(newhyp);
        end if;
        Sample(newsps,hyp,mapsps,mapsps_last);
        Update(f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
        Deep_Clear(mapsps);
        done := Is_Factorization(tol,f.all,grid);
        exit when (done or (cnt >= threshold) or (nf = 1));
        tmp := sl;                                -- from new to old lists
        if use_laurent
         then Sampling_Laurent_Machine.Change_Slices(newhyp);
         else Sampling_Machine.Change_Slices(newhyp);
        end if;
        while not Is_Null(tmp) loop
          prevsps := Head_Of(tmp);
          k0 := Keys(rvk,len,prevsps);
          prevhyp := Hyperplane_Sections(Head_Of(prevsps));
          Sample(newsps,prevhyp,mapsps,mapsps_last);
          Update(f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
          Deep_Clear(mapsps);
          done := Is_Factorization(tol,f.all,grid);
          exit when (done or (cnt >= threshold) or (nf = 1));
          tmp := Tail_Of(tmp);
        end loop;
        exit when (done or (cnt >= threshold) or (nf = 1));
        Append(sl,sl_last,newsps);
      end;
      if use_laurent
       then Sampling_Laurent_Machine.Change_Slices(hyp);
       else Sampling_Machine.Change_Slices(hyp);
      end if;
    end loop;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( file : in file_type;
                  grid : in Array_of_Standard_Sample_Lists;
                  dim,threshold : in natural32; tol : in double_float;
                  f : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    sps : constant Standard_Sample_List := grid(0);
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    len : constant integer32 := integer32(Length_Of(sps));
    rvk : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    k0 : Standard_Complex_Vectors.Vector(1..len);
    nf : natural32 := natural32(len);
    hyp : constant Standard_Complex_VecVecs.VecVec(1..integer32(dim))
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : Standard_Complex_VecVecs.VecVec(1..integer32(dim));
    mapsps,mapsps_last,prevsps : Standard_Sample_List;
    sl,sl_last,tmp : Standard_Sample_Grid;
    cnt,nit : natural32 := 0;
    done : boolean;

  begin
    loop
      k0 := Keys(rvk,len,sps);
      newhyp := Random_Hyperplanes(dim,natural32(n)); 
      declare                                    -- create new sample list
        newsps,newsps_last : Standard_Sample_List;
      begin
        Sample(sps,newhyp,newsps,newsps_last);
        if use_laurent
         then Sampling_Laurent_Machine.Change_Slices(newhyp);
         else Sampling_Machine.Change_Slices(newhyp);
        end if;
        Sample(newsps,hyp,mapsps,mapsps_last);
        Update(file,f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
        Deep_Clear(mapsps);
        done := Is_Factorization(file,tol,f.all,grid);
        exit when (done or (cnt >= threshold) or (nf = 1));
        tmp := sl;                                -- from new to old lists
        if use_laurent
         then Sampling_Laurent_Machine.Change_Slices(newhyp);
         else Sampling_Machine.Change_Slices(newhyp);
        end if;
        while not Is_Null(tmp) loop
          prevsps := Head_Of(tmp);
          k0 := Keys(rvk,len,prevsps);
          prevhyp := Hyperplane_Sections(Head_Of(prevsps));
          Sample(newsps,prevhyp,mapsps,mapsps_last);
          Update(file,f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
          Deep_Clear(mapsps);
          done := Is_Factorization(file,tol,f.all,grid);
          exit when (done or (cnt >= threshold) or (nf = 1));
          tmp := Tail_Of(tmp);
        end loop;
        exit when (done or (cnt >= threshold) or (nf = 1));
        Append(sl,sl_last,newsps);
      end;
      if use_laurent
       then Sampling_Laurent_Machine.Change_Slices(hyp);
       else Sampling_Machine.Change_Slices(hyp);
      end if;
    end loop;
    put(file,"Number of monodromy loops : "); put(file,nit,1);
    new_line(file);
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( grid : in Array_of_DoblDobl_Sample_Lists;
                  dim,threshold : in natural32; tol : in double_float;
                  f : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    sps : constant DoblDobl_Sample_List := grid(0);
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    len : constant integer32 := integer32(Length_Of(sps));
    rvk : constant DoblDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    k0 : DoblDobl_Complex_Vectors.Vector(1..len);
    nf : natural32 := natural32(len);
    hyp : constant DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim));
    mapsps,mapsps_last,prevsps : DoblDobl_Sample_List;
    sl,sl_last,tmp : DoblDobl_Sample_Grid;
    cnt,nit : natural32 := 0;
    done : boolean;

  begin
    loop
      k0 := Keys(rvk,len,sps);
      newhyp := Random_Hyperplanes(dim,natural32(n)); 
      declare -- create new sample list
        newsps,newsps_last : DoblDobl_Sample_List;
      begin
        Sample(sps,newhyp,newsps,newsps_last);
        if use_laurent
         then DoblDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else DoblDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        Sample(newsps,hyp,mapsps,mapsps_last);
        Update(f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
        Deep_Clear(mapsps);
        done := Is_Factorization(tol,f.all,grid);
        exit when (done or (cnt >= threshold) or (nf = 1));
        tmp := sl;                                -- from new to old lists
        if use_laurent
         then DoblDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else DoblDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        while not Is_Null(tmp) loop
          prevsps := Head_Of(tmp);
          k0 := Keys(rvk,len,prevsps);
          prevhyp := Hyperplane_Sections(Head_Of(prevsps));
          Sample(newsps,prevhyp,mapsps,mapsps_last);
          Update(f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
          Deep_Clear(mapsps);
          done := Is_Factorization(tol,f.all,grid);
          exit when (done or (cnt >= threshold) or (nf = 1));
          tmp := Tail_Of(tmp);
        end loop;
        exit when (done or (cnt >= threshold) or (nf = 1));
        Append(sl,sl_last,newsps);
      end;
      if use_laurent
       then DoblDobl_Sampling_Laurent_Machine.Change_Slices(hyp);
       else DoblDobl_Sampling_Machine.Change_Slices(hyp);
      end if;
    end loop;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( file : in file_type;
                  grid : in Array_of_DoblDobl_Sample_Lists;
                  dim,threshold : in natural32; tol : in double_float;
                  f : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    sps : constant DoblDobl_Sample_List := grid(0);
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    len : constant integer32 := integer32(Length_Of(sps));
    rvk : constant DoblDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    k0 : DoblDobl_Complex_Vectors.Vector(1..len);
    nf : natural32 := natural32(len);
    hyp : constant DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : DoblDobl_Complex_VecVecs.VecVec(1..integer32(dim));
    mapsps,mapsps_last,prevsps : DoblDobl_Sample_List;
    sl,sl_last,tmp : DoblDobl_Sample_Grid;
    cnt,nit : natural32 := 0;
    done : boolean;

  begin
    loop
      k0 := Keys(rvk,len,sps);
      newhyp := Random_Hyperplanes(dim,natural32(n)); 
      declare                                    -- create new sample list
        newsps,newsps_last : DoblDobl_Sample_List;
      begin
        Sample(sps,newhyp,newsps,newsps_last);
        if use_laurent
         then DoblDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else DoblDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        Sample(newsps,hyp,mapsps,mapsps_last);
        Update(file,f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
        Deep_Clear(mapsps);
        done := Is_Factorization(file,tol,f.all,grid);
        exit when (done or (cnt >= threshold) or (nf = 1));
        tmp := sl;                                -- from new to old lists
        if use_laurent
         then DoblDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else DoblDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        while not Is_Null(tmp) loop
          prevsps := Head_Of(tmp);
          k0 := Keys(rvk,len,prevsps);
          prevhyp := Hyperplane_Sections(Head_Of(prevsps));
          Sample(newsps,prevhyp,mapsps,mapsps_last);
          Update(file,f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
          Deep_Clear(mapsps);
          done := Is_Factorization(file,tol,f.all,grid);
          exit when (done or (cnt >= threshold) or (nf = 1));
          tmp := Tail_Of(tmp);
        end loop;
        exit when (done or (cnt >= threshold) or (nf = 1));
        Append(sl,sl_last,newsps);
      end;
      if use_laurent
       then DoblDobl_Sampling_Laurent_Machine.Change_Slices(hyp);
       else DoblDobl_Sampling_Machine.Change_Slices(hyp);
      end if;
    end loop;
    put(file,"Number of monodromy loops : "); put(file,nit,1);
    new_line(file);
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( grid : in Array_of_QuadDobl_Sample_Lists;
                  dim,threshold : in natural32; tol : in double_float;
                  f : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    sps : constant QuadDobl_Sample_List := grid(0);
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    len : constant integer32 := integer32(Length_Of(sps));
    rvk : constant QuadDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    k0 : QuadDobl_Complex_Vectors.Vector(1..len);
    nf : natural32 := natural32(len);
    hyp : constant QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim));
    mapsps,mapsps_last,prevsps : QuadDobl_Sample_List;
    sl,sl_last,tmp : QuadDobl_Sample_Grid;
    cnt,nit : natural32 := 0;
    done : boolean;

  begin
    loop
      k0 := Keys(rvk,len,sps);
      newhyp := Random_Hyperplanes(dim,natural32(n)); 
      declare -- create new sample list
        newsps,newsps_last : QuadDobl_Sample_List;
      begin
        Sample(sps,newhyp,newsps,newsps_last);
        if use_laurent
         then QuadDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else QuadDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        Sample(newsps,hyp,mapsps,mapsps_last);
        Update(f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
        Deep_Clear(mapsps);
        done := Is_Factorization(tol,f.all,grid);
        exit when (done or (cnt >= threshold) or (nf = 1));
        tmp := sl;                                -- from new to old lists
        if use_laurent
         then QuadDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else QuadDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        while not Is_Null(tmp) loop
          prevsps := Head_Of(tmp);
          k0 := Keys(rvk,len,prevsps);
          prevhyp := Hyperplane_Sections(Head_Of(prevsps));
          Sample(newsps,prevhyp,mapsps,mapsps_last);
          Update(f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
          Deep_Clear(mapsps);
          done := Is_Factorization(tol,f.all,grid);
          exit when (done or (cnt >= threshold) or (nf = 1));
          tmp := Tail_Of(tmp);
        end loop;
        exit when (done or (cnt >= threshold) or (nf = 1));
        Append(sl,sl_last,newsps);
      end;
      if use_laurent
       then QuadDobl_Sampling_Laurent_Machine.Change_Slices(hyp);
       else QuadDobl_Sampling_Machine.Change_Slices(hyp);
      end if;
    end loop;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( file : in file_type;
                  grid : in Array_of_QuadDobl_Sample_Lists;
                  dim,threshold : in natural32; tol : in double_float;
                  f : in Standard_Natural_VecVecs.Link_to_VecVec ) is

    sps : constant QuadDobl_Sample_List := grid(0);
    n : constant integer32 := Number_of_Variables(Head_Of(sps));
    len : constant integer32 := integer32(Length_Of(sps));
    rvk : constant QuadDobl_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    k0 : QuadDobl_Complex_Vectors.Vector(1..len);
    nf : natural32 := natural32(len);
    hyp : constant QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim))
        := Hyperplane_Sections(Head_Of(sps));
    newhyp,prevhyp : QuadDobl_Complex_VecVecs.VecVec(1..integer32(dim));
    mapsps,mapsps_last,prevsps : QuadDobl_Sample_List;
    sl,sl_last,tmp : QuadDobl_Sample_Grid;
    cnt,nit : natural32 := 0;
    done : boolean;

  begin
    loop
      k0 := Keys(rvk,len,sps);
      newhyp := Random_Hyperplanes(dim,natural32(n)); 
      declare                                    -- create new sample list
        newsps,newsps_last : QuadDobl_Sample_List;
      begin
        Sample(sps,newhyp,newsps,newsps_last);
        if use_laurent
         then QuadDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else QuadDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        Sample(newsps,hyp,mapsps,mapsps_last);
        Update(file,f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
        Deep_Clear(mapsps);
        done := Is_Factorization(file,tol,f.all,grid);
        exit when (done or (cnt >= threshold) or (nf = 1));
        tmp := sl;                                -- from new to old lists
        if use_laurent
         then QuadDobl_Sampling_Laurent_Machine.Change_Slices(newhyp);
         else QuadDobl_Sampling_Machine.Change_Slices(newhyp);
        end if;
        while not Is_Null(tmp) loop
          prevsps := Head_Of(tmp);
          k0 := Keys(rvk,len,prevsps);
          prevhyp := Hyperplane_Sections(Head_Of(prevsps));
          Sample(newsps,prevhyp,mapsps,mapsps_last);
          Update(file,f,nf,cnt,nit,tol,k0,rvk,len,mapsps);
          Deep_Clear(mapsps);
          done := Is_Factorization(file,tol,f.all,grid);
          exit when (done or (cnt >= threshold) or (nf = 1));
          tmp := Tail_Of(tmp);
        end loop;
        exit when (done or (cnt >= threshold) or (nf = 1));
        Append(sl,sl_last,newsps);
      end;
      if use_laurent
       then QuadDobl_Sampling_Laurent_Machine.Change_Slices(hyp);
       else QuadDobl_Sampling_Machine.Change_Slices(hyp);
      end if;
    end loop;
    put(file,"Number of monodromy loops : "); put(file,nit,1);
    new_line(file);
  end Monodromy_Breakup;

-- DRIVER ROUTINES :

  procedure Factor ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     dim : in natural32;
                     grid : in Array_of_Standard_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := false;
    f := Init_Factors(d);
    Monodromy_Breakup(grid,dim,t,tol,f);
  end Factor;

  procedure Factor ( file : in file_type;
                     p : in Standard_Complex_Poly_Systems.Poly_Sys;
                     dim : in natural32;
                     grid : in Array_of_Standard_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    timer : Timing_Widget;
    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := false;
    new_line(file);
    put_line(file,"MONODROMY GROUP BREAKS UP INTO IRREDUCIBLE COMPONENTS");
    new_line(file);
    tstart(timer);
    f := Init_Factors(d);
    Monodromy_Breakup(file,grid,dim,t,tol,f);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Monodromy Factorization");
    new_line(file);
  end Factor;

  procedure Factor ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                     dim : in natural32;
                     grid : in Array_of_DoblDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := false;
    f := Init_Factors(d);
    Monodromy_Breakup(grid,dim,t,tol,f);
  end Factor;

  procedure Factor ( file : in file_type;
                     p : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                     dim : in natural32;
                     grid : in Array_of_DoblDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    timer : Timing_Widget;
    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := false;
    new_line(file);
    put_line(file,"MONODROMY GROUP BREAKS UP INTO IRREDUCIBLE COMPONENTS");
    new_line(file);
    tstart(timer);
    f := Init_Factors(d);
    Monodromy_Breakup(file,grid,dim,t,tol,f);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Monodromy Factorization");
    new_line(file);
  end Factor;

  procedure Factor ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                     dim : in natural32;
                     grid : in Array_of_QuadDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := false;
    f := Init_Factors(d);
    Monodromy_Breakup(grid,dim,t,tol,f);
  end Factor;

  procedure Factor ( file : in file_type;
                     p : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                     dim : in natural32;
                     grid : in Array_of_QuadDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    timer : Timing_Widget;
    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := false;
    new_line(file);
    put_line(file,"MONODROMY GROUP BREAKS UP INTO IRREDUCIBLE COMPONENTS");
    new_line(file);
    tstart(timer);
    f := Init_Factors(d);
    Monodromy_Breakup(file,grid,dim,t,tol,f);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Monodromy Factorization");
    new_line(file);
  end Factor;

  procedure Factor ( p : in Standard_Complex_Laur_Systems.Laur_Sys;
                     dim : in natural32;
                     grid : in Array_of_Standard_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := true;
    f := Init_Factors(d);
    Monodromy_Breakup(grid,dim,t,tol,f);
  end Factor;

  procedure Factor ( file : in file_type;
                     p : in Standard_Complex_Laur_Systems.Laur_Sys;
                     dim : in natural32;
                     grid : in Array_of_Standard_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    timer : Timing_Widget;
    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := true;
    new_line(file);
    put_line(file,"MONODROMY GROUP BREAKS UP INTO IRREDUCIBLE COMPONENTS");
    new_line(file);
    tstart(timer);
    f := Init_Factors(d);
    Monodromy_Breakup(file,grid,dim,t,tol,f);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Monodromy Factorization");
    new_line(file);
  end Factor;

  procedure Factor ( p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                     dim : in natural32;
                     grid : in Array_of_DoblDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := true;
    f := Init_Factors(d);
    Monodromy_Breakup(grid,dim,t,tol,f);
  end Factor;

  procedure Factor ( file : in file_type;
                     p : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                     dim : in natural32;
                     grid : in Array_of_DoblDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    timer : Timing_Widget;
    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := true;
    new_line(file);
    put_line(file,"MONODROMY GROUP BREAKS UP INTO IRREDUCIBLE COMPONENTS");
    new_line(file);
    tstart(timer);
    f := Init_Factors(d);
    Monodromy_Breakup(file,grid,dim,t,tol,f);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Monodromy Factorization");
    new_line(file);
  end Factor;

  procedure Factor ( p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                     dim : in natural32;
                     grid : in Array_of_QuadDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := true;
    f := Init_Factors(d);
    Monodromy_Breakup(grid,dim,t,tol,f);
  end Factor;

  procedure Factor ( file : in file_type;
                     p : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                     dim : in natural32;
                     grid : in Array_of_QuadDobl_Sample_Lists;
                     f : out Standard_Natural_VecVecs.Link_to_VecVec ) is

    timer : Timing_Widget;
    d : constant natural32 := Length_Of(grid(grid'first));
    t : constant natural32 := 10;
    tol : constant double_float := 1.0E-8;

  begin
    use_laurent := true;
    new_line(file);
    put_line(file,"MONODROMY GROUP BREAKS UP INTO IRREDUCIBLE COMPONENTS");
    new_line(file);
    tstart(timer);
    f := Init_Factors(d);
    Monodromy_Breakup(file,grid,dim,t,tol,f);
    tstop(timer);
    new_line(file);
    print_times(file,timer,"Monodromy Factorization");
    new_line(file);
  end Factor;

begin
  use_laurent := false;
end Monodromy_Component_Breakup;
