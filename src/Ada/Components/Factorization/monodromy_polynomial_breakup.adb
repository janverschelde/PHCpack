with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Natural_Vectors;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with DoblDobl_Random_Vectors;            use DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;            use QuadDobl_Random_Vectors;
with Standard_Complex_VecLists;
with DoblDobl_Complex_VecLists;
with QuadDobl_Complex_VecLists;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_SysFun;
with Standard_Lined_Hypersurfaces;       use Standard_Lined_Hypersurfaces;
with DoblDobl_Lined_Hypersurfaces;       use DoblDobl_Lined_Hypersurfaces;
with QuadDobl_Lined_Hypersurfaces;       use QuadDobl_Lined_Hypersurfaces;
with Monodromy_Partitions;               use Monodromy_Partitions;

package body Monodromy_Polynomial_Breakup is

  tol : constant double_float := 1.0E-8;
  threshold : constant natural32 := 10;

-- AUXILIARY OPERATIONS :

  function Good_Map ( mp : Standard_Natural_Vectors.Vector )
                    return boolean is

  -- DESCRIPTION :
  --   A good map contains no zero elements.

  begin
    for i in mp'range loop
      if mp(i) = 0
       then return false;
      end if;
    end loop;
    return true;
  end Good_Map;

  function Eval_Sys ( n : integer32;
                      p : Standard_Complex_Polynomials.Poly;
                      ep : Standard_Complex_Poly_Functions.Eval_Poly )
                    return Standard_Complex_Poly_SysFun.Eval_Poly_Sys is

  -- DESCRIPTION :
  --   The array on return has range 0..n, where n is the number
  --   of variables in the polynomial p.  The 0-th entry is ep,
  --   while the i-th entry is the i-th derivative of p.
  --   This version is for standard double precision.

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Functions;
    use Standard_Complex_Poly_SysFun;

    res : Eval_Poly_Sys(0..n);
    dp : Poly;

  begin
    res(0) := ep;
    for i in 1..n loop
      dp := Diff(p,i);
      res(i) := Create(dp);
      Clear(dp);
    end loop;
    return res;
  end Eval_Sys;

  function Eval_Sys ( n : integer32;
                      p : DoblDobl_Complex_Polynomials.Poly;
                      ep : DoblDobl_Complex_Poly_Functions.Eval_Poly )
                    return DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys is

  -- DESCRIPTION :
  --   The array on return has range 0..n, where n is the number
  --   of variables in the polynomial p.  The 0-th entry is ep,
  --   while the i-th entry is the i-th derivative of p.
  --   This version is for double double precision.

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Functions;
    use DoblDobl_Complex_Poly_SysFun;

    res : Eval_Poly_Sys(0..n);
    dp : Poly;

  begin
    res(0) := ep;
    for i in 1..n loop
      dp := Diff(p,i);
      res(i) := Create(dp);
      Clear(dp);
    end loop;
    return res;
  end Eval_Sys;

  function Eval_Sys ( n : integer32;
                      p : QuadDobl_Complex_Polynomials.Poly;
                      ep : QuadDobl_Complex_Poly_Functions.Eval_Poly )
                    return QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys is

  -- DESCRIPTION :
  --   The array on return has range 0..n, where n is the number
  --   of variables in the polynomial p.  The 0-th entry is ep,
  --   while the i-th entry is the i-th derivative of p.
  --   This version is for quad double precision.

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Functions;
    use QuadDobl_Complex_Poly_SysFun;

    res : Eval_Poly_Sys(0..n);
    dp : Poly;

  begin
    res(0) := ep;
    for i in 1..n loop
      dp := Diff(p,i);
      res(i) := Create(dp);
      Clear(dp);
    end loop;
    return res;
  end Eval_Sys;

  procedure Append_Line 
              ( bl,bl_last,vl,vl_last,tl,tl_last
                  : in out Standard_Complex_VecLists.List;
                b,v,t : in Standard_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   A line consists of a base point b, direction v and roots t.

  -- ON ENTRY :
  --   (bl,bl_last)    first and last element of list of base points;
  --   (vl,vl_last)    first and last element of list of directions;
  --   (tl,tl_last)    first and last element of list of roots;
  --   (b,v,t)         new line.

  -- ON RETURN :
  --   (bl,bl_last)    updated list of base points;
  --   (vl,vl_last)    updated list of directions;
  --   (tl,tl_last)    updated list of roots.

    use Standard_Complex_VecLists;

  begin
    Append(bl,bl_last,b);
    Append(vl,vl_last,v);
    Append(tl,tl_last,t);
  end Append_Line;

  procedure Append_Line 
              ( bl,bl_last,vl,vl_last,tl,tl_last
                  : in out DoblDobl_Complex_VecLists.List;
                b,v,t : in DoblDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   A line consists of a base point b, direction v and roots t.

  -- ON ENTRY :
  --   (bl,bl_last)    first and last element of list of base points;
  --   (vl,vl_last)    first and last element of list of directions;
  --   (tl,tl_last)    first and last element of list of roots;
  --   (b,v,t)         new line.

  -- ON RETURN :
  --   (bl,bl_last)    updated list of base points;
  --   (vl,vl_last)    updated list of directions;
  --   (tl,tl_last)    updated list of roots.

    use DoblDobl_Complex_VecLists;

  begin
    Append(bl,bl_last,b);
    Append(vl,vl_last,v);
    Append(tl,tl_last,t);
  end Append_Line;

  procedure Append_Line 
              ( bl,bl_last,vl,vl_last,tl,tl_last
                  : in out QuadDobl_Complex_VecLists.List;
                b,v,t : in QuadDobl_Complex_Vectors.Vector ) is

  -- DESCRIPTION :
  --   A line consists of a base point b, direction v and roots t.

  -- ON ENTRY :
  --   (bl,bl_last)    first and last element of list of base points;
  --   (vl,vl_last)    first and last element of list of directions;
  --   (tl,tl_last)    first and last element of list of roots;
  --   (b,v,t)         new line.

  -- ON RETURN :
  --   (bl,bl_last)    updated list of base points;
  --   (vl,vl_last)    updated list of directions;
  --   (tl,tl_last)    updated list of roots.

    use QuadDobl_Complex_VecLists;

  begin
    Append(bl,bl_last,b);
    Append(vl,vl_last,v);
    Append(tl,tl_last,t);
  end Append_Line;

  procedure First_Loop
              ( fp : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out Standard_Complex_VecLists.List;
                b0,v0,t0 : in Standard_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb : in out natural32;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Executes the first monodromy loop, from the line (b0,v0) with
  --   the roots in t0.  Initializes the lists (see the Append_Line
  --   above for their meaning) with three lines and their roots.
  --   On return is the map connecting points on the same factor.
  --   If fail, then the map was found to contain zeros.

    use Standard_Complex_Vectors;
    use Standard_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b1 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v1 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    b2 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v2 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    t1,t2 : Vector(t0'range);
    ft,dt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);

  begin
   -- b1(1) := Create(0.0);
   -- b2(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v1(i) := Create(0.0);
   --   v2(i) := Create(0.0);
   -- end loop;
    Silent_Hypersurface_Sampler(fp,b0,v0,b1,v1,t);
    Silent_Refiner(fp,b1,v1,t,ft,dt,eps,max);
    t1 := t;
    Silent_Hypersurface_Sampler(fp,b1,v1,b2,v2,t);
    Silent_Refiner(fp,b2,v2,t,ft,dt,eps,max);
    t2 := t;
    Silent_Hypersurface_Sampler(fp,b2,v2,b0,v0,t);
    Silent_Refiner(fp,b0,v0,t,ft,dt,eps,max);
    mp := Map(t0,t,tol);
    fail := not Good_Map(mp);
    if not fail then
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b0,v0,t0);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b1,v1,t1);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b2,v2,t2);
      nb := natural32(deco'last);
      Add_Map(deco,nb,mp);
    end if;
  end First_Loop;

  procedure First_Loop
              ( fp : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out DoblDobl_Complex_VecLists.List;
                b0,v0,t0 : in DoblDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb : in out natural32;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Executes the first monodromy loop, from the line (b0,v0) with
  --   the roots in t0.  Initializes the lists (see the Append_Line
  --   above for their meaning) with three lines and their roots.
  --   On return is the map connecting points on the same factor.
  --   If fail, then the map was found to contain zeros.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b1 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v1 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    b2 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v2 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    t1,t2 : Vector(t0'range);
    ft,dt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);

  begin
   -- b1(1) := Create(0.0);
   -- b2(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v1(i) := Create(0.0);
   --   v2(i) := Create(0.0);
   -- end loop;
    Silent_Hypersurface_Sampler(fp,b0,v0,b1,v1,t);
    Silent_Refiner(fp,b1,v1,t,ft,dt,eps,max);
    t1 := t;
    Silent_Hypersurface_Sampler(fp,b1,v1,b2,v2,t);
    Silent_Refiner(fp,b2,v2,t,ft,dt,eps,max);
    t2 := t;
    Silent_Hypersurface_Sampler(fp,b2,v2,b0,v0,t);
    Silent_Refiner(fp,b0,v0,t,ft,dt,eps,max);
    mp := Map(t0,t,tol);
    fail := not Good_Map(mp);
    if not fail then
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b0,v0,t0);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b1,v1,t1);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b2,v2,t2);
      nb := natural32(deco'last);
      Add_Map(deco,nb,mp);
    end if;
  end First_Loop;

  procedure First_Loop
              ( fp : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out QuadDobl_Complex_VecLists.List;
                b0,v0,t0 : in QuadDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb : in out natural32;
                fail : out boolean ) is

  -- DESCRIPTION :
  --   Executes the first monodromy loop, from the line (b0,v0) with
  --   the roots in t0.  Initializes the lists (see the Append_Line
  --   above for their meaning) with three lines and their roots.
  --   On return is the map connecting points on the same factor.
  --   If fail, then the map was found to contain zeros.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b1 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v1 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    b2 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v2 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    t1,t2 : Vector(t0'range);
    ft,dt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);

  begin
   -- b1(1) := Create(0.0);
   -- b2(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v1(i) := Create(0.0);
   --   v2(i) := Create(0.0);
   -- end loop;
    Silent_Hypersurface_Sampler(fp,b0,v0,b1,v1,t);
    Silent_Refiner(fp,b1,v1,t,ft,dt,eps,max);
    t1 := t;
    Silent_Hypersurface_Sampler(fp,b1,v1,b2,v2,t);
    Silent_Refiner(fp,b2,v2,t,ft,dt,eps,max);
    t2 := t;
    Silent_Hypersurface_Sampler(fp,b2,v2,b0,v0,t);
    Silent_Refiner(fp,b0,v0,t,ft,dt,eps,max);
    mp := Map(t0,t,tol);
    fail := not Good_Map(mp);
    if not fail then
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b0,v0,t0);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b1,v1,t1);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b2,v2,t2);
      nb := natural32(deco'last);
      Add_Map(deco,nb,mp);
    end if;
  end First_Loop;

  procedure First_Loop
              ( file : in file_type;
                fp : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                output : in boolean;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out Standard_Complex_VecLists.List;
                b0,v0,t0 : in Standard_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Executes the first monodromy loop, from the line (b0,v0) with
  --   the roots in t0.  Initializes the lists (see the Append_Line
  --   above for their meaning) with three lines and their roots.
  --   On return is the updated connectivity information.

    use Standard_Complex_Vectors;
    use Standard_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b1 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v1 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    b2 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v2 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    t1,t2 : Vector(t0'range);
    ft,dt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);

  begin
   -- b1(1) := Create(0.0);
   -- b2(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v1(i) := Create(0.0);
   --   v2(i) := Create(0.0);
   -- end loop;
    Reporting_Hypersurface_Sampler(file,fp,b0,v0,b1,v1,output,t);
   -- Reporting_Refiner(file,fp,b1,v1,t,ft,dt,eps,max);
    Silent_Refiner(fp,b1,v1,t,ft,dt,eps,max);
    t1 := t;
    Reporting_Hypersurface_Sampler(file,fp,b1,v1,b2,v2,output,t);
   -- Reporting_Refiner(file,fp,b2,v2,t,ft,dt,eps,max);
    Silent_Refiner(fp,b2,v2,t,ft,dt,eps,max);
    t2 := t;
    Reporting_Hypersurface_Sampler(file,fp,b2,v2,b0,v0,output,t);
   -- Reporting_Refiner(file,fp,b0,v0,t,ft,dt,eps,max);
    Silent_Refiner(fp,b0,v0,t,ft,dt,eps,max);
    mp := Map(t0,t,tol);
    Write_Map(file,mp);
    fail := not Good_Map(mp);
    if fail then
      put_line(file,"The map is not good for processing.");
    else
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b0,v0,t0);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b1,v1,t1);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b2,v2,t2);
      nb := natural32(deco'last);
      Add_Map(deco,nb,mp);
    end if;
  end First_Loop;

  procedure First_Loop
              ( file : in file_type;
                fp : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                output : in boolean;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out DoblDobl_Complex_VecLists.List;
                b0,v0,t0 : in DoblDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Executes the first monodromy loop, from the line (b0,v0) with
  --   the roots in t0.  Initializes the lists (see the Append_Line
  --   above for their meaning) with three lines and their roots.
  --   On return is the updated connectivity information.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b1 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v1 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    b2 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v2 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    t1,t2 : Vector(t0'range);
    ft,dt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);

  begin
   -- b1(1) := Create(0.0);
   -- b2(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v1(i) := Create(0.0);
   --   v2(i) := Create(0.0);
   -- end loop;
    Reporting_Hypersurface_Sampler(file,fp,b0,v0,b1,v1,output,t);
   -- Reporting_Refiner(file,fp,b1,v1,t,ft,dt,eps,max);
    Silent_Refiner(fp,b1,v1,t,ft,dt,eps,max);
    t1 := t;
    Reporting_Hypersurface_Sampler(file,fp,b1,v1,b2,v2,output,t);
   -- Reporting_Refiner(file,fp,b2,v2,t,ft,dt,eps,max);
    Silent_Refiner(fp,b2,v2,t,ft,dt,eps,max);
    t2 := t;
    Reporting_Hypersurface_Sampler(file,fp,b2,v2,b0,v0,output,t);
   -- Reporting_Refiner(file,fp,b0,v0,t,ft,dt,eps,max);
    Silent_Refiner(fp,b0,v0,t,ft,dt,eps,max);
    mp := Map(t0,t,tol);
    Write_Map(file,mp);
    fail := not Good_Map(mp);
    if fail then
      put_line(file,"The map is not good for processing.");
    else
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b0,v0,t0);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b1,v1,t1);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b2,v2,t2);
      nb := natural32(deco'last);
      Add_Map(deco,nb,mp);
    end if;
  end First_Loop;

  procedure First_Loop
              ( file : in file_type;
                fp : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                output : in boolean;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out QuadDobl_Complex_VecLists.List;
                b0,v0,t0 : in QuadDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Executes the first monodromy loop, from the line (b0,v0) with
  --   the roots in t0.  Initializes the lists (see the Append_Line
  --   above for their meaning) with three lines and their roots.
  --   On return is the updated connectivity information.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b1 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v1 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    b2 : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v2 : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    t1,t2 : Vector(t0'range);
    ft,dt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);

  begin
   -- b1(1) := Create(0.0);
   -- b2(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v1(i) := Create(0.0);
   --   v2(i) := Create(0.0);
   -- end loop;
    Reporting_Hypersurface_Sampler(file,fp,b0,v0,b1,v1,output,t);
   -- Reporting_Refiner(file,fp,b1,v1,t,ft,dt,eps,max);
    Silent_Refiner(fp,b1,v1,t,ft,dt,eps,max);
    t1 := t;
    Reporting_Hypersurface_Sampler(file,fp,b1,v1,b2,v2,output,t);
   -- Reporting_Refiner(file,fp,b2,v2,t,ft,dt,eps,max);
    Silent_Refiner(fp,b2,v2,t,ft,dt,eps,max);
    t2 := t;
    Reporting_Hypersurface_Sampler(file,fp,b2,v2,b0,v0,output,t);
   -- Reporting_Refiner(file,fp,b0,v0,t,ft,dt,eps,max);
    Silent_Refiner(fp,b0,v0,t,ft,dt,eps,max);
    mp := Map(t0,t,tol);
    Write_Map(file,mp);
    fail := not Good_Map(mp);
    if fail then
      put_line(file,"The map is not good for processing.");
    else
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b0,v0,t0);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b1,v1,t1);
      Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b2,v2,t2);
      nb := natural32(deco'last);
      Add_Map(deco,nb,mp);
    end if;
  end First_Loop;

  procedure New_Loop
              ( fp : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                d : in natural32;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out Standard_Complex_VecLists.List;
                b0,v0,t0 : in Standard_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb,scnt : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Attempts to find a new monodromy group action.

    use Standard_Complex_Vectors;
    use Standard_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    wb,wv : Vector(b'range);
    ft,dt,wt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);
    tmpbl : List := Tail_Of(bl);
    tmpvl : List := Tail_Of(vl);
    tmptl : List := Tail_Of(tl);
    prvnb : natural32 := nb;
    gm : boolean;

  begin
   -- b(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v(i) := Create(0.0);
   -- end loop;
    fail := true;
    Silent_Hypersurface_Sampler(fp,b0,v0,b,v,t);
    Silent_Refiner(fp,b,v,t,ft,dt,eps,max);
    while not Is_Null(tmpbl) loop
      wb := Head_Of(tmpbl).all;
      wv := Head_Of(tmpvl).all;
      wt := t;
      Silent_Hypersurface_Sampler(fp,b,v,wb,wv,wt);
      Silent_Refiner(fp,wb,wv,wt,ft,dt,eps,max);
      mp := Map(Head_Of(tmptl).all,wt,tol);
      gm := Good_Map(mp);
      if gm then
        Add_Map(deco,nb,mp);
        fail := false;
      else
        return;
      end if;
      exit when (nb = 1);
      if gm then
        if prvnb = nb then
          scnt := scnt + 1;
        else
          prvnb := nb;
          scnt := 0;
        end if;
      end if;
      exit when (scnt >= threshold);
      tmpbl := Tail_Of(tmpbl);
      tmpvl := Tail_Of(tmpvl);
      tmptl := Tail_Of(tmptl);
    end loop;
    if not fail
     then Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b,v,t);
    end if;
  end New_Loop;

  procedure New_Loop
              ( fp : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                d : in natural32;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out DoblDobl_Complex_VecLists.List;
                b0,v0,t0 : in DoblDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb,scnt : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Attempts to find a new monodromy group action.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    wb,wv : Vector(b'range);
    ft,dt,wt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);
    tmpbl : List := Tail_Of(bl);
    tmpvl : List := Tail_Of(vl);
    tmptl : List := Tail_Of(tl);
    prvnb : natural32 := nb;
    gm : boolean;

  begin
   -- b(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v(i) := Create(0.0);
   -- end loop;
    fail := true;
    Silent_Hypersurface_Sampler(fp,b0,v0,b,v,t);
    Silent_Refiner(fp,b,v,t,ft,dt,eps,max);
    while not Is_Null(tmpbl) loop
      wb := Head_Of(tmpbl).all;
      wv := Head_Of(tmpvl).all;
      wt := t;
      Silent_Hypersurface_Sampler(fp,b,v,wb,wv,wt);
      Silent_Refiner(fp,wb,wv,wt,ft,dt,eps,max);
      mp := Map(Head_Of(tmptl).all,wt,tol);
      gm := Good_Map(mp);
      if gm then
        Add_Map(deco,nb,mp);
        fail := false;
      else
        return;
      end if;
      exit when (nb = 1);
      if gm then
        if prvnb = nb then
          scnt := scnt + 1;
        else
          prvnb := nb;
          scnt := 0;
        end if;
      end if;
      exit when (scnt >= threshold);
      tmpbl := Tail_Of(tmpbl);
      tmpvl := Tail_Of(tmpvl);
      tmptl := Tail_Of(tmptl);
    end loop;
    if not fail
     then Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b,v,t);
    end if;
  end New_Loop;

  procedure New_Loop
              ( fp : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                d : in natural32;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out QuadDobl_Complex_VecLists.List;
                b0,v0,t0 : in QuadDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb,scnt : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Attempts to find a new monodromy group action.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-13;
    max : constant natural32 := 4;
    b : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    wb,wv : Vector(b'range);
    ft,dt,wt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);
    tmpbl : List := Tail_Of(bl);
    tmpvl : List := Tail_Of(vl);
    tmptl : List := Tail_Of(tl);
    prvnb : natural32 := nb;
    gm : boolean;

  begin
   -- b(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v(i) := Create(0.0);
   -- end loop;
    fail := true;
    Silent_Hypersurface_Sampler(fp,b0,v0,b,v,t);
    Silent_Refiner(fp,b,v,t,ft,dt,eps,max);
    while not Is_Null(tmpbl) loop
      wb := Head_Of(tmpbl).all;
      wv := Head_Of(tmpvl).all;
      wt := t;
      Silent_Hypersurface_Sampler(fp,b,v,wb,wv,wt);
      Silent_Refiner(fp,wb,wv,wt,ft,dt,eps,max);
      mp := Map(Head_Of(tmptl).all,wt,tol);
      gm := Good_Map(mp);
      if gm then
        Add_Map(deco,nb,mp);
        fail := false;
      else
        return;
      end if;
      exit when (nb = 1);
      if gm then
        if prvnb = nb then
          scnt := scnt + 1;
        else
          prvnb := nb;
          scnt := 0;
        end if;
      end if;
      exit when (scnt >= threshold);
      tmpbl := Tail_Of(tmpbl);
      tmpvl := Tail_Of(tmpvl);
      tmptl := Tail_Of(tmptl);
    end loop;
    if not fail
     then Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b,v,t);
    end if;
  end New_Loop;

  procedure New_Loop
              ( file : in file_type; d : in natural32;
                fp : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                output : in boolean;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out Standard_Complex_VecLists.List;
                b0,v0,t0 : in Standard_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb,scnt : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Attempts to find a new monodromy group action.

    use Standard_Complex_Vectors;
    use Standard_Complex_VecLists;

    eps : constant double_float := 1.0E-14;
    max : constant natural32 := 4;
    b : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    wb,wv : Vector(b'range);
    ft,dt,wt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);
    tmpbl : List := Tail_Of(bl);
    tmpvl : List := Tail_Of(vl);
    tmptl : List := Tail_Of(tl);
    prvnb : natural32 := nb;
    gm : boolean;

  begin
   -- b(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v(i) := Create(0.0);
   -- end loop;
    Reporting_Hypersurface_Sampler(file,fp,b0,v0,b,v,output,t);
   -- Reporting_Refiner(file,fp,b,v,t,ft,dt,eps,max);
    Silent_Refiner(fp,b,v,t,ft,dt,eps,max);
    fail := true;
    while not Is_Null(tmpbl) loop
      wb := Head_Of(tmpbl).all;
      wv := Head_Of(tmpvl).all;
      wt := t;
      Reporting_Hypersurface_Sampler(file,fp,b,v,wb,wv,output,wt);
     -- Reporting_Refiner(file,fp,wb,wv,wt,ft,dt,eps,max);
      Silent_Refiner(fp,wb,wv,wt,ft,dt,eps,max);
      mp := Map(Head_Of(tmptl).all,wt,tol);
      Write_Map(file,mp);
      gm := Good_Map(mp);
      if gm then
        Add_Map(deco,nb,mp);
        fail := false;
      else
        return;
      end if;
      exit when (nb = 1);
      if gm then
        if prvnb = nb then
          scnt := scnt + 1;
        else
          prvnb := nb;
          scnt := 0;
        end if;
      end if;
      exit when (scnt > threshold);
      tmpbl := Tail_Of(tmpbl);
      tmpvl := Tail_Of(tmpvl);
      tmptl := Tail_Of(tmptl);
    end loop;
    if not fail
     then Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b,v,t);
    end if;
  end New_Loop;

  procedure New_Loop
              ( file : in file_type; d : in natural32;
                fp : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                output : in boolean;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out DoblDobl_Complex_VecLists.List;
                b0,v0,t0 : in DoblDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb,scnt : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Attempts to find a new monodromy group action.

    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-14;
    max : constant natural32 := 4;
    b : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    wb,wv : Vector(b'range);
    ft,dt,wt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);
    tmpbl : List := Tail_Of(bl);
    tmpvl : List := Tail_Of(vl);
    tmptl : List := Tail_Of(tl);
    prvnb : natural32 := nb;
    gm : boolean;

  begin
   -- b(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v(i) := Create(0.0);
   -- end loop;
    Reporting_Hypersurface_Sampler(file,fp,b0,v0,b,v,output,t);
   -- Reporting_Refiner(file,fp,b,v,t,ft,dt,eps,max);
    Silent_Refiner(fp,b,v,t,ft,dt,eps,max);
    fail := true;
    while not Is_Null(tmpbl) loop
      wb := Head_Of(tmpbl).all;
      wv := Head_Of(tmpvl).all;
      wt := t;
      Reporting_Hypersurface_Sampler(file,fp,b,v,wb,wv,output,wt);
     -- Reporting_Refiner(file,fp,wb,wv,wt,ft,dt,eps,max);
      Silent_Refiner(fp,wb,wv,wt,ft,dt,eps,max);
      mp := Map(Head_Of(tmptl).all,wt,tol);
      Write_Map(file,mp);
      gm := Good_Map(mp);
      if gm then
        Add_Map(deco,nb,mp);
        fail := false;
      else
        return;
      end if;
      exit when (nb = 1);
      if gm then
        if prvnb = nb then
          scnt := scnt + 1;
        else
          prvnb := nb;
          scnt := 0;
        end if;
      end if;
      exit when (scnt > threshold);
      tmpbl := Tail_Of(tmpbl);
      tmpvl := Tail_Of(tmpvl);
      tmptl := Tail_Of(tmptl);
    end loop;
    if not fail
     then Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b,v,t);
    end if;
  end New_Loop;

  procedure New_Loop
              ( file : in file_type; d : in natural32;
                fp : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                output : in boolean;
                bl,bl_last,vl,vl_last,tl,tl_last
                  : in out QuadDobl_Complex_VecLists.List;
                b0,v0,t0 : in QuadDobl_Complex_Vectors.Vector;
                deco : in out Standard_Natural_VecVecs.Link_to_VecVec;
                nb,scnt : in out natural32; fail : out boolean ) is

  -- DESCRIPTION :
  --   Attempts to find a new monodromy group action.

    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_VecLists;

    eps : constant double_float := 1.0E-14;
    max : constant natural32 := 4;
    b : constant Vector(b0'range) := Random_Vector(b0'first,b0'last);
    v : constant Vector(v0'range) := Random_Vector(v0'first,v0'last);
    t : Vector(t0'range) := t0;
    wb,wv : Vector(b'range);
    ft,dt,wt : Vector(t'range);
    mp : Standard_Natural_Vectors.Vector(t'range);
    tmpbl : List := Tail_Of(bl);
    tmpvl : List := Tail_Of(vl);
    tmptl : List := Tail_Of(tl);
    prvnb : natural32 := nb;
    gm : boolean;

  begin
   -- b(1) := Create(0.0);
   -- for i in v0'first+1..v0'last loop
   --   v(i) := Create(0.0);
   -- end loop;
    Reporting_Hypersurface_Sampler(file,fp,b0,v0,b,v,output,t);
   -- Reporting_Refiner(file,fp,b,v,t,ft,dt,eps,max);
    Silent_Refiner(fp,b,v,t,ft,dt,eps,max);
    fail := true;
    while not Is_Null(tmpbl) loop
      wb := Head_Of(tmpbl).all;
      wv := Head_Of(tmpvl).all;
      wt := t;
      Reporting_Hypersurface_Sampler(file,fp,b,v,wb,wv,output,wt);
     -- Reporting_Refiner(file,fp,wb,wv,wt,ft,dt,eps,max);
      Silent_Refiner(fp,wb,wv,wt,ft,dt,eps,max);
      mp := Map(Head_Of(tmptl).all,wt,tol);
      Write_Map(file,mp);
      gm := Good_Map(mp);
      if gm then
        Add_Map(deco,nb,mp);
        fail := false;
      else
        return;
      end if;
      exit when (nb = 1);
      if gm then
        if prvnb = nb then
          scnt := scnt + 1;
        else
          prvnb := nb;
          scnt := 0;
        end if;
      end if;
      exit when (scnt > threshold);
      tmpbl := Tail_Of(tmpbl);
      tmpvl := Tail_Of(tmpvl);
      tmptl := Tail_Of(tmptl);
    end loop;
    if not fail
     then Append_Line(bl,bl_last,vl,vl_last,tl,tl_last,b,v,t);
    end if;
  end New_Loop;

-- TARGET PROCEDURES :

  procedure Monodromy_Breakup
                ( file : in file_type; 
                  p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in Standard_Complex_Vectors.Vector;
                  output : in boolean;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

    use Standard_Natural_VecVecs;
    use Standard_Complex_VecLists;
    use Standard_Complex_Poly_SysFun;

    cnt,nb : natural32 := 0;
    d : constant natural32 := natural32(gp'length); 
    fp : constant Eval_Poly_Sys(0..b'last) := Eval_Sys(b'last,p,ep);
    bl,bl_last,vl,vl_last,tl,tl_last : List;
    fail : boolean := true;

  begin
    while fail loop
      First_Loop(file,fp,output,bl,bl_last,vl,vl_last,tl,tl_last,b,v,
                 gp,deco,nb,fail);
    end loop;
    if deco /= null
     then nb := Number_of_Factors(deco.all);
     else nb := 0;
    end if;
    if nb > 1 then
      loop
        New_Loop(file,d,fp,output,bl,bl_last,vl,vl_last,tl,tl_last,
                 b,v,gp,deco,nb,cnt,fail);
        nb := Number_of_Factors(deco.all);
        exit when (nb = 1) or (cnt >= threshold);
      end loop;
    end if;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( file : in file_type; 
                  p : in DoblDobl_Complex_Polynomials.Poly;
                  ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in DoblDobl_Complex_Vectors.Vector;
                  output : in boolean;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

    use Standard_Natural_VecVecs;
    use DoblDobl_Complex_VecLists;
    use DoblDobl_Complex_Poly_SysFun;

    cnt,nb : natural32 := 0;
    d : constant natural32 := natural32(gp'length); 
    fp : constant Eval_Poly_Sys(0..b'last) := Eval_Sys(b'last,p,ep);
    bl,bl_last,vl,vl_last,tl,tl_last : List;
    fail : boolean := true;

  begin
    while fail loop
      First_Loop(file,fp,output,bl,bl_last,vl,vl_last,tl,tl_last,b,v,
                 gp,deco,nb,fail);
    end loop;
    if deco /= null
     then nb := Number_of_Factors(deco.all);
     else nb := 0;
    end if;
    if nb > 1 then
      loop
        New_Loop(file,d,fp,output,bl,bl_last,vl,vl_last,tl,tl_last,
                 b,v,gp,deco,nb,cnt,fail);
        nb := Number_of_Factors(deco.all);
        exit when (nb = 1) or (cnt >= threshold);
      end loop;
    end if;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( file : in file_type; 
                  p : in QuadDobl_Complex_Polynomials.Poly;
                  ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in QuadDobl_Complex_Vectors.Vector;
                  output : in boolean;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

    use Standard_Natural_VecVecs;
    use QuadDobl_Complex_VecLists;
    use QuadDobl_Complex_Poly_SysFun;

    cnt,nb : natural32 := 0;
    d : constant natural32 := natural32(gp'length); 
    fp : constant Eval_Poly_Sys(0..b'last) := Eval_Sys(b'last,p,ep);
    bl,bl_last,vl,vl_last,tl,tl_last : List;
    fail : boolean := true;

  begin
    while fail loop
      First_Loop(file,fp,output,bl,bl_last,vl,vl_last,tl,tl_last,b,v,
                 gp,deco,nb,fail);
    end loop;
    if deco /= null
     then nb := Number_of_Factors(deco.all);
     else nb := 0;
    end if;
    if nb > 1 then
      loop
        New_Loop(file,d,fp,output,bl,bl_last,vl,vl_last,tl,tl_last,
                 b,v,gp,deco,nb,cnt,fail);
        nb := Number_of_Factors(deco.all);
        exit when (nb = 1) or (cnt >= threshold);
      end loop;
    end if;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( p : in Standard_Complex_Polynomials.Poly;
                  ep : in Standard_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in Standard_Complex_Vectors.Vector;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

    use Standard_Natural_VecVecs;
    use Standard_Complex_VecLists;
    use Standard_Complex_Poly_SysFun;

    cnt,nb : natural32 := 0;
    d : constant natural32 := natural32(gp'length);
    fp : constant Eval_Poly_Sys(0..b'last) := Eval_Sys(b'last,p,ep);
    bl,bl_last,vl,vl_last,tl,tl_last : List;
    fail : boolean := true;

  begin
    while fail loop
      First_Loop(fp,bl,bl_last,vl,vl_last,tl,tl_last,b,v,gp,deco,nb,fail);
    end loop;
    if deco /= null
     then nb := Number_of_Factors(deco.all);
     else nb := 0;
    end if;
    if nb > 1 then
      loop
        New_Loop(fp,d,bl,bl_last,vl,vl_last,tl,tl_last,b,v,
                 gp,deco,nb,cnt,fail);
        nb := Number_of_Factors(deco.all);
        exit when (nb = 1) or (cnt >= threshold);
      end loop;
    end if;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( p : in DoblDobl_Complex_Polynomials.Poly;
                  ep : in DoblDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in DoblDobl_Complex_Vectors.Vector;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

    use Standard_Natural_VecVecs;
    use DoblDobl_Complex_VecLists;
    use DoblDobl_Complex_Poly_SysFun;

    cnt,nb : natural32 := 0;
    d : constant natural32 := natural32(gp'length);
    fp : constant Eval_Poly_Sys(0..b'last) := Eval_Sys(b'last,p,ep);
    bl,bl_last,vl,vl_last,tl,tl_last : List;
    fail : boolean := true;

  begin
    while fail loop
      First_Loop(fp,bl,bl_last,vl,vl_last,tl,tl_last,b,v,gp,deco,nb,fail);
    end loop;
    if deco /= null
     then nb := Number_of_Factors(deco.all);
     else nb := 0;
    end if;
    if nb > 1 then
      loop
        New_Loop(fp,d,bl,bl_last,vl,vl_last,tl,tl_last,b,v,
                 gp,deco,nb,cnt,fail);
        nb := Number_of_Factors(deco.all);
        exit when (nb = 1) or (cnt >= threshold);
      end loop;
    end if;
  end Monodromy_Breakup;

  procedure Monodromy_Breakup
                ( p : in QuadDobl_Complex_Polynomials.Poly;
                  ep : in QuadDobl_Complex_Poly_Functions.Eval_Poly;
                  b,v,gp : in QuadDobl_Complex_Vectors.Vector;
                  deco : in out Standard_Natural_VecVecs.Link_to_VecVec ) is

    use Standard_Natural_VecVecs;
    use QuadDobl_Complex_VecLists;
    use QuadDobl_Complex_Poly_SysFun;

    cnt,nb : natural32 := 0;
    d : constant natural32 := natural32(gp'length);
    fp : constant Eval_Poly_Sys(0..b'last) := Eval_Sys(b'last,p,ep);
    bl,bl_last,vl,vl_last,tl,tl_last : List;
    fail : boolean := true;

  begin
    while fail loop
      First_Loop(fp,bl,bl_last,vl,vl_last,tl,tl_last,b,v,gp,deco,nb,fail);
    end loop;
    if deco /= null
     then nb := Number_of_Factors(deco.all);
     else nb := 0;
    end if;
    if nb > 1 then
      loop
        New_Loop(fp,d,bl,bl_last,vl,vl_last,tl,tl_last,b,v,
                 gp,deco,nb,cnt,fail);
        nb := Number_of_Factors(deco.all);
        exit when (nb = 1) or (cnt >= threshold);
      end loop;
    end if;
  end Monodromy_Breakup;

end Monodromy_Polynomial_Breakup;
