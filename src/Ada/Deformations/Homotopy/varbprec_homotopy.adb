with String_Splitters;                   use String_Splitters;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Multprec_Floating_Numbers;
with Symbol_Table;
with Standard_Complex_Poly_Strings;
with DoblDobl_Complex_Poly_Strings;
with QuadDobl_Complex_Poly_Strings;
with Multprec_Complex_Poly_Strings;
with Standard_Homotopy;
with DoblDobl_Homotopy;
with QuadDobl_Homotopy;
with Multprec_Homotopy;

package body Varbprec_Homotopy is

-- INTERNAL DATA :

  start,target : Link_to_Array_of_Strings;
  st_start,st_target : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  dd_start,dd_target : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
  qd_start,qd_target : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
  mp_start,mp_target : Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
  gamma : Standard_Complex_Numbers.Complex_Number;
  exp4t : natural32;
  standard_homotopy_initialized : boolean; -- is standard homotopy defined ?
  dobldobl_homotopy_initialized : boolean; -- is dobldobl homotopy defined ?
  quaddobl_homotopy_initialized : boolean; -- is quaddobl homotopy defined ?
  multprec_homotopy_numbsize : natural32;
  -- size of the numbers in the multiprecision homotopy

-- AUXILIARY CONSTRUCTOR :

  procedure Initialize_Standard_Homotopy is

  -- DESCRIPTION :
  --   Initializes the standard homotopy, using the string representations
  --   for the start and target polynomial systems.

    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_Strings;

    nvr : natural32;

  begin
    Standard_Homotopy.Clear;
    if start /= null then
      nvr := natural32(start'last);
      if Symbol_Table.Number < nvr + 1
       then Symbol_Table.Init(nvr+1);
      end if;
      st_start := new Poly_Sys'(Parse(nvr,start.all));
      if target /= null
       then st_target := new Poly_Sys'(Parse(nvr,target.all));
      end if;
      Standard_Homotopy.Create(st_target.all,st_start.all,exp4t,gamma);
      standard_homotopy_initialized := true;
    end if;
  end Initialize_Standard_Homotopy;

  procedure Initialize_DoblDobl_Homotopy is

  -- DESCRIPTION :
  --   Initializes the homotopy for double double precision,
  --   using the string representations
  --   for the start and target polynomial systems.

    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_Strings;

    nvr : natural32;
    regamma : constant double_float
            := Standard_Complex_Numbers.REAL_PART(gamma);
    imgamma : constant double_float
            := Standard_Complex_Numbers.IMAG_PART(gamma);
    ddregamma : constant double_double
              := Double_Double_Numbers.create(regamma);
    ddimgamma : constant double_double
              := Double_Double_Numbers.create(imgamma);
    ddgamma : constant DoblDobl_Complex_Numbers.Complex_Number
            := DoblDobl_Complex_Numbers.Create(ddregamma,ddimgamma);

  begin
    DoblDobl_Homotopy.Clear;
    if start /= null then
      nvr := natural32(start'last);
      if Symbol_Table.Number < nvr + 1
       then Symbol_Table.Init(nvr+1);
      end if;
      dd_start := new Poly_Sys'(Parse(nvr,start.all));
      if target /= null
       then dd_target := new Poly_Sys'(Parse(nvr,target.all));
      end if;
      DoblDobl_Homotopy.Create(dd_target.all,dd_start.all,exp4t,ddgamma);
      dobldobl_homotopy_initialized := true;
    end if;
  end Initialize_DoblDobl_Homotopy;

  procedure Initialize_QuadDobl_Homotopy is

  -- DESCRIPTION :
  --   Initializes the homotopy for quad double precision,
  --   using the string representations
  --   for the start and target polynomial systems.

    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_Strings;

    nvr : natural32;
    regamma : constant double_float
            := Standard_Complex_Numbers.REAL_PART(gamma);
    imgamma : constant double_float
            := Standard_Complex_Numbers.IMAG_PART(gamma);
    qdregamma : constant quad_double
              := Quad_Double_Numbers.create(regamma);
    qdimgamma : constant quad_double
              := Quad_Double_Numbers.create(imgamma);
    qdgamma : constant QuadDobl_Complex_Numbers.Complex_Number
            := QuadDobl_Complex_Numbers.Create(qdregamma,qdimgamma);

  begin
    QuadDobl_Homotopy.Clear;
    if start /= null then
      nvr := natural32(start'last);
      if Symbol_Table.Number < nvr + 1
       then Symbol_Table.Init(nvr+1);
      end if;
      qd_start := new Poly_Sys'(Parse(nvr,start.all));
      if target /= null
       then qd_target := new Poly_Sys'(Parse(nvr,target.all));
      end if;
      QuadDobl_Homotopy.Create(qd_target.all,qd_start.all,exp4t,qdgamma);
      quaddobl_homotopy_initialized := true;
    end if;
  end Initialize_QuadDobl_Homotopy;

  procedure Initialize_Multprec_Homotopy ( size : in natural32 ) is

  -- DESCRIPTION :
  --   Parses start and target system with respect to the precision
  --   as defined by the size of the numbers given in size.

    use Multprec_Complex_Poly_Systems;
    use Multprec_Complex_Poly_Strings;

    nvr : natural32;
    regamma : constant double_float
            := Standard_Complex_Numbers.REAL_PART(gamma);
    imgamma : constant double_float
            := Standard_Complex_Numbers.IMAG_PART(gamma);
    mpregamma : Multprec_Floating_Numbers.Floating_Number
              := Multprec_Floating_Numbers.create(regamma);
    mpimgamma : Multprec_Floating_Numbers.Floating_Number
              := Multprec_Floating_Numbers.create(imgamma);
    mpgamma : Multprec_Complex_Numbers.Complex_Number
            := Multprec_Complex_Numbers.Create(mpregamma,mpimgamma);

  begin
    Multprec_Homotopy.Clear;
    Multprec_Complex_Poly_Systems.Clear(mp_start);
    Multprec_Complex_Poly_Systems.Clear(mp_target);
    if start /= null then
      nvr := natural32(start'last);
      if Symbol_Table.Number < nvr + 1
       then Symbol_Table.Init(nvr+1);
      end if;
      mp_start := new Poly_Sys'(Parse(nvr,size,start.all));
      if target /= null
       then mp_target := new Poly_Sys'(Parse(nvr,size,target.all));
      end if;
      Multprec_Homotopy.Create(mp_target.all,mp_start.all,exp4t,mpgamma);
      multprec_homotopy_numbsize := size;
    end if;
  end Initialize_Multprec_Homotopy;

-- CONSTRUCTORS :

  procedure Create ( p,q : in Link_to_Array_of_Strings; k : in natural32;
                     a : in Standard_Complex_Numbers.Complex_Number ) is
  begin
    gamma := a;
    exp4t := k;
    if start /= null
     then Clear(start);
    end if;
    start := new Array_of_Strings(q'range);
    for i in p'range loop
      start(i) := new string'(q(i).all);
    end loop;
    if target /= null
     then Clear(target);
    end if;
    target := new Array_of_Strings(p'range);
    for i in p'range loop
      target(i) := new string'(p(i).all);
    end loop;
  end Create;

-- SELECTORS :

  function Standard_Homotopy_System
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys is

    use Standard_Complex_Poly_Systems;
    res : Link_to_Poly_Sys;

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := new Poly_Sys'(Standard_Homotopy.Homotopy_System);
    end if;
    return res;
  end Standard_Homotopy_System;

  function DoblDobl_Homotopy_System
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys is

    use DoblDobl_Complex_Poly_Systems;
    res : Link_to_Poly_Sys;

  begin
    if not dobldobl_homotopy_initialized
     then Initialize_DoblDobl_Homotopy;
    end if;
    if dobldobl_homotopy_initialized
     then res := new Poly_Sys'(DoblDobl_Homotopy.Homotopy_System);
    end if;
    return res;
  end DoblDobl_Homotopy_System;

  function QuadDobl_Homotopy_System
             return QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys is

    use QuadDobl_Complex_Poly_Systems;
    res : Link_to_Poly_Sys;

  begin
    if not quaddobl_homotopy_initialized
     then Initialize_QuadDobl_Homotopy;
    end if;
    if quaddobl_homotopy_initialized
     then res := new Poly_Sys'(QuadDobl_Homotopy.Homotopy_System);
    end if;
    return res;
  end QuadDobl_Homotopy_System;

  function Multprec_Homotopy_System ( deci : natural32 )
             return Multprec_Complex_Poly_Systems.Link_to_Poly_Sys is

    use Multprec_Complex_Poly_Systems;
    res : Link_to_Poly_Sys;
    size : constant natural32
         := Multprec_Floating_Numbers.Decimal_to_Size(deci);

  begin
    if multprec_homotopy_numbsize /= size
     then Initialize_Multprec_Homotopy(size);
    end if;
    if multprec_homotopy_numbsize = size
     then res := new Poly_Sys'(Multprec_Homotopy.Homotopy_System);
    end if;
    return res;
  end Multprec_Homotopy_System;

  function Eval ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := Standard_Homotopy.Eval(x,t);
    end if;
    return res;
  end Eval;

  function Eval ( x : DoblDobl_Complex_Vectors.Vector;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(x'range);

  begin
    if not dobldobl_homotopy_initialized
     then Initialize_DoblDobl_Homotopy;
    end if;
    if dobldobl_homotopy_initialized
     then res := DoblDobl_Homotopy.Eval(x,t);
    end if;
    return res;
  end Eval;

  function Eval ( x : QuadDobl_Complex_Vectors.Vector;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(x'range);

  begin
    if not quaddobl_homotopy_initialized
     then Initialize_QuadDobl_Homotopy;
    end if;
    if quaddobl_homotopy_initialized
     then res := QuadDobl_Homotopy.Eval(x,t);
    end if;
    return res;
  end Eval;

  function Eval ( x : Multprec_Complex_Vectors.Vector;
                  t : Multprec_Complex_Numbers.Complex_Number;
                  d : natural32 )
                return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(x'range);
    size : constant natural32 := Multprec_Floating_Numbers.Decimal_to_Size(d);

  begin
    if multprec_homotopy_numbsize /= size
     then Initialize_Multprec_Homotopy(size);
    end if;
    if multprec_homotopy_numbsize = size
     then res := Multprec_Homotopy.Eval(x,t);
    end if;
    return res;
  end Eval;

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := Standard_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

  function Diff ( x : DoblDobl_Complex_Vectors.Vector;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(x'range);

  begin
    if not dobldobl_homotopy_initialized
     then Initialize_DoblDobl_Homotopy;
    end if;
    if dobldobl_homotopy_initialized
     then res := DoblDobl_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

  function Diff ( x : QuadDobl_Complex_Vectors.Vector;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(x'range);

  begin
    if not quaddobl_homotopy_initialized
     then Initialize_QuadDobl_Homotopy;
    end if;
    if quaddobl_homotopy_initialized
     then res := QuadDobl_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

  function Diff ( x : Multprec_Complex_Vectors.Vector;
                  t : Multprec_Complex_Numbers.Complex_Number;
                  d : natural32 )
                return Multprec_Complex_Vectors.Vector is

    res : Multprec_Complex_Vectors.Vector(x'range);
    size : constant natural32 := Multprec_Floating_Numbers.Decimal_to_Size(d);

  begin
    if multprec_homotopy_numbsize /= size
     then Initialize_Multprec_Homotopy(size);
    end if;
    if multprec_homotopy_numbsize = size
     then res := Multprec_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

  function Diff ( x : Standard_Complex_Vectors.Vector;
                  t : Standard_Complex_Numbers.Complex_Number )
                return Standard_Complex_Matrices.Matrix is

    res : Standard_Complex_Matrices.Matrix(x'range,x'range);

  begin
    if not standard_homotopy_initialized
     then Initialize_Standard_Homotopy;
    end if;
    if standard_homotopy_initialized
     then res := Standard_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

  function Diff ( x : DoblDobl_Complex_Vectors.Vector;
                  t : DoblDobl_Complex_Numbers.Complex_Number )
                return DoblDobl_Complex_Matrices.Matrix is

    res : DoblDobl_Complex_Matrices.Matrix(x'range,x'range);

  begin
    if not dobldobl_homotopy_initialized
     then Initialize_DoblDobl_Homotopy;
    end if;
    if dobldobl_homotopy_initialized
     then res := DoblDobl_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

  function Diff ( x : QuadDobl_Complex_Vectors.Vector;
                  t : QuadDobl_Complex_Numbers.Complex_Number )
                return QuadDobl_Complex_Matrices.Matrix is

    res : QuadDobl_Complex_Matrices.Matrix(x'range,x'range);

  begin
    if not quaddobl_homotopy_initialized
     then Initialize_QuadDobl_Homotopy;
    end if;
    if quaddobl_homotopy_initialized
     then res := QuadDobl_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;


  function Diff ( x : Multprec_Complex_Vectors.Vector;
                  t : Multprec_Complex_Numbers.Complex_Number;
                  d : natural32 )
                return Multprec_Complex_Matrices.Matrix is

    res : Multprec_Complex_Matrices.Matrix(x'range,x'range);
    size : constant natural32 := Multprec_Floating_Numbers.Decimal_to_Size(d);

  begin
    if multprec_homotopy_numbsize /= size
     then Initialize_Multprec_Homotopy(size);
    end if;
    if multprec_homotopy_numbsize = size
     then res := Multprec_Homotopy.Diff(x,t);
    end if;
    return res;
  end Diff;

-- DESTRUCTOR :

  procedure Clear is
  begin
    if start /= null
     then Clear(start);
    end if;
    if target /= null
     then Clear(target);
    end if;
    if standard_homotopy_initialized then
      Standard_Homotopy.Clear;
      Standard_Complex_Poly_Systems.Clear(st_start);
      Standard_Complex_Poly_Systems.Clear(st_target);
    end if;
    if dobldobl_homotopy_initialized then
      DoblDobl_Homotopy.Clear;
      DoblDobl_Complex_Poly_Systems.Clear(dd_start);
      DoblDobl_Complex_Poly_Systems.Clear(dd_target);
    end if;
    if quaddobl_homotopy_initialized then
      QuadDobl_Homotopy.Clear;
      QuadDobl_Complex_Poly_Systems.Clear(qd_start);
      QuadDobl_Complex_Poly_Systems.Clear(qd_target);
    end if;
    if multprec_homotopy_numbsize /= 0 then
      Multprec_Homotopy.Clear;
      Multprec_Complex_Poly_Systems.Clear(mp_start);
      Multprec_Complex_Poly_Systems.Clear(mp_target);
    end if;
  end Clear;

begin
  start := null;
  target := null;
  standard_homotopy_initialized := false;
  dobldobl_homotopy_initialized := false;
  quaddobl_homotopy_initialized := false;
  multprec_homotopy_numbsize := 0;
end Varbprec_Homotopy;
