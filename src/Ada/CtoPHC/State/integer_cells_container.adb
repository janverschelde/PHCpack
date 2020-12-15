-- for exception handlers only :
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;
-- normal dependencies :
with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Standard_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_Functions;
with Exponent_Vectors;                   use Exponent_Vectors;
with Random_Coefficient_Systems;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_JacoMats;
with Standard_Poly_Laur_Convertors;
with DoblDobl_Complex_Laur_Systems;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_JacoMats;
with DoblDobl_Poly_Laur_Convertors;
with QuadDobl_Complex_Laur_Systems;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_JacoMats;
with QuadDobl_Poly_Laur_Convertors;
with Transforming_Laurent_Systems;
with Standard_Simpomial_Solvers;
with DoblDobl_Simpomial_Solvers;
with QuadDobl_Simpomial_Solvers;
with Floating_Polyhedral_Continuation;   use Floating_Polyhedral_Continuation;
with DoblDobl_Polyhedral_Continuation;   use DoblDobl_Polyhedral_Continuation;
with QuadDobl_Polyhedral_Continuation;   use QuadDobl_Polyhedral_Continuation;
with Drivers_for_Static_Lifting;
with Black_Mixed_Volume_Computations;
with PHCpack_Operations;

package body Integer_Cells_Container is

-- DATA STRUCTURES :

  mix : Standard_Integer_Vectors.Link_to_Vector;
  lifsup,lifsup_last : Link_to_Array_of_Lists;
  cells,cells_last : Mixed_Subdivision;
  st_rndcffsys : Standard_Complex_Poly_Systems.Link_to_Poly_Sys;
  dd_rndcffsys : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
  qd_rndcffsys : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys;
  st_lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
  dd_lq : DoblDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
  qd_lq : QuadDobl_Complex_Laur_Systems.Link_to_Laur_Sys;
  st_hq : Standard_Complex_Laur_SysFun.Link_to_Eval_Coeff_Laur_Sys;
  dd_hq : DoblDobl_Complex_Laur_SysFun.Link_to_Eval_Coeff_Laur_Sys;
  qd_hq : QuadDobl_Complex_Laur_SysFun.Link_to_Eval_Coeff_Laur_Sys;
  st_homty_expvec : Link_to_Exponent_Vectors_Array;
  dd_homty_expvec : Link_to_Exponent_Vectors_Array;
  qd_homty_expvec : Link_to_Exponent_Vectors_Array;
  st_homty_coeffv : Standard_Complex_VecVecs.Link_to_VecVec;
  dd_homty_coeffv : DoblDobl_Complex_VecVecs.Link_to_VecVec;
  qd_homty_coeffv : QuadDobl_Complex_VecVecs.Link_to_VecVec;
  st_homty_jacmat : Standard_Complex_Laur_JacoMats.Link_to_Eval_Coeff_Jaco_Mat;
  dd_homty_jacmat : DoblDobl_Complex_Laur_JacoMats.Link_to_Eval_Coeff_Jaco_Mat;
  qd_homty_jacmat : QuadDobl_Complex_Laur_JacoMats.Link_to_Eval_Coeff_Jaco_Mat;
  st_homty_mulfac : Standard_Complex_Laur_JacoMats.Link_to_Mult_Factors;
  dd_homty_mulfac : DoblDobl_Complex_Laur_JacoMats.Link_to_Mult_Factors;
  qd_homty_mulfac : QuadDobl_Complex_Laur_JacoMats.Link_to_Mult_Factors;
  st_start_sols,st_target_sols,st_target_last :
     Standard_Complex_Solutions.Link_to_Array_of_Solution_Lists;
  dd_start_sols,dd_target_sols,dd_target_last :
     DoblDobl_Complex_Solutions.Link_to_Array_of_Solution_Lists;
  qd_start_sols,qd_target_sols,qd_target_last :
     QuadDobl_Complex_Solutions.Link_to_Array_of_Solution_Lists;

-- AUXILIARIES :

  function Retrieve ( L : List; k : natural32 )
                    return Standard_Integer_Vectors.Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the k-th point of the list L.
  --   If k > Length_Of(L), then null is returned.

    tmp : List := L;

  begin
    for i in 1..(k-1) loop
      exit when Is_Null(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    if not Is_Null(tmp)
     then return Head_Of(tmp);
     else return null;
    end if;
  end Retrieve;

  function Position ( L : List; x : Standard_Integer_Vectors.Vector )
                    return natural32 is

  -- DESCRIPTION :
  --   Returns 0 if x does not occur in the list L,
  --   otherwise returns the position of x in the list L.

    res : natural32 := 0;
    tmp : List := L;
    p : Standard_Integer_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      res := res + 1;
      p := Head_Of(tmp);
      if Standard_Integer_Vectors.Equal(x,p.all)
       then return res;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end Position;

-- CREATORS :

  procedure Initialize_Supports ( nbr : in natural32 ) is
  begin
    if lifsup /= null then
      Deep_Clear(lifsup);
      if lifsup = null then
        lifsup := new Array_of_Lists(1..integer32(nbr));
        lifsup_last := new Array_of_Lists(1..integer32(nbr));
      end if;
      lifsup_last := null;
    end if;
  end Initialize_Supports;

  procedure Initialize 
              ( mixture : in Standard_Integer_Vectors.Link_to_Vector;
                lifting : in Link_to_Array_of_Lists;
                mcc : in Mixed_Subdivision ) is
  begin
    mix := mixture;
    lifsup := lifting;
    cells := mcc;
  end Initialize;

  procedure Initialize 
              ( mixture : in Standard_Integer_Vectors.Link_to_Vector ) is
  begin
    mix := mixture;
  end Initialize;

  procedure Initialize ( lifting : in Link_to_Array_of_Lists ) is
  begin
    lifsup := lifting;
  end Initialize;

  procedure Initialize ( mcc : in Mixed_Subdivision ) is
  begin
    cells := mcc;
  end Initialize;

  procedure Make_Subdivision is

    use Standard_Integer_Vectors;
    use Drivers_for_Static_Lifting;

  begin
    if lifsup /= null then
      if mix = null then
        declare
          dim : constant integer32 := integer32(Dimension_of_Supports);
          mixtype : constant Standard_Integer_Vectors.Vector(1..dim)
                  := (1..dim => 1);
        begin
          Integer_Create_Mixed_Cells(dim,mixtype,lifsup.all,cells);
        end;
      else
        declare
          dim : constant integer32 := integer32(Dimension_of_Supports);
        begin
          Integer_Create_Mixed_Cells(dim,mix.all,lifsup.all,cells);
        end;
      end if;
    end if;
  end Make_Subdivision;

  procedure Generate_Random_Standard_Coefficient_System is

    n : constant natural32 := Integer_Cells_Container.Dimension-1;
    q : constant Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n))
      := Random_Coefficient_Systems.Create(n,mix.all,lifsup.all);

  begin
    st_rndcffsys := new Standard_Complex_Poly_Systems.Poly_Sys'(q);
  end Generate_Random_Standard_Coefficient_System;

  procedure Generate_Random_DoblDobl_Coefficient_System is

    n : constant natural32 := Integer_Cells_Container.Dimension-1;
    q : constant DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n))
      := Random_Coefficient_Systems.Create(n,mix.all,lifsup.all);

  begin
    dd_rndcffsys := new DoblDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end Generate_Random_DoblDobl_Coefficient_System;

  procedure Generate_Random_QuadDobl_Coefficient_System is

    n : constant natural32 := Integer_Cells_Container.Dimension-1;
    q : constant QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n))
      := Random_Coefficient_Systems.Create(n,mix.all,lifsup.all);

  begin
    qd_rndcffsys := new QuadDobl_Complex_Poly_Systems.Poly_Sys'(q);
  end Generate_Random_QuadDobl_Coefficient_System;

  procedure Initialize_Random_Standard_Coefficient_System
              ( q : in Standard_Complex_Poly_Systems.Poly_Sys ) is
  begin
    st_rndcffsys := new Standard_Complex_Poly_Systems.Poly_Sys(q'range);
    for i in q'range loop
      Standard_Complex_Polynomials.Copy(q(i),st_rndcffsys(i));
    end loop;
  end Initialize_Random_Standard_Coefficient_System;

  procedure Initialize_Random_DoblDobl_Coefficient_System
              ( q : in DoblDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    dd_rndcffsys := new DoblDobl_Complex_Poly_Systems.Poly_Sys(q'range);
    for i in q'range loop
      DoblDobl_Complex_Polynomials.Copy(q(i),dd_rndcffsys(i));
    end loop;
  end Initialize_Random_DoblDobl_Coefficient_System;

  procedure Initialize_Random_QuadDobl_Coefficient_System
              ( q : in QuadDobl_Complex_Poly_Systems.Poly_Sys ) is
  begin
    qd_rndcffsys := new QuadDobl_Complex_Poly_Systems.Poly_Sys(q'range);
    for i in q'range loop
      QuadDobl_Complex_Polynomials.Copy(q(i),qd_rndcffsys(i));
    end loop;
  end Initialize_Random_QuadDobl_Coefficient_System;

  procedure Standard_Polyhedral_Homotopy is

    use Standard_Complex_Solutions;
    use Standard_Complex_Laur_Functions;
    use Standard_Complex_Laur_Systems;
    use Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_JacoMats;
    use Standard_Poly_Laur_Convertors;

    q : Standard_Complex_Poly_Systems.Link_to_Poly_Sys renames st_rndcffsys;
    use Standard_Complex_Vectors,Standard_Complex_VecVecs;
    len : constant integer32
        := integer32(Integer_Cells_Container.Length);

  begin
    if st_lq /= null
     then Clear(st_lq);
    end if;
    st_lq := new Laur_Sys'(Polynomial_to_Laurent_System(q.all));
    if st_hq /= null
     then Clear(st_hq);
    end if;
    st_hq := new Eval_Coeff_Laur_Sys'(Create(st_lq.all));
    if st_homty_expvec /= null
     then Clear(st_homty_expvec);
    end if;
    st_homty_expvec := new Exponent_Vectors_Array'(Create(q.all));
    if st_homty_coeffv /= null then
      for i in st_homty_coeffv'range loop
        if st_homty_coeffv(i) /= null
         then Clear(st_homty_coeffv(i));
        end if;
      end loop;
    end if;
    st_homty_coeffv := new Standard_Complex_VecVecs.VecVec(q'range);
    for i in q'range loop
      declare
        c : constant Standard_Complex_Vectors.Vector := Coeff(st_lq(i));
      begin
        st_homty_coeffv(i) := new Standard_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          st_homty_coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    if st_homty_jacmat /= null
     then Clear(st_homty_jacmat);
    end if;
    st_homty_jacmat := new Eval_Coeff_Jaco_Mat(q'range,q'first..q'last+1);
    if st_homty_mulfac /= null
     then Clear(st_homty_mulfac);
    end if;
    st_homty_mulfac
      := new Mult_Factors(st_homty_jacmat'range(1),st_homty_jacmat'range(2));
    Create(st_lq.all,st_homty_jacmat.all,st_homty_mulfac.all);
    if st_start_sols /= null
     then Clear(st_start_sols);
    end if;
    st_start_sols := new Array_of_Solution_Lists(1..len);
    if st_target_sols /= null
     then Clear(st_target_sols);
    end if;
    st_target_sols := new Array_of_Solution_Lists(1..len);
    st_target_last := new Array_of_Solution_Lists(1..len);
  end Standard_Polyhedral_Homotopy;

  procedure DoblDobl_Polyhedral_Homotopy is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Laur_Functions;
    use DoblDobl_Complex_Laur_Systems;
    use DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_JacoMats;
    use DoblDobl_Poly_Laur_Convertors;

    q : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys renames dd_rndcffsys;
    use DoblDobl_Complex_Vectors,DoblDobl_Complex_VecVecs;
    len : constant integer32
        := integer32(Integer_Cells_Container.Length);

  begin
    if dd_lq /= null
     then Clear(dd_lq);
    end if;
    dd_lq := new Laur_Sys'(Polynomial_to_Laurent_System(q.all));
    if dd_hq /= null
     then Clear(dd_hq);
    end if;
    dd_hq := new Eval_Coeff_Laur_Sys'(Create(dd_lq.all));
    if dd_homty_expvec /= null
     then Clear(dd_homty_expvec);
    end if;
    dd_homty_expvec := new Exponent_Vectors_Array'(Create(q.all));
    if dd_homty_coeffv /= null then
      for i in dd_homty_coeffv'range loop
        if dd_homty_coeffv(i) /= null
         then Clear(dd_homty_coeffv(i));
        end if;
      end loop;
    end if;
    dd_homty_coeffv := new DoblDobl_Complex_VecVecs.VecVec(q'range);
    for i in q'range loop
      declare
        c : constant DoblDobl_Complex_Vectors.Vector := Coeff(dd_lq(i));
      begin
        dd_homty_coeffv(i) := new DoblDobl_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          dd_homty_coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    if dd_homty_jacmat /= null
     then Clear(dd_homty_jacmat);
    end if;
    dd_homty_jacmat := new Eval_Coeff_Jaco_Mat(q'range,q'first..q'last+1);
    if dd_homty_mulfac /= null
     then Clear(dd_homty_mulfac);
    end if;
    dd_homty_mulfac
      := new Mult_Factors(dd_homty_jacmat'range(1),dd_homty_jacmat'range(2));
    Create(dd_lq.all,dd_homty_jacmat.all,dd_homty_mulfac.all);
    if dd_start_sols /= null
     then Clear(dd_start_sols);
    end if;
    dd_start_sols := new Array_of_Solution_Lists(1..len);
    if dd_target_sols /= null
     then Clear(dd_target_sols);
    end if;
    dd_target_sols := new Array_of_Solution_Lists(1..len);
    dd_target_last := new Array_of_Solution_Lists(1..len);
  end DoblDobl_Polyhedral_Homotopy;

  procedure QuadDobl_Polyhedral_Homotopy is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Laur_Functions;
    use QuadDobl_Complex_Laur_Systems;
    use QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_JacoMats;
    use QuadDobl_Poly_Laur_Convertors;

    q : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys renames qd_rndcffsys;
    use QuadDobl_Complex_Vectors,QuadDobl_Complex_VecVecs;
    len : constant integer32
        := integer32(Integer_Cells_Container.Length);

  begin
    if qd_lq /= null
     then Clear(qd_lq);
    end if;
    qd_lq := new Laur_Sys'(Polynomial_to_Laurent_System(q.all));
    if qd_hq /= null
     then Clear(qd_hq);
    end if;
    qd_hq := new Eval_Coeff_Laur_Sys'(Create(qd_lq.all));
    if qd_homty_expvec /= null
     then Clear(qd_homty_expvec);
    end if;
    qd_homty_expvec := new Exponent_Vectors_Array'(Create(q.all));
    if qd_homty_coeffv /= null then
      for i in qd_homty_coeffv'range loop
        if qd_homty_coeffv(i) /= null
         then Clear(qd_homty_coeffv(i));
        end if;
      end loop;
    end if;
    qd_homty_coeffv := new QuadDobl_Complex_VecVecs.VecVec(q'range);
    for i in q'range loop
      declare
        c : constant QuadDobl_Complex_Vectors.Vector := Coeff(qd_lq(i));
      begin
        qd_homty_coeffv(i) := new QuadDobl_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          qd_homty_coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    if qd_homty_jacmat /= null
     then Clear(qd_homty_jacmat);
    end if;
    qd_homty_jacmat := new Eval_Coeff_Jaco_Mat(q'range,q'first..q'last+1);
    if qd_homty_mulfac /= null
     then Clear(qd_homty_mulfac);
    end if;
    qd_homty_mulfac
      := new Mult_Factors(qd_homty_jacmat'range(1),qd_homty_jacmat'range(2));
    Create(qd_lq.all,qd_homty_jacmat.all,qd_homty_mulfac.all);
    if qd_start_sols /= null
     then Clear(qd_start_sols);
    end if;
    qd_start_sols := new Array_of_Solution_Lists(1..len);
    if qd_target_sols /= null
     then Clear(qd_target_sols);
    end if;
    qd_target_sols := new Array_of_Solution_Lists(1..len);
    qd_target_last := new Array_of_Solution_Lists(1..len);
  end QuadDobl_Polyhedral_Homotopy;

-- SELECTORS :

  function Length return natural32 is
  begin
    return Length_Of(cells);
  end Length;

  function Dimension return natural32 is

    use Standard_Integer_Vectors;

    mic : Mixed_Cell;

  begin
    if Is_Null(cells) then
      return 0;
    else
      mic := Head_Of(cells);
      if mic.nor = null
       then return 0;
       else return natural32(mic.nor'last);
      end if;
    end if;
  end Dimension;

  function Dimension_of_Supports return natural32 is

    res : natural32 := 0;
    pts : List;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    if lifsup /= null then
      pts := lifsup(lifsup'first);
      if not Is_Null(pts) then
        lpt := Head_Of(pts);
        res := natural32(lpt'last) - 1;
      end if;
    end if;
    return res;
  end Dimension_of_Supports;

  function Type_of_Mixture return Standard_Integer_Vectors.Link_to_Vector is
  begin
    return mix;
  end Type_of_Mixture;

  function Lifted_Supports return Link_to_Array_of_Lists is
  begin
    return lifsup;
  end Lifted_Supports;

  procedure Retrieve ( k : in natural32; mic : out Mixed_Cell;
                       fail : out boolean ) is

    tmp : Mixed_Subdivision := cells;
    cnt : natural32 := 0;

  begin
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      if cnt = k then
        mic := Head_Of(tmp);
        fail := false;
        return;
      else
        tmp := Tail_Of(tmp);
      end if;
    end loop;
    fail := true;
  end Retrieve;

  procedure Retrieve_Mixed_Cell
             ( k : in natural32; fail : out boolean;
               cnt,lab : out Standard_Integer_Vectors.Link_to_Vector;
               normal : out Standard_Integer_Vectors.Link_to_Vector ) is

    mic : Mixed_Cell;
    sum,ind : integer32 := 0;
    tmp : List;

  begin
    Retrieve(k,mic,fail);
    if not fail then
      normal := mic.nor;
      cnt := new Standard_Integer_Vectors.Vector(mic.pts'range);
      for i in mic.pts'range loop
        cnt(i) := integer32(Length_Of(mic.pts(i)));
        sum := sum + cnt(i);
      end loop;
      lab := new Standard_Integer_Vectors.Vector(1..sum);
      for i in mic.pts'range loop
        tmp := mic.pts(i);
        while not Is_Null(tmp) loop
          ind := ind + 1;
          lab(ind) := integer32(Position(lifsup(i),Head_Of(tmp).all));
          tmp := Tail_Of(tmp);
        end loop;
      end loop;
    end if;
  end Retrieve_Mixed_Cell;

  function Retrieve return Mixed_Subdivision is
  begin
    return cells;
  end Retrieve;

  function Retrieve_Random_Standard_Coefficient_System
             return Standard_Complex_Poly_Systems.Poly_Sys is
  begin
    return st_rndcffsys.all;
  end Retrieve_Random_Standard_Coefficient_System;

  function Retrieve_Random_Standard_Coefficient_System
             return Standard_Complex_Poly_Systems.Link_to_Poly_Sys is
  begin
    return st_rndcffsys;
  end Retrieve_Random_Standard_Coefficient_System;

  function Retrieve_Random_DoblDobl_Coefficient_System
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is
  begin
    return dd_rndcffsys.all;
  end Retrieve_Random_DoblDobl_Coefficient_System;

  function Retrieve_Random_DoblDobl_Coefficient_System
             return DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys is
  begin
    return dd_rndcffsys;
  end Retrieve_Random_DoblDobl_Coefficient_System;

  function Retrieve_Random_QuadDobl_Coefficient_System
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is
  begin
    return qd_rndcffsys.all;
  end Retrieve_Random_QuadDobl_Coefficient_System;

  function Retrieve_Random_QuadDobl_Coefficient_System
             return QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys is
  begin
    return qd_rndcffsys;
  end Retrieve_Random_QuadDobl_Coefficient_System;

  function Retrieve_Standard_Start_Solution
             ( k,i : natural32 )
             return Standard_Complex_Solutions.Link_to_Solution is
  begin
    if Standard_Complex_Solutions.Is_Null(st_start_sols(integer32(k)))
     then return null;
     else return Standard_Complex_Solutions.Retrieve
                   (st_start_sols(integer32(k)),i);
    end if;
  end Retrieve_Standard_Start_Solution;

  function Retrieve_Standard_Target_Solution
             ( k,i : natural32 )
             return Standard_Complex_Solutions.Link_to_Solution is
  begin
    if Standard_Complex_Solutions.Is_Null(st_target_sols(integer32(k)))
     then return null;
     else return Standard_Complex_Solutions.Retrieve
                   (st_target_sols(integer32(k)),i);
    end if;
  end Retrieve_Standard_Target_Solution;

  function Retrieve_DoblDobl_Start_Solution
             ( k,i : natural32 )
             return DoblDobl_Complex_Solutions.Link_to_Solution is
  begin
    if DoblDobl_Complex_Solutions.Is_Null(dd_start_sols(integer32(k)))
     then return null;
     else return DoblDobl_Complex_Solutions.Retrieve
                   (dd_start_sols(integer32(k)),i);
    end if;
  end Retrieve_DoblDobl_Start_Solution;

  function Retrieve_DoblDobl_Target_Solution
             ( k,i : natural32 )
             return DoblDobl_Complex_Solutions.Link_to_Solution is
  begin
    if DoblDobl_Complex_Solutions.Is_Null(dd_target_sols(integer32(k)))
     then return null;
     else return DoblDobl_Complex_Solutions.Retrieve
                   (dd_target_sols(integer32(k)),i);
    end if;
  end Retrieve_DoblDobl_Target_Solution;

  function Retrieve_QuadDobl_Start_Solution
             ( k,i : natural32 )
             return QuadDobl_Complex_Solutions.Link_to_Solution is
  begin
    if QuadDobl_Complex_Solutions.Is_Null(qd_start_sols(integer32(k)))
     then return null;
     else return QuadDobl_Complex_Solutions.Retrieve
                   (qd_start_sols(integer32(k)),i);
    end if;
  end Retrieve_QuadDobl_Start_Solution;

  function Retrieve_QuadDobl_Target_Solution
             ( k,i : natural32 )
             return QuadDobl_Complex_Solutions.Link_to_Solution is
  begin
    if QuadDobl_Complex_Solutions.Is_Null(qd_target_sols(integer32(k)))
     then return null;
     else return QuadDobl_Complex_Solutions.Retrieve
                   (qd_target_sols(integer32(k)),i);
    end if;
  end Retrieve_QuadDobl_Target_Solution;

-- CONSTRUCTOR :

  function Append_to_Support
             ( k : in natural32;
               x : in Standard_Integer_Vectors.Vector ) 
             return boolean is

    r : integer32;
    use Standard_Integer_Vectors;

  begin
    if mix = null then
      return true;
    else 
      r := mix'last;
      if integer32(k) > r then
        return true;
      else
        if lifsup = null then
          lifsup := new Array_of_Lists(1..r);
          lifsup_last := new Array_of_Lists(1..r);
        end if;
        Append(lifsup(integer32(k)),lifsup_last(integer32(k)),x);
        return false;
      end if;
    end if;
  end Append_to_Support;

  procedure Append ( mic : in Mixed_Cell ) is
  begin
    Append(cells,cells_last,mic);
  end Append;

  procedure Append_Mixed_Cell
             ( cnt,lab : in Standard_Integer_Vectors.Vector;
               normal : in Standard_Integer_Vectors.Vector ) is

    mic : Mixed_Cell;
    last : Array_of_Lists(cnt'range);
    ind : integer32;
    lpt : Standard_Integer_Vectors.Link_to_Vector;

  begin
    mic.nor := new Standard_Integer_Vectors.Vector'(normal);
    mic.pts := new Array_of_Lists(cnt'range);
    ind := lab'first;
    for i in cnt'range loop
      for j in 1..cnt(i) loop
        lpt := Retrieve(lifsup(i),natural32(lab(ind)));
        ind := ind + 1;
        Append(mic.pts(i),last(i),lpt);
      end loop;
    end loop;
    Integer_Cells_Container.Append(mic);
  end Append_Mixed_Cell;

  procedure Solve_Standard_Start_System
              ( k : in natural32; mv : out natural32 ) is

    use Standard_Complex_Solutions;
    use Standard_Poly_Laur_Convertors;

    mic : Mixed_Cell;
    tol_zero : constant double_float := 1.0E-12;
    fail,zero_y : boolean;

  begin
   -- put("Before retrieval of mixed cell k = "); put(k,1); put_line(" ...");
    Retrieve(k,mic,fail);
    if fail then
      mv := 0;
    else
     -- put("... retrieval of cell "); put(k,1); put_line(" succeeded.");
      declare
        q : Standard_Complex_Poly_Systems.Link_to_Poly_Sys
            renames st_rndcffsys;
        sub_q : Standard_Complex_Poly_Systems.Poly_Sys(q'range);
        lau_q : Standard_Complex_Laur_Systems.Laur_Sys(q'range);
        sols : Solution_List;
      begin
        sub_q := Select_Terms(q.all,mix.all,mic.pts.all);
        lau_q := Polynomial_to_Laurent_System(sub_q);
        Transforming_Laurent_Systems.Shift(lau_q);
       -- put_line("The supported start system : "); put_line(sub_q);
        Standard_Simpomial_Solvers.Solve(lau_q,tol_zero,sols,fail,zero_y);
        if fail then
          mv := 0;
         -- put_line("The system is not a fewnomial system!");
        else
          mv := Length_Of(sols);
          st_start_sols(integer32(k)) := sols;
        end if;
        Standard_Complex_Poly_Systems.Clear(sub_q);
        Standard_Complex_Laur_Systems.Clear(lau_q);
      exception
        when others =>
          put_line("Exception raised for supported subsystem :");
          put_line(sub_q);
          raise;
      end;
    end if;
  exception
    when others =>
      put_line("Exception raised for mixed cell with supports :");
      put(mic.pts.all);
      raise;
  end Solve_Standard_Start_System;

  procedure Solve_DoblDobl_Start_System
              ( k : in natural32; mv : out natural32 ) is

    use DoblDobl_Complex_Solutions;
    use DoblDobl_Poly_Laur_Convertors;

    mic : Mixed_Cell;
    tol_zero : constant double_double := create(1.0E-12);
    fail,zero_y : boolean;

  begin
   -- put("Before retrieval of mixed cell k = "); put(k,1); put_line(" ...");
    Retrieve(k,mic,fail);
    if fail then
      mv := 0;
    else
     -- put("... retrieval of cell "); put(k,1); put_line(" succeeded.");
      declare
        q : DoblDobl_Complex_Poly_Systems.Link_to_Poly_Sys
            renames dd_rndcffsys;
        sub_q : DoblDobl_Complex_Poly_Systems.Poly_Sys(q'range);
        lau_q : DoblDobl_Complex_Laur_Systems.Laur_Sys(q'range);
        sols : Solution_List;
      begin
        sub_q := Select_Terms(q.all,mix.all,mic.pts.all);
        lau_q := Polynomial_to_Laurent_System(sub_q);
        Transforming_Laurent_Systems.Shift(lau_q);
       -- put_line("The supported start system : "); put_line(sub_q);
        DoblDobl_Simpomial_Solvers.Solve(lau_q,tol_zero,sols,fail,zero_y);
        if fail then
          mv := 0;
         -- put_line("The system is not a fewnomial system!");
        else
          mv := Length_Of(sols);
          dd_start_sols(integer32(k)) := sols;
        end if;
        DoblDobl_Complex_Poly_Systems.Clear(sub_q);
        DoblDobl_Complex_Laur_Systems.Clear(lau_q);
      exception
        when others =>
          put_line("Exception raised for supported subsystem :");
          put_line(sub_q);
          raise;
      end;
    end if;
  exception
    when others =>
      put_line("Exception raised for mixed cell with supports :");
      put(mic.pts.all);
      raise;
  end Solve_DoblDobl_Start_System;

  procedure Solve_QuadDobl_Start_System
              ( k : in natural32; mv : out natural32 ) is

    use QuadDobl_Complex_Solutions;
    use QuadDobl_Poly_Laur_Convertors;

    mic : Mixed_Cell;
    tol_zero : constant quad_double := create(1.0E-12);
    fail,zero_y : boolean;

  begin
   -- put("Before retrieval of mixed cell k = "); put(k,1); put_line(" ...");
    Retrieve(k,mic,fail);
    if fail then
      mv := 0;
    else
     -- put("... retrieval of cell "); put(k,1); put_line(" succeeded.");
      declare
        q : QuadDobl_Complex_Poly_Systems.Link_to_Poly_Sys
            renames qd_rndcffsys;
        sub_q : QuadDobl_Complex_Poly_Systems.Poly_Sys(q'range);
        lau_q : QuadDobl_Complex_Laur_Systems.Laur_Sys(q'range);
        sols : Solution_List;
      begin
        sub_q := Select_Terms(q.all,mix.all,mic.pts.all);
        lau_q := Polynomial_to_Laurent_System(sub_q);
        Transforming_Laurent_Systems.Shift(lau_q);
       -- put_line("The supported start system : "); put_line(sub_q);
        QuadDobl_Simpomial_Solvers.Solve(lau_q,tol_zero,sols,fail,zero_y);
        if fail then
          mv := 0;
         -- put_line("The system is not a fewnomial system!");
        else
          mv := Length_Of(sols);
          qd_start_sols(integer32(k)) := sols;
        end if;
        QuadDobl_Complex_Poly_Systems.Clear(sub_q);
        QuadDobl_Complex_Laur_Systems.Clear(lau_q);
      exception
        when others =>
          put_line("Exception raised for supported subsystem :");
          put_line(sub_q);
          raise;
      end;
    end if;
  exception
    when others =>
      put_line("Exception raised for mixed cell with supports :");
      put(mic.pts.all);
      raise;
  end Solve_QuadDobl_Start_System;

  procedure Track_Standard_Solution_Path ( k,i,otp : in natural32 ) is

    use Standard_Complex_Solutions;

    fail : boolean;
    mic : Mixed_Cell;
    ls : constant Link_to_Solution := Retrieve_Standard_Start_Solution(k,i);
    ns : constant Link_to_Solution := new Solution'(ls.all);
    s : Solution_List;

  begin
    Retrieve(k,mic,fail);
    Construct(ns,s);
    if otp = 0 then
      null;
     -- Mixed_Continuation
     --   (mix.all,lifsup.all,st_hq.all,st_homty_coeffv.all,st_homty_expvec.all,
     --    st_homty_jacmat.all,st_homty_mulfac.all,mic.nor.all,s);
    else
      if PHCpack_Operations.Is_File_Defined then
        null;
       -- Mixed_Continuation(PHCpack_Operations.Output_File,mix.all,
       --    lifsup.all,st_hq.all,st_homty_coeffv.all,st_homty_expvec.all,
       --    st_homty_jacmat.all,st_homty_mulfac.all,mic.nor.all,s);
      else
        null;
       -- Mixed_Continuation(standard_output,mix.all,
       --    lifsup.all,st_hq.all,st_homty_coeffv.all,st_homty_expvec.all,
       --    st_homty_jacmat.all,st_homty_mulfac.all,mic.nor.all,s);
      end if;
    end if;
    Append(st_target_sols(integer32(k)),st_target_last(integer32(k)),
           Head_Of(s).all);
    Clear(s);
  end Track_Standard_Solution_Path;

  procedure Track_DoblDobl_Solution_Path ( k,i,otp : in natural32 ) is

    use DoblDobl_Complex_Solutions;

    fail : boolean;
    mic : Mixed_Cell;
    ls : constant Link_to_Solution := Retrieve_DoblDobl_Start_Solution(k,i);
    ns : constant Link_to_Solution := new Solution'(ls.all);
    s : Solution_List;

  begin
    Retrieve(k,mic,fail);
    Construct(ns,s);
    if otp = 0 then
      null;
     -- Mixed_Continuation
     --   (mix.all,lifsup.all,dd_hq.all,dd_homty_coeffv.all,dd_homty_expvec.all,
     --    dd_homty_jacmat.all,dd_homty_mulfac.all,mic.nor.all,s);
    else
      if PHCpack_Operations.Is_File_Defined then
        null;
       -- Mixed_Continuation(PHCpack_Operations.Output_File,mix.all,
       --    lifsup.all,dd_hq.all,dd_homty_coeffv.all,dd_homty_expvec.all,
       --    dd_homty_jacmat.all,dd_homty_mulfac.all,mic.nor.all,s);
      else
        null;
       -- Mixed_Continuation(standard_output,mix.all,
       --    lifsup.all,dd_hq.all,dd_homty_coeffv.all,dd_homty_expvec.all,
       --    dd_homty_jacmat.all,dd_homty_mulfac.all,mic.nor.all,s);
      end if;
    end if;
    Append(dd_target_sols(integer32(k)),dd_target_last(integer32(k)),
           Head_Of(s).all);
    Clear(s);
  end Track_DoblDobl_Solution_Path;

  procedure Track_QuadDobl_Solution_Path ( k,i,otp : in natural32 ) is

    use QuadDobl_Complex_Solutions;

    fail : boolean;
    mic : Mixed_Cell;
    ls : constant Link_to_Solution := Retrieve_QuadDobl_Start_Solution(k,i);
    ns : constant Link_to_Solution := new Solution'(ls.all);
    s : Solution_List;

  begin
    Retrieve(k,mic,fail);
    Construct(ns,s);
    if otp = 0 then
      null;
     -- Mixed_Continuation
     --   (mix.all,lifsup.all,qd_hq.all,qd_homty_coeffv.all,qd_homty_expvec.all,
     --    qd_homty_jacmat.all,qd_homty_mulfac.all,mic.nor.all,s);
    else
      if PHCpack_Operations.Is_File_Defined then
        null;
       -- Mixed_Continuation(PHCpack_Operations.Output_File,mix.all,
       --    lifsup.all,qd_hq.all,qd_homty_coeffv.all,qd_homty_expvec.all,
       --    qd_homty_jacmat.all,qd_homty_mulfac.all,mic.nor.all,s);
      else
        null;
       -- Mixed_Continuation(standard_output,mix.all,
       --    lifsup.all,qd_hq.all,qd_homty_coeffv.all,qd_homty_expvec.all,
       --    qd_homty_jacmat.all,qd_homty_mulfac.all,mic.nor.all,s);
      end if;
    end if;
    Append(qd_target_sols(integer32(k)),qd_target_last(integer32(k)),
           Head_Of(s).all);
    Clear(s);
  end Track_QuadDobl_Solution_Path;

  function Mixed_Volume return natural32 is

    use Black_Mixed_Volume_Computations;

    res : natural32 := 0;
    n : constant natural32 := Dimension_of_Supports;
    q : Standard_Complex_Laur_Systems.Laur_Sys(1..integer32(n))
      := Random_Coefficient_Systems.Create(n,mix.all,lifsup.all);
   -- perm,iprm : Standard_Integer_Vectors.Link_to_Vector;

  begin
    Deep_Clear(lifsup); Clear(cells);
   -- Black_Box_Mixed_Volume_Computation(q,mix,perm,iprm,lifsup,cells,res);
    Standard_Complex_Laur_Systems.Clear(q);
    return res;
  end Mixed_Volume;

-- DESTRUCTORS :

  procedure Clear_Cell_Data is
  begin
    Deep_Clear(lifsup);
    lifsup_last := null;
    Deep_Clear(cells);
    Standard_Integer_Vectors.Clear(mix);
  end Clear_Cell_Data;

  procedure Clear_Standard_Data is
  begin
    Clear(st_homty_expvec);
    Standard_Complex_Poly_Systems.Clear(st_rndcffsys);
    Standard_Complex_Laur_Systems.Clear(st_lq);
    Standard_Complex_Laur_SysFun.Clear(st_hq);
    Standard_Complex_VecVecs.Deep_Clear(st_homty_coeffv);
    Standard_Complex_Laur_JacoMats.Clear(st_homty_jacmat);
    Standard_Complex_Laur_JacoMats.Clear(st_homty_mulfac);
    Standard_Complex_Solutions.Clear(st_start_sols);
    Standard_Complex_Solutions.Clear(st_target_sols);
  end Clear_Standard_Data;

  procedure Clear_DoblDobl_Data is
  begin
    Clear(dd_homty_expvec);
    DoblDobl_Complex_Poly_Systems.Clear(dd_rndcffsys);
    DoblDobl_Complex_Laur_Systems.Clear(dd_lq);
    DoblDobl_Complex_Laur_SysFun.Clear(dd_hq);
    DoblDobl_Complex_VecVecs.Deep_Clear(dd_homty_coeffv);
    DoblDobl_Complex_Laur_JacoMats.Clear(dd_homty_jacmat);
    DoblDobl_Complex_Laur_JacoMats.Clear(dd_homty_mulfac);
    DoblDobl_Complex_Solutions.Clear(dd_start_sols);
    DoblDobl_Complex_Solutions.Clear(dd_target_sols);
  end Clear_DoblDobl_Data;

  procedure Clear_QuadDobl_Data is
  begin
    Clear(qd_homty_expvec);
    QuadDobl_Complex_Poly_Systems.Clear(qd_rndcffsys);
    QuadDobl_Complex_Laur_Systems.Clear(qd_lq);
    QuadDobl_Complex_Laur_SysFun.Clear(qd_hq);
    QuadDobl_Complex_VecVecs.Deep_Clear(qd_homty_coeffv);
    QuadDobl_Complex_Laur_JacoMats.Clear(qd_homty_jacmat);
    QuadDobl_Complex_Laur_JacoMats.Clear(qd_homty_mulfac);
    QuadDobl_Complex_Solutions.Clear(qd_start_sols);
    QuadDobl_Complex_Solutions.Clear(qd_target_sols);
  end Clear_QuadDobl_Data;

  procedure Clear is
  begin
    Clear_Cell_Data;
    Clear_Standard_Data;
    Clear_DoblDobl_Data;
    Clear_QuadDobl_Data;
  end Clear;

end Integer_Cells_Container;
