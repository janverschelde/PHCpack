-- for exception handlers only :
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Arrays_of_Floating_Vector_Lists_io; use Arrays_of_Floating_Vector_Lists_io;
-- normal dependencies :
with text_io;                            use text_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Lists_of_Floating_Vectors;          use Lists_of_Floating_Vectors;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Laur_Functions;    use Standard_Complex_Laur_Functions;
with Exponent_Vectors;                   use Exponent_Vectors;
with Random_Coefficient_Systems;
with Supports_of_Polynomial_Systems;     use Supports_of_Polynomial_Systems;
with Standard_Complex_Laur_Systems;      use Standard_Complex_Laur_Systems;
with Standard_Complex_Laur_SysFun;       use Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;     use Standard_Complex_Laur_Jacomats;
with Standard_Poly_Laur_Convertors;      use Standard_Poly_Laur_Convertors;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Standard_Simpomial_Solvers;
with Floating_Polyhedral_Continuation;   use Floating_Polyhedral_Continuation;
with PHCpack_Operations;

package body Cells_Container is

-- DATA STRUCTURES :

  mix : Standard_Integer_Vectors.Link_to_Vector;
  lifsup,lifsup_last : Link_to_Array_of_Lists;
  cells,cells_last : Mixed_Subdivision;
  rndcffsys : Link_to_Poly_Sys;
  lq : Standard_Complex_Laur_Systems.Link_to_Laur_Sys;
  hq : Standard_Complex_Laur_SysFun.Link_to_Eval_Coeff_Laur_Sys;
  homty_expvec : Link_to_Exponent_Vectors_Array;
  homty_coeffv : Standard_Complex_VecVecs.Link_to_VecVec;
  homty_jacmat : Link_to_Eval_Coeff_Jaco_Mat;
  homty_mulfac : Link_to_Mult_Factors;
  start_sols,target_sols,target_last : Link_to_Array_of_Solution_Lists;

-- AUXILIARIES :

  function Retrieve ( L : List; k : natural32 )
                    return Standard_Floating_Vectors.Link_to_Vector is

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

  function Position ( L : List; x : Standard_Floating_Vectors.Vector )
                    return natural32 is

  -- DESCRIPTION :
  --   Returns 0 if x does not occur in the list L,
  --   otherwise returns the position of x in the list L.

    res : natural32 := 0;
    tmp : List := L;
    p : Standard_Floating_Vectors.Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      res := res + 1;
      p := Head_Of(tmp);
      if Standard_Floating_Vectors.Equal(x,p.all)
       then return res;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return 0;
  end Position;

-- CREATORS :

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

  procedure Generate_Random_Coefficient_System is

    n : constant natural32 := Cells_Container.Dimension-1;
    q : constant Poly_Sys(1..integer32(n))
      := Random_Coefficient_Systems.Create(n,mix.all,lifsup.all);

  begin
    rndcffsys := new Poly_Sys'(q);
  end Generate_Random_Coefficient_System;

  procedure Initialize_Random_Coefficient_System ( q : in Poly_Sys ) is
  begin
    rndcffsys := new Poly_Sys(q'range);
    for i in q'range loop
      Copy(q(i),rndcffsys(i));
    end loop;
  end Initialize_Random_Coefficient_System;

  procedure Create_Polyhedral_Homotopy is

    q : Link_to_Poly_Sys renames rndcffsys;
    use Standard_Complex_Vectors,Standard_Complex_VecVecs;
    len : constant integer32 := integer32(Cells_Container.Length);

  begin
    if lq /= null
     then Clear(lq);
    end if;
    lq := new Laur_Sys'(Polynomial_to_Laurent_System(q.all));
    if hq /= null
     then Clear(hq);
    end if;
    hq := new Eval_Coeff_Laur_Sys'(Create(lq.all));
    if homty_expvec /= null
     then Clear(homty_expvec);
    end if;
    homty_expvec := new Exponent_Vectors_Array'(Create(q.all));
    if homty_coeffv /= null then
      for i in homty_coeffv'range loop
        if homty_coeffv(i) /= null
         then Clear(homty_coeffv(i));
        end if;
      end loop;
    end if;
    homty_coeffv := new Standard_Complex_VecVecs.VecVec(q'range);
    for i in q'range loop
      declare
        c : constant Standard_Complex_Vectors.Vector := Coeff(lq(i));
      begin
        homty_coeffv(i) := new Standard_Complex_Vectors.Vector(c'range);
        for j in c'range loop
          homty_coeffv(i)(j) := c(j);
        end loop;
      end;
    end loop;
    if homty_jacmat /= null
     then Clear(homty_jacmat);
    end if;
    homty_jacmat := new Eval_Coeff_Jaco_Mat(q'range,q'first..q'last+1);
    if homty_mulfac /= null
     then Clear(homty_mulfac);
    end if;
    homty_mulfac
      := new Mult_Factors(homty_jacmat'range(1),homty_jacmat'range(2));
    Create(lq.all,homty_jacmat.all,homty_mulfac.all);
    if start_sols /= null
     then Clear(start_sols);
    end if;
    start_sols := new Array_of_Solution_Lists(1..len);
    if target_sols /= null
     then Clear(target_sols);
    end if;
    target_sols := new Array_of_Solution_Lists(1..len);
    target_last := new Array_of_Solution_Lists(1..len);
  end Create_Polyhedral_Homotopy;

-- SELECTORS :

  function Length return natural32 is
  begin
    return Length_Of(cells);
  end Length;

  function Dimension return natural32 is

    use Standard_Floating_Vectors;

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
               normal : out Standard_Floating_Vectors.Link_to_Vector ) is

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

  function Retrieve_Random_Coefficient_System return Poly_Sys is
  begin
    return rndcffsys.all;
  end Retrieve_Random_Coefficient_System;

  function Retrieve_Random_Coefficient_System return Link_to_Poly_Sys is
  begin
    return rndcffsys;
  end Retrieve_Random_Coefficient_System;

  function Retrieve_Start_Solution
             ( k,i : natural32 ) return Link_to_Solution is
  begin
    if Is_Null(start_sols(integer32(k)))
     then return null;
     else return Standard_Complex_Solutions.Retrieve
                   (start_sols(integer32(k)),i);
    end if;
  end Retrieve_Start_Solution;

  function Retrieve_Target_Solution
             ( k,i : natural32 ) return Link_to_Solution is
  begin
    if Is_Null(target_sols(integer32(k)))
     then return null;
     else return Standard_Complex_Solutions.Retrieve
                   (target_sols(integer32(k)),i);
    end if;
  end Retrieve_Target_Solution;

-- CONSTRUCTOR :

  function Append_to_Support
             ( k : in natural32;
               x : in Standard_Floating_Vectors.Vector ) 
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
               normal : in Standard_Floating_Vectors.Vector ) is

    mic : Mixed_Cell;
    last : Array_of_Lists(cnt'range);
    ind : integer32;
    lpt : Standard_Floating_Vectors.Link_to_Vector;

  begin
    mic.nor := new Standard_Floating_Vectors.Vector'(normal);
    mic.pts := new Array_of_Lists(cnt'range);
    ind := lab'first;
    for i in cnt'range loop
      for j in 1..cnt(i) loop
        lpt := Retrieve(lifsup(i),natural32(lab(ind)));
        ind := ind + 1;
        Append(mic.pts(i),last(i),lpt);
      end loop;
    end loop;
    Cells_Container.Append(mic);
  end Append_Mixed_Cell;

  procedure Solve_Start_System ( k : in natural32; mv : out natural32 ) is

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
        q : Link_to_Poly_Sys renames rndcffsys;
        sub_q : Poly_Sys(q'range);
        lau_q : Laur_Sys(q'range);
        sols : Solution_List;
      begin
        sub_q := Select_Terms(q.all,mix.all,mic.pts.all);
        lau_q := Polynomial_to_Laurent_System(sub_q);
        Shift(lau_q);
       -- put_line("The supported start system : "); put_line(sub_q);
        Standard_Simpomial_Solvers.Solve(lau_q,tol_zero,sols,fail,zero_y);
        if fail then
          mv := 0;
         -- put_line("The system is not a fewnomial system!");
        else
          mv := Length_Of(sols);
          start_sols(integer32(k)) := sols;
        end if;
        Clear(sub_q); Clear(lau_q);
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
  end Solve_Start_System;

  procedure Track_Solution_Path ( k,i,otp : in natural32 ) is

    fail : boolean;
    mic : Mixed_Cell;
    ls : constant Link_to_Solution := Retrieve_Start_Solution(k,i);
    ns : constant Link_to_Solution := new Solution'(ls.all);
    s : Solution_List;

  begin
    Retrieve(k,mic,fail);
    Construct(ns,s);
    if otp = 0 then
      Mixed_Continuation
        (mix.all,lifsup.all,hq.all,homty_coeffv.all,homty_expvec.all,
         homty_jacmat.all,homty_mulfac.all,mic.nor.all,s);
    else
      if PHCpack_Operations.Is_File_Defined then
        Mixed_Continuation(PHCpack_Operations.Output_File,
           mix.all,lifsup.all,hq.all,homty_coeffv.all,homty_expvec.all,
           homty_jacmat.all,homty_mulfac.all,mic.nor.all,s);
      else
        Mixed_Continuation(standard_output,
           mix.all,lifsup.all,hq.all,homty_coeffv.all,homty_expvec.all,
           homty_jacmat.all,homty_mulfac.all,mic.nor.all,s);
      end if;
    end if;
    Append(target_sols(integer32(k)),target_last(integer32(k)),Head_Of(s).all);
    Clear(s);
  end Track_Solution_Path;

-- DESTRUCTOR :

  procedure Clear is
  begin
    Deep_Clear(lifsup);
    Deep_Clear(cells);
    Standard_Integer_Vectors.Clear(mix);
    Clear(rndcffsys);
    Standard_Complex_Laur_Systems.Clear(lq);
    Standard_Complex_Laur_SysFun.Clear(hq);
    Clear(homty_expvec);
    Standard_Complex_VecVecs.Deep_Clear(homty_coeffv);
    Clear(homty_jacmat);
    Clear(homty_mulfac);
    Clear(start_sols); Clear(target_sols);
  end Clear;

end Cells_Container;
