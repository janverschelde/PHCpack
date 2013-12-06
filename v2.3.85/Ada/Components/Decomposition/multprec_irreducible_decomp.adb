with unchecked_deallocation;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Integer_VecVecs;
with Standard_Random_Vectors;            use Standard_Random_Vectors;
with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;       use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;        use Multprec_Complex_Vectors_io;
with Multprec_Complex_Vector_Tools;      use Multprec_Complex_Vector_Tools;
with Multprec_Complex_Polynomials;       use Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_SysFun;       use Multprec_Complex_Poly_SysFun;
with Multprec_Linear_Projections;        use Multprec_Linear_Projections;
with Multprec_Polynomial_Interpolators;  use Multprec_Polynomial_Interpolators;
with Multprec_Breakup_Components;        use Multprec_Breakup_Components;
with Filtered_Points;

package body Multprec_Irreducible_Decomp is

-- AUXILIARY :

  procedure Count_Components
              ( file : in file_type;
                p : in out Multprec_Complex_Poly_Systems.Link_to_Poly_Sys;
                tol : in double_float ) is

    dis : Floating_Number;
    found : array(p'range) of boolean := (p'range => false);
    cnt : integer32;

  begin
    for i in p'range loop                        -- compute mutual distances
      if not found(i) then
        for j in i+1..p'last loop
          if not found(j) then
            if Degree(p(i)) = Degree(p(j)) then
              dis := Distance(p(i),p(j));
              put(file,"Distance between p("); put(file,i,1);
              put(file,") and p("); put(file,j,1);
              put(file,") is "); put(file,dis,3); new_line(file);
              if dis < tol
               then found(j) := true;
              end if;
            end if;
          end if;
        end loop;
      end if;
    end loop;
    cnt := 0;                               -- count the number of components
    for i in found'range loop
      if not found(i)
       then cnt := cnt+1;
      end if;
    end loop;
    put(file,"Found "); put(file,cnt,1);
    put(file," component");
    if cnt = 1
     then put_line(file,".");
     else put_line(file,"s.");
    end if;
    declare                                         -- collect the equations
      equ : Multprec_Complex_Poly_Systems.Poly_Sys(1..cnt);
      ind : integer32 := 0;
    begin
      for i in found'range loop
        if not found(i) then
          ind := ind+1;
          Copy(p(i),equ(ind));
        end if;
      end loop;
      Multprec_Complex_Poly_Systems.Clear(p);
      p := new Multprec_Complex_Poly_Systems.Poly_Sys'(equ);
    end;
  end Count_Components;

-- CREATOR :

  procedure Create_Hypersurfaces
              ( file : in file_type; n,d,size,itp : in natural32;
                skewproj : in boolean;
                embsys : in Standard_Complex_Poly_Systems.Poly_Sys;
                orgsys : in Multprec_Complex_Poly_Systems.Poly_Sys;
                genpts : in Standard_Complex_Solutions.Solution_List;
                surf : out Link_to_Hypersurfaces; timing : out Duration ) is

  -- NOTE :
  --   The random projection operators are standard double floats to make
  --   the evaluation of standard number solutions efficient.

    res : Hypersurfaces(integer32(n),integer32(d));
    tol : double_float := 1.0E-8;
    timer : Timing_Widget;

  begin
    Standard_Complex_Solutions.Copy(genpts,res.pts);
    for i in res.hyp'range loop
      res.hyp(i) := new Multprec_Complex_Vectors.Vector'
                          (Create(Random_Vector(0,integer32(n))));
    end loop;
    tstart(timer);
    case itp is
      when 1 =>
        declare
          equs : constant Multprec_Complex_Poly_Systems.Poly_Sys
               := Massive_Interpolate
                    (file,embsys,orgsys,genpts,res.hyp,d,size);
        begin
          res.equ := new Multprec_Complex_Poly_Systems.Poly_Sys'(equs);
        end;
        Count_Components(file,res.equ,tol);
      when 2 =>
        declare
          equs : constant Multprec_Complex_Poly_Systems.Poly_Sys
               := Incremental_Interpolate
                    (file,embsys,orgsys,genpts,res.hyp,d,size);
        begin
          res.equ := new Multprec_Complex_Poly_Systems.Poly_Sys'(equs);
        end;
        Count_Components(file,res.equ,tol);
      when 3 =>
        declare
          len : constant integer32 
              := integer32(Standard_Complex_Solutions.Length_Of(genpts));
          subspaces : Multprec_Complex_Poly_Systems.Array_of_Poly_Sys(1..len);
          pivots : Standard_Integer_VecVecs.VecVec(1..len);
          basepts : Multprec_Complex_VecVecs.Array_of_VecVecs(1..len);
          basecard : Standard_Natural_Vectors.Vector(1..len) := (1..len => 0);
        begin
          Dynamic_Interpolate
            (file,embsys,orgsys,genpts,d,size,res.hyp,skewproj,
             subspaces,pivots,basepts,basecard,res.equ);
        end;
      when others => null;
    end case;
    tstop(timer);
    surf := new Hypersurfaces'(res);
    timing := Elapsed_User_Time(timer);
  end Create_Hypersurfaces;

-- AUXILIARY TO SELECTORS :

  function Zero_Component ( v : Multprec_Complex_Vectors.Vector;
                            tol : double_float ) return integer32 is

  -- DESCRIPTION :
  --   If |v(i)| < tol, then returns i, otherwise returns v'first-1.
 
  begin
    for i in v'range loop
      if AbsVal(v(i)) < tol
       then return i;
      end if;
    end loop;
    return v'first-1;
  end Zero_Component;

-- SELECTORS :

  procedure Component_Test
              ( file : in file_type; level : in natural32;
                surf : in Hypersurfaces; tol : in double_float;
                sols : in Standard_Complex_Solutions.Solution_List;
                solslvl : in natural32; num : out natural32;
                filpts : in out List;
                rem_first,rem_last
                  : in out Standard_Complex_Solutions.Solution_List ) is

    n : constant natural32 := Number_of_Unknowns(surf.equ(surf.equ'first));
    projsols : constant VecVec := Evaluate(surf.hyp,sols,integer32(n));
    eva : Multprec_Complex_Vectors.Vector(surf.equ'range);
    cnt : natural32 := 0;
    tmp : Standard_Complex_Solutions.Solution_List := sols;
    ls : Standard_Complex_Solutions.Link_to_Solution;
    ind : integer32;

  begin
    put(file,"COMPONENT TEST AT LEVEL "); put(file,level,1);
    put_line(file," :");
    for i in projsols'range loop
      ls := Standard_Complex_Solutions.Head_Of(tmp);
      eva := Eval(surf.equ.all,projsols(i).all);
      put(file," eva("); put(file,i,1); put(file,") : ");
      put(file,eva,3); new_line(file);
      put(file,"Solution "); put(file,i,1);
      ind := Zero_Component(eva,tol);
      if ind >= eva'first then
        put(file," lies on component ");
        put(file,ind,1); put_line(file,".");
        cnt := cnt + 1;
        Filtered_Points.Update(filpts,integer32(level),ind,integer32(solslvl));
      else
        put_line(file," fails the component test.");
        Standard_Complex_Solutions.Append(rem_first,rem_last,ls.all);
      end if;
      Multprec_Complex_Vectors.Clear(eva);
      tmp := Standard_Complex_Solutions.Tail_Of(tmp);
    end loop;
    num := cnt;
  end Component_Test;

  procedure Component_Test
              ( file : in file_type; soco : in Solution_Components;
                tol : in double_float;
                sols : in Standard_Complex_Solutions.Solution_List;
                level : in natural32; num : out natural32;
                filpts : in out List;
                rem_sols : out Standard_Complex_Solutions.Solution_List ) is

    cnt,icnt : natural32 := 0;
    wrk,rem_first,rem_last : Standard_Complex_Solutions.Solution_List;
    first_time : boolean := true;

  begin
    Standard_Complex_Solutions.Copy(sols,wrk);
    for i in reverse 1..soco'last loop
      if soco(i) /= null then
        if first_time then
          first_time := false;
        else
          Standard_Complex_Solutions.Copy(rem_first,wrk);
          Standard_Complex_Solutions.Clear(rem_first);
          rem_last := rem_first;
        end if;
        Component_Test(file,natural32(i),soco(i).all,tol,wrk,level,
                       icnt,filpts,rem_first,rem_last);
        cnt := cnt + icnt;
      end if;
    end loop;
    num := cnt;
    rem_sols := rem_first;
  end Component_Test;

-- DESTRUCTORS :

  procedure Clear ( h : in out Hypersurfaces ) is
  begin
    Standard_Complex_Solutions.Clear(h.pts);
    Multprec_Complex_Poly_Systems.Clear(h.equ);
  end Clear;

  procedure Clear ( h : in out Link_to_Hypersurfaces ) is

    procedure free is
      new unchecked_deallocation(Hypersurfaces,Link_to_Hypersurfaces);

  begin
    if h /= null
     then Clear(h.all); free(h);
    end if;
  end Clear;

  procedure Clear ( s : in out Solution_Components ) is
  begin
    for i in s'range loop
      Clear(s(i));
    end loop;
  end Clear;

end Multprec_Irreducible_Decomp;
