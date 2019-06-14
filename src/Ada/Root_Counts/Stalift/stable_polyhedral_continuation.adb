with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Integer_Vectors_io;
with Standard_Complex_VecVecs;
with DoblDobl_Complex_VecVecs;
with QuadDobl_Complex_VecVecs;
with Standard_Complex_Laur_Functions;
with Standard_Complex_Laur_SysFun;
with Standard_Complex_Laur_Jacomats;
with DoblDobl_Complex_Laur_Functions;
with DoblDobl_Complex_Laur_SysFun;
with DoblDobl_Complex_Laur_Jacomats;
with QuadDobl_Complex_Laur_Functions;
with QuadDobl_Complex_Laur_SysFun;
with QuadDobl_Complex_Laur_Jacomats;
with Exponent_Vectors;                  use Exponent_Vectors;
with Standard_Complex_Laur_Systems_io;  use Standard_Complex_Laur_Systems_io;
with DoblDobl_Complex_Laur_Systems_io;  use DoblDobl_Complex_Laur_Systems_io;
with QuadDobl_Complex_Laur_Systems_io;  use QuadDobl_Complex_Laur_Systems_io;
with Standard_Stable_Homotopies;
with DoblDobl_Stable_Homotopies;
with QuadDobl_Stable_Homotopies;
with Floating_Mixed_Subdivisions_io;    use Floating_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;          use Mixed_Volume_Computation;
with Floating_Polyhedral_Continuation;  use Floating_Polyhedral_Continuation;
with DoblDobl_Polyhedral_Continuation;  use DoblDobl_Polyhedral_Continuation;
with QuadDobl_Polyhedral_Continuation;  use QuadDobl_Polyhedral_Continuation;

package body Stable_Polyhedral_Continuation is

-- PREPARATORY FUNCTIONS :

  function Eliminate_Zeroes
             ( v : Standard_Floating_Vectors.Vector;
               z : Standard_Integer_Vectors.Vector; nbz : integer32 )
             return Standard_Floating_Vectors.Vector is
  begin
    if nbz <= 0 then
      return v;
    else
      declare
        res : Standard_Floating_Vectors.Vector(v'first..v'last-nbz);
        ind : integer32 := v'first-1;
      begin
        for i in z'range loop
          if z(i) > 0 then
            ind := ind + 1;
            res(ind) := v(i);
          end if;
        end loop;
        if v'last > z'last then
          res(res'last) := v(v'last); -- copy lifting value
        end if;
        return res;
      end;
    end if;
  end Eliminate_Zeroes;

  function Vanish_by_Zeroes
             ( x : Standard_Floating_Vectors.Vector;
               z : Standard_Integer_Vectors.Vector; nbz : integer32 )
             return boolean is
  begin
    if nbz <= 0 then
      return false;
    else
      for i in z'range loop
        if ((z(i) = 0) and then (x(i) /= 0.0))
         then return true;
        end if;
      end loop;
      return false;
    end if;
  end Vanish_by_Zeroes;

  function Substitute_Zeroes
             ( pts : Lists_of_Floating_Vectors.List; nbz : integer32;
               z : Standard_Integer_Vectors.Vector )
             return Lists_of_Floating_Vectors.List is

    use Lists_of_Floating_Vectors;

  begin
    if nbz <= 0 then
      return pts;
    else
      declare
        res,res_last,tmp : List;
        lpt : Standard_Floating_Vectors.Link_to_Vector;
      begin
        tmp := pts;
        while not Is_Null(tmp) loop
          lpt := Head_Of(tmp);
          if not Vanish_by_Zeroes(lpt.all,z,nbz) then
            Append(res,res_last,Eliminate_Zeroes(lpt.all,z,nbz));
          end if;
          tmp := Tail_Of(tmp);
        end loop;
        return res;
      end;
    end if;
  end Substitute_Zeroes;

  function Substitute_Zeroes
             ( pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
               nbz : integer32; z : Standard_Integer_Vectors.Vector )
             return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Floating_Vector_Lists.Array_of_Lists(pts'range);

  begin
    for i in pts'range loop
      res(i) := Substitute_Zeroes(pts(i),nbz,z);
    end loop;
    return res;
  end Substitute_Zeroes;

  function Substitute_Zeroes
             ( mic : Mixed_Cell; nbz : integer32;
               z : Standard_Integer_Vectors.Vector ) return Mixed_Cell is

    use Arrays_of_Floating_Vector_Lists;

  begin
    if nbz <= 0 then
      return mic;
    else
      declare
        res : Mixed_Cell;
      begin
        res.nor := new Standard_Floating_Vectors.Vector'
                         (Eliminate_Zeroes(mic.nor.all,z,nbz));
        res.pts := new Array_of_Lists(mic.pts'range);
        for i in mic.pts'range loop
          res.pts(i) := Substitute_Zeroes(mic.pts(i),nbz,z);
        end loop;
        return res;
      end;
    end if;
  end Substitute_Zeroes;

  function Filter ( pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
                  return Arrays_of_Floating_Vector_Lists.Array_of_Lists is

    use Arrays_of_Floating_Vector_Lists;
    use Lists_of_Floating_Vectors;

    res : Array_of_Lists(pts'range);
    cnt : integer32 := pts'first - 1;

  begin
    for i in pts'range loop
      if Length_Of(pts(i)) > 1 then
        cnt := cnt + 1;
        res(cnt) := pts(i);
      end if;
    end loop;
    if cnt >= pts'first
     then return res(res'first..cnt);
     else return res(res'first..res'first);
    end if;
  end Filter;

  function Filter ( pts : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                    mix : Standard_Integer_Vectors.Link_to_Vector )
                  return Standard_Integer_Vectors.Link_to_Vector is

    res : Standard_Integer_Vectors.Link_to_Vector;
    newmix : Standard_Integer_Vectors.Vector(mix'range);
    cnt : integer32 := mix'first - 1;

  begin
    for i in pts'range loop
      if Lists_of_Floating_Vectors.Length_Of(pts(i)) > 1 then
        cnt := cnt + 1;
        newmix(cnt) := mix(i);
      end if;
    end loop;
    if cnt >= newmix'first then
      res := new Standard_Integer_Vectors.Vector'(newmix(newmix'first..cnt));
    else
      res := new Standard_Integer_Vectors.Vector'
                    (newmix(newmix'first..newmix'first));
    end if;
    return res;
  end Filter;

  function Filter ( mic : Mixed_Cell ) return Mixed_Cell is

    use Arrays_of_Floating_Vector_Lists;

    res : Mixed_Cell;

  begin
    res.nor := mic.nor;
    res.pts := new Array_of_Lists'(Filter(mic.pts.all));
    res.sub := mic.sub;
    return res;
  end Filter;

-- DEFINING POLYHEDRAL HOMOTOPIES :

  procedure Assign_Multiplicity
               ( s : in out Standard_Complex_Solutions.Solution_List;
                 vol : in natural32 ) is

    use Standard_Complex_Solutions;
    m : integer32;
    tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    if vol > Length_Of(s) then
      m := integer32(vol/Length_Of(s));
      tmp := s;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        ls.m := m;
        Set_Head(tmp,ls);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Assign_Multiplicity;

  procedure Assign_Multiplicity
               ( s : in out DoblDobl_Complex_Solutions.Solution_List;
                 vol : in natural32 ) is

    use DoblDobl_Complex_Solutions;
    m : integer32;
    tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    if vol > Length_Of(s) then
      m := integer32(vol/Length_Of(s));
      tmp := s;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        ls.m := m;
        Set_Head(tmp,ls);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Assign_Multiplicity;

  procedure Assign_Multiplicity
               ( s : in out QuadDobl_Complex_Solutions.Solution_List;
                 vol : in natural32 ) is

    use QuadDobl_Complex_Solutions;
    m : integer32;
    tmp : Solution_List;
    ls : Link_to_Solution;

  begin
    if vol > Length_Of(s) then
      m := integer32(vol/Length_Of(s));
      tmp := s;
      while not Is_Null(tmp) loop
        ls := Head_Of(tmp);
        ls.m := m;
        Set_Head(tmp,ls);
        tmp := Tail_Of(tmp);
      end loop;
    end if;
  end Assign_Multiplicity;

  function Create ( q : in Standard_Complex_Laur_Systems.Laur_Sys )
                  return Standard_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficient vectors of the polynomials in q.

    res : Standard_Complex_VecVecs.VecVec(q'range);

  begin
    for i in q'range loop
      declare
        coeff_q : constant Standard_Complex_Vectors.Vector
                := Standard_Complex_Laur_Functions.Coeff(q(i));
      begin
        res(i) := new Standard_Complex_Vectors.Vector(coeff_q'range);
        for k in coeff_q'range loop
          res(i)(k) := coeff_q(k);
        end loop;
      end;
    end loop;
    return res;
  end Create;

  function Create ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys )
                  return DoblDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficient vectors of the polynomials in q.

    res : DoblDobl_Complex_VecVecs.VecVec(q'range);

  begin
    for i in q'range loop
      declare
        coeff_q : constant DoblDobl_Complex_Vectors.Vector
                := DoblDobl_Complex_Laur_Functions.Coeff(q(i));
      begin
        res(i) := new DoblDobl_Complex_Vectors.Vector(coeff_q'range);
        for k in coeff_q'range loop
          res(i)(k) := coeff_q(k);
        end loop;
      end;
    end loop;
    return res;
  end Create;

  function Create ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys )
                  return QuadDobl_Complex_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the coefficient vectors of the polynomials in q.

    res : QuadDobl_Complex_VecVecs.VecVec(q'range);

  begin
    for i in q'range loop
      declare
        coeff_q : constant QuadDobl_Complex_Vectors.Vector
                := QuadDobl_Complex_Laur_Functions.Coeff(q(i));
      begin
        res(i) := new QuadDobl_Complex_Vectors.Vector(coeff_q'range);
        for k in coeff_q'range loop
          res(i)(k) := coeff_q(k);
        end loop;
      end;
    end loop;
    return res;
  end Create;

  function Number_of_Zero_Polynomials
             ( q : Standard_Complex_Laur_Systems.Laur_Sys )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the number of zero polynomials in q.

    use Standard_Complex_Laurentials;

    res : integer32 := 0;

  begin
    for i in q'range loop
      if q(i) = Null_Poly
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Zero_Polynomials;

  function Number_of_Zero_Polynomials
             ( q : DoblDobl_Complex_Laur_Systems.Laur_Sys )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the number of zero polynomials in q.

    use DoblDobl_Complex_Laurentials;

    res : integer32 := 0;

  begin
    for i in q'range loop
      if q(i) = Null_Poly
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Zero_Polynomials;

  function Number_of_Zero_Polynomials
             ( q : QuadDobl_Complex_Laur_Systems.Laur_Sys )
             return integer32 is

  -- DESCRIPTION :
  --   Returns the number of zero polynomials in q.

    use QuadDobl_Complex_Laurentials;

    res : integer32 := 0;

  begin
    for i in q'range loop
      if q(i) = Null_Poly
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Number_of_Zero_Polynomials;

  -- function Number_of_Invalid_Supports
  --            ( s : Arrays_of_Floating_Vector_Lists.Array_of_Lists )
  --            return integer32 is

  -- DESCRIPTION :
  --   Returns the number of supports in s with less than two points.

  --  res : integer32 := 0;
  
  -- begin
  --   for i in s'range loop
  --     if Lists_of_Floating_Vectors.Length_Of(s(i)) < 2
  --      then res := res + 1;
  --     end if;
  --   end loop;
  --   return res;
  -- end Number_of_Invalid_Supports;

  function Number_of_Invalid_Supports
              ( s : Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mix : Standard_Integer_Vectors.Link_to_Vector )
              return integer32 is

  -- DESCRIPTION :
  --   Returns the number of supports in s with less than two points.
  --   The support are counted with the type of mixture,
  --   provided in mix.

    res : integer32 := 0;

  begin
    for i in s'range loop
      if Lists_of_Floating_Vectors.Length_Of(s(i)) < 2
       then res := res + mix(i);
      end if;
    end loop;
    return res;
  end Number_of_Invalid_Supports;

  procedure Silent_Polyhedral_Continuation
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Arrays_of_Floating_Vector_Lists,Standard_Complex_Solutions;
    use Standard_Complex_Laur_Systems,Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_Jacomats,Standard_Stable_Homotopies;

   -- n : constant natural32 := natural32(q'last);
    sq : Laur_Sys(q'range) := Substitute_Zeroes(q,ztp,nbz);

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 1,");
      put_line("for Laurent systems in double precision ...");
    end if;
   -- if Number_of_Zero_Polynomials(sq) /= nbz then
   --   put_line("Degenerate subsystem, there are no isolated solutions.");
   -- else
    if Number_of_Zero_Polynomials(sq) = nbz then
     -- put_line("Nondegenerate subsystem, there are isolated solutions.");
     -- put("nbz = "); put(nbz,1); new_line;
      declare
        fq : constant Laur_Sys(q'first..q'last-nbz) := Filter(sq);
        smic : constant Mixed_Cell := Substitute_Zeroes(mic,nbz,ztp);
        fmic : constant Mixed_Cell := Filter(smic);
        slif : constant Array_of_Lists(lif'range)
             := Substitute_Zeroes(lif,nbz,ztp);
        flif : constant Array_of_Lists := Filter(slif);
        fmix : Standard_Integer_Vectors.Link_to_Vector;
        h : Eval_Coeff_Laur_Sys(fq'range);
        c : Standard_Complex_VecVecs.VecVec(h'range);
        e : Exponent_Vectors_Array(h'range);
        j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
        m : Mult_Factors(j'range(1),j'range(2));
        sols1,sols1_last : Solution_List;
      begin
       -- put_line("The original lifted points : "); put(lif);
       -- put_line("After substituting zeroes : "); put(slif);
        if Number_of_Invalid_Supports(slif,mix) = nbz then
          fmix := Filter(slif,mix);
         -- put_line("The lifted points after substitution : "); put(flif);
         -- put("The new type of mixture : ");
         -- Standard_Integer_Vectors_io.put(fmix.all); new_line;
         -- put_line("The system after substituting zeroes : "); put_line(fq);
         -- put_line("The cell after elimination :"); put(n,mix.all,fmic);
          if Number_of_Invalid_Supports(flif,fmix) = 0
            and then Number_of_Zero_Polynomials(fq) = 0 then
              h := Create(fq); c := Create(fq); e := Create(fq);
              Create(fq,j,m);
              Mixed_Solve(fq,flif,h,c,e,j,m,fmix.all,fmic,sols1,sols1_last);
              sols := Insert_Zeroes(sols1,ztp);
              if not Is_Null(sols)
               then Assign_Multiplicity(sols,vol);
              end if;
              Deep_Clear(sols1);
              Clear(h); Clear(j); Clear(m); Clear(e);
              Standard_Complex_VecVecs.Clear(c);
          end if;
        end if;
      exception 
        when others => put_line("exception in declare block...");
                       put("nbz = "); put(nbz,1); new_line;
                       put_line("The system after filtering : "); put(fq);
                       raise;
      end;
    end if;
    Clear(sq);
  exception
    when others => put_line("exception in silent polyhedral continuation");
                   raise;
  end Silent_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Arrays_of_Floating_Vector_Lists,DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Laur_Systems,DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_Jacomats,DoblDobl_Stable_Homotopies;

   -- n : constant natural := q'last;
    sq : Laur_Sys(q'range) := Substitute_Zeroes(q,ztp,nbz);

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 2,");
      put_line("for Laurent systems in double double precision ...");
    end if;
   -- if Number_of_Zero_Polynomials(sq) /= nbz then
   --   put_line("Degenerate subsystem, there are no isolated solutions.");
   -- else
    if Number_of_Zero_Polynomials(sq) = nbz then
     -- put_line("Nondegenerate subsystem, there are isolated solutions.");
      declare
        fq : constant Laur_Sys(q'first..q'last-nbz) := Filter(sq);
        smic : constant Mixed_Cell := Substitute_Zeroes(mic,nbz,ztp);
        fmic : constant Mixed_Cell := Filter(smic);
        slif : constant Array_of_Lists(lif'range)
             := Substitute_Zeroes(lif,nbz,ztp);
        flif : constant Array_of_Lists := Filter(slif);
        fmix : Standard_Integer_Vectors.Link_to_Vector;
        h : Eval_Coeff_Laur_Sys(fq'range);
        c : DoblDobl_Complex_VecVecs.VecVec(h'range);
        e : Exponent_Vectors_Array(h'range);
        j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
        m : Mult_Factors(j'range(1),j'range(2));
        sols1,sols1_last : Solution_List;
      begin
       -- put_line("The original lifted points : "); put(lif);
       -- put_line("After substituting zeroes : "); put(slif);
       -- if Number_of_Invalid_Supports(slif) = nbz then
        if Number_of_Invalid_Supports(slif,mix) = nbz then
          fmix := Filter(slif,mix);
         -- put_line("The lifted points after substitution : "); put(flif);
         -- put("The new type of mixture : ");
         -- Standard_Integer_Vectors_io.put(fmix.all); new_line;
         -- put_line("The system after substituting zeroes : "); put_line(fq);
         -- put_line("The cell after elimination :"); put(n,mix.all,fmic);
         -- if Number_of_Invalid_Supports(flif) = 0 and then
          if Number_of_Invalid_Supports(flif,fmix) = 0 and then
            Number_of_Zero_Polynomials(fq) = 0 then
            h := Create(fq); c := Create(fq); e := Create(fq);
            Create(fq,j,m);
            Mixed_Solve(fq,flif,h,c,e,j,m,fmix.all,fmic,sols1,sols1_last);
            sols := Insert_Zeroes(sols1,ztp);
            if not Is_Null(sols)
             then Assign_Multiplicity(sols,vol);
            end if;
            Deep_Clear(sols1);
            Clear(h); Clear(j); Clear(m); Clear(e);
            DoblDobl_Complex_VecVecs.Clear(c);
          end if;
        end if;
      exception 
        when others => put_line("exception in declare block..."); raise;
      end;
    end if;
    Clear(sq);
  exception
    when others => put_line("exception in silent polyhedral continuation");
                   raise;
  end Silent_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Arrays_of_Floating_Vector_Lists,QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Laur_Systems,QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_Jacomats,QuadDobl_Stable_Homotopies;

   -- n : constant natural := q'last;
    sq : Laur_Sys(q'range) := Substitute_Zeroes(q,ztp,nbz);

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 3,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
   -- if Number_of_Zero_Polynomials(sq) /= nbz then
   --   put_line("Degenerate subsystem, there are no isolated solutions.");
   -- else
    if Number_of_Zero_Polynomials(sq) = nbz then
     -- put_line("Nondegenerate subsystem, there are isolated solutions.");
      declare
        fq : constant Laur_Sys(q'first..q'last-nbz) := Filter(sq);
        smic : constant Mixed_Cell := Substitute_Zeroes(mic,nbz,ztp);
        fmic : constant Mixed_Cell := Filter(smic);
        slif : constant Array_of_Lists(lif'range)
             := Substitute_Zeroes(lif,nbz,ztp);
        flif : constant Array_of_Lists := Filter(slif);
        fmix : Standard_Integer_Vectors.Link_to_Vector;
        h : Eval_Coeff_Laur_Sys(fq'range);
        c : QuadDobl_Complex_VecVecs.VecVec(h'range);
        e : Exponent_Vectors_Array(h'range);
        j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
        m : Mult_Factors(j'range(1),j'range(2));
        sols1,sols1_last : Solution_List;
      begin
       -- put_line("The original lifted points : "); put(lif);
       -- put_line("After substituting zeroes : "); put(slif);
       -- if Number_of_Invalid_Supports(slif) = nbz then
        if Number_of_Invalid_Supports(slif,mix) = nbz then
          fmix := Filter(slif,mix);
         -- put_line("The lifted points after substitution : "); put(flif);
         -- put("The new type of mixture : ");
         -- Standard_Integer_Vectors_io.put(fmix.all); new_line;
         -- put_line("The system after substituting zeroes : "); put_line(fq);
         -- put_line("The cell after elimination :"); put(n,mix.all,fmic);
         -- if Number_of_Invalid_Supports(flif) = 0 and then
          if Number_of_Invalid_Supports(flif,fmix) = 0 and then
            Number_of_Zero_Polynomials(fq) = 0 then
            h := Create(fq); c := Create(fq); e := Create(fq);
            Create(fq,j,m);
            Mixed_Solve(fq,flif,h,c,e,j,m,fmix.all,fmic,sols1,sols1_last);
            sols := Insert_Zeroes(sols1,ztp);
            if not Is_Null(sols)
             then Assign_Multiplicity(sols,vol);
            end if;
            Deep_Clear(sols1);
            Clear(h); Clear(j); Clear(m); Clear(e);
            QuadDobl_Complex_VecVecs.Clear(c);
          end if;
        end if;
      exception 
        when others => put_line("exception in declare block..."); raise;
      end;
    end if;
    Clear(sq);
  exception
    when others => put_line("exception in silent polyhedral continuation");
                   raise;
  end Silent_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Arrays_of_Floating_Vector_Lists,Standard_Complex_Solutions;
    use Standard_Complex_Laur_Systems,Standard_Complex_Laur_SysFun;
    use Standard_Complex_Laur_Jacomats,Standard_Stable_Homotopies;

    n : constant integer32 := q'last;
    sq : Laur_Sys(q'range) := Substitute_Zeroes(q,ztp,nbz);

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 1,");
      put_line("for Laurent systems in double precision ...");
    end if;
    if Number_of_Zero_Polynomials(sq) /= nbz then
      put_line(file,"Degenerate subsystem, there are no isolated solutions.");
      put_line(file,"the Laurent system after substitution : ");
      put(file,sq);
    else
      put_line(file,"Nondegenerate subsystem, there are isolated solutions.");
      declare
        fq : constant Laur_Sys(q'first..q'last-nbz) := Filter(sq);
        smic : constant Mixed_Cell := Substitute_Zeroes(mic,nbz,ztp);
        fmic : constant Mixed_Cell := Filter(smic);
        slif : constant Array_of_Lists(lif'range)
             := Substitute_Zeroes(lif,nbz,ztp);
        flif : constant Array_of_Lists := Filter(slif);
        fmix : Standard_Integer_Vectors.Link_to_Vector;
        h : Eval_Coeff_Laur_Sys(fq'range);
        c : Standard_Complex_VecVecs.VecVec(h'range);
        e : Exponent_Vectors_Array(h'range);
        j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
        m : Mult_Factors(j'range(1),j'range(2));
        sols1,sols1_last : Solution_List;
        nbzp : integer32;
        nbis : constant integer32 -- := Number_of_Invalid_Supports(slif);
             := Number_of_Invalid_Supports(slif,mix);
      begin
        if nbis /= nbz then
          put(file,"The number of invalid supports : ");
          put(file,nbis,1); put_line(file," degenerate case ...");
        else
          fmix := Filter(slif,mix);
          put_line(file,"The original lifted points : "); put(file,lif);
          put_line(file,"The lifted points after substitution : ");
          put(file,flif);
          put(file,"The new type of mixture : ");
          Standard_Integer_Vectors_io.put(file,fmix.all); new_line(file);
          put_line(file,"The system after substituting zeroes : ");
          put_line(file,fq);
          nbzp := Number_of_Zero_Polynomials(fq);
          if nbzp > 0 then
            put(file,"Number of zero polynomials : "); put(file,nbzp,1);
            put_line(file," degenerate case.");
          else
            h := Create(fq); c := Create(fq); e := Create(fq);
            put_line(file,"The cell after elimination :");
            put(file,natural32(n),mix.all,fmic);
            Create(fq,j,m);
            Mixed_Solve(file,fq,flif,h,c,e,j,m,fmix.all,fmic,sols1,sols1_last);
            sols := Insert_Zeroes(sols1,ztp);
            if not Is_Null(sols)
             then Assign_Multiplicity(sols,vol);
            end if;
            Deep_Clear(sols1);
            Clear(j); Clear(m);
            Clear(h); Clear(e);
            Standard_Complex_VecVecs.Clear(c);
          end if;
        end if;
      end;
    end if;
    Clear(sq);
  exception
    when others =>
      put_line("Exception raised in Reporting_Polyhedral_Continuation,");
      put_line("in the package Stable_Polyhedral_Continuation.");
      raise;
  end Reporting_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Arrays_of_Floating_Vector_Lists,DoblDobl_Complex_Solutions;
    use DoblDobl_Complex_Laur_Systems,DoblDobl_Complex_Laur_SysFun;
    use DoblDobl_Complex_Laur_Jacomats,DoblDobl_Stable_Homotopies;

    n : constant integer32 := q'last;
    sq : Laur_Sys(q'range) := Substitute_Zeroes(q,ztp,nbz);

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 2,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    if Number_of_Zero_Polynomials(sq) /= nbz then
      put_line(file,"Degenerate subsystem, there are no isolated solutions.");
      put_line(file,"the Laurent system after substitution : ");
      put(file,sq);
    else
      put_line(file,"Nondegenerate subsystem, there are isolated solutions.");
      declare
        fq : constant Laur_Sys(q'first..q'last-nbz) := Filter(sq);
        smic : constant Mixed_Cell := Substitute_Zeroes(mic,nbz,ztp);
        fmic : constant Mixed_Cell := Filter(smic);
        slif : constant Array_of_Lists(lif'range)
             := Substitute_Zeroes(lif,nbz,ztp);
        flif : constant Array_of_Lists := Filter(slif);
        fmix : Standard_Integer_Vectors.Link_to_Vector;
        h : Eval_Coeff_Laur_Sys(fq'range);
        c : DoblDobl_Complex_VecVecs.VecVec(h'range);
        e : Exponent_Vectors_Array(h'range);
        j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
        m : Mult_Factors(j'range(1),j'range(2));
        sols1,sols1_last : Solution_List;
        nbzp : integer32;
        nbis : constant integer32 -- := Number_of_Invalid_Supports(slif);
             := Number_of_Invalid_Supports(slif,mix);
      begin
        if nbis /= nbz then
          put(file,"nbz = "); put(file,nbz,1); new_line(file);
          put(file,"The number of invalid supports : ");
          put(file,nbis,1); put_line(file," degenerate case ...");
        else
          fmix := Filter(slif,mix);
          put_line(file,"The original lifted points : "); put(file,lif);
          put_line(file,"The lifted points after substitution : ");
          put(file,flif);
          put(file,"The new type of mixture : ");
          Standard_Integer_Vectors_io.put(file,fmix.all); new_line(file);
          put_line(file,"The system after substituting zeroes : ");
          put_line(file,fq);
          nbzp := Number_of_Zero_Polynomials(fq);
          if nbzp > 0 then
            put(file,"Number of zero polynomials : "); put(file,nbzp,1);
            put_line(file," degenerate case.");
          else
            h := Create(fq); c := Create(fq); e := Create(fq);
            put_line(file,"The cell after elimination :");
            put(file,natural32(n),mix.all,fmic);
            Create(fq,j,m);
            Mixed_Solve(file,fq,flif,h,c,e,j,m,fmix.all,fmic,sols1,sols1_last);
            sols := Insert_Zeroes(sols1,ztp);
            if not Is_Null(sols)
             then Assign_Multiplicity(sols,vol);
            end if;
            Deep_Clear(sols1);
            Clear(j); Clear(m);
            Clear(h); Clear(e);
            DoblDobl_Complex_VecVecs.Clear(c);
          end if;
        end if;
      end;
    end if;
    Clear(sq);
  end Reporting_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in Mixed_Cell; nbz : in integer32; vol : in natural32;
                ztp : in Standard_Integer_Vectors.Vector;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Arrays_of_Floating_Vector_Lists,QuadDobl_Complex_Solutions;
    use QuadDobl_Complex_Laur_Systems,QuadDobl_Complex_Laur_SysFun;
    use QuadDobl_Complex_Laur_Jacomats,QuadDobl_Stable_Homotopies;

    n : constant integer32 := q'last;
    sq : Laur_Sys(q'range) := Substitute_Zeroes(q,ztp,nbz);

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 3,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    if Number_of_Zero_Polynomials(sq) /= nbz then
      put_line(file,"Degenerate subsystem, there are no isolated solutions.");
      put_line(file,"the Laurent system after substitution : ");
      put(file,sq);
    else
      put_line(file,"Nondegenerate subsystem, there are isolated solutions.");
      declare
        fq : constant Laur_Sys(q'first..q'last-nbz) := Filter(sq);
        smic : constant Mixed_Cell := Substitute_Zeroes(mic,nbz,ztp);
        fmic : constant Mixed_Cell := Filter(smic);
        slif : constant Array_of_Lists(lif'range)
             := Substitute_Zeroes(lif,nbz,ztp);
        flif : constant Array_of_Lists := Filter(slif);
        fmix : Standard_Integer_Vectors.Link_to_Vector;
        h : Eval_Coeff_Laur_Sys(fq'range);
        c : QuadDobl_Complex_VecVecs.VecVec(h'range);
        e : Exponent_Vectors_Array(h'range);
        j : Eval_Coeff_Jaco_Mat(h'range,h'first..h'last+1);
        m : Mult_Factors(j'range(1),j'range(2));
        sols1,sols1_last : Solution_List;
        nbzp : integer32;
        nbis : constant integer32 -- := Number_of_Invalid_Supports(slif);
             := Number_of_Invalid_Supports(slif,mix);
      begin
        if nbis /= nbz then
          put_line(file,"The number of invalid supports : ");
          put(file,nbis,1); put_line(file," degenerate case...");
        else
          fmix := Filter(slif,mix);
          put_line(file,"The original lifted points : "); put(file,lif);
          put_line(file,"The lifted points after substitution : ");
          put(file,flif);
          put(file,"The new type of mixture : ");
          Standard_Integer_Vectors_io.put(file,fmix.all); new_line(file);
          put_line(file,"The system after substituting zeroes : ");
          put_line(file,fq);
          nbzp := Number_of_Zero_Polynomials(fq);
          if nbzp > 0 then
            put(file,"Number of zero polynomials : "); put(file,nbzp,1);
            put_line(file," degenerate case.");
          else
            h := Create(fq); c := Create(fq); e := Create(fq);
            put_line(file,"The cell after elimination :");
            put(file,natural32(n),mix.all,fmic);
            Create(fq,j,m);
            Mixed_Solve(file,fq,flif,h,c,e,j,m,fmix.all,fmic,sols1,sols1_last);
            sols := Insert_Zeroes(sols1,ztp);
            if not Is_Null(sols)
             then Assign_Multiplicity(sols,vol);
            end if;
            Deep_Clear(sols1);
            Clear(j); Clear(m);
            Clear(h); Clear(e);
            QuadDobl_Complex_VecVecs.Clear(c);
          end if;
        end if;
      end;
    end if;
    Clear(sq);
  end Reporting_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
              ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; -- k : in integer32;
                sols : out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions,Standard_Stable_Homotopies;

    ztp : constant Standard_Integer_Vectors.Vector(q'range)
        := Zero_Type(mic.nor.all,b,mic.pts.all);
    nbz : constant integer32 := Number_of_Zeroes(ztp);
    vol : natural32;
    ls : Link_to_Solution;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 4,");
      put_line("for Laurent systems in double precision ...");
    end if;
   -- put("Type of Cell "); put(k,1); put(" : ");
   -- Standard_Integer_Vectors_io.put(ztp); 
   -- put("  #zeroes : "); put(nbz,1); new_line;
    Mixed_Volume(q'last,mix.all,mic,vol);
    if nbz = q'last then
     -- put("Origin is solution with multiplicity "); put(vol,1); new_line;
      ls := new Solution'(Origin(q'last,integer32(vol)));
      Construct(ls,sols);
    else
     -- put("Solution with "); put(nbz,1); put(" zeroes and volume ");
     -- put(vol,1); new_line;
      Silent_Polyhedral_Continuation
        (q,mix,lif,mic,nbz,vol,ztp,sols,verbose=>verbose-1);
    end if;
  exception
    when others => put_line("exception in main silent polyhedral continuation");
                   raise;
  end Silent_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
              ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; -- k : in integer32;
                sols : out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions,DoblDobl_Stable_Homotopies;

    ztp : constant Standard_Integer_Vectors.Vector(q'range)
        := Zero_Type(mic.nor.all,b,mic.pts.all);
    nbz : constant integer32 := Number_of_Zeroes(ztp);
    vol : natural32;
    ls : Link_to_Solution;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 5,");
      put_line("for Laurent systems in double double precision ...");
    end if;
   -- put("Type of Cell "); put(k,1); put(" : ");
   -- Standard_Integer_Vectors_io.put(ztp); 
   -- put("  #zeroes : "); put(nbz,1); new_line;
    Mixed_Volume(q'last,mix.all,mic,vol);
    if nbz = q'last then
     -- put("Origin is solution with multiplicity "); put(vol,1); new_line;
      ls := new Solution'(Origin(q'last,integer32(vol)));
      Construct(ls,sols);
    else
     -- put("Solution with "); put(nbz,1); put(" zeroes and volume ");
     -- put(vol,1); new_line;
      Silent_Polyhedral_Continuation
        (q,mix,lif,mic,nbz,vol,ztp,sols,verbose=>verbose-1);
    end if;
  exception
    when others => put_line("exception in main silent polyhedral continuation");
                   raise;
  end Silent_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
              ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; -- k : in integer32;
                sols : out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions,QuadDobl_Stable_Homotopies;

    ztp : constant Standard_Integer_Vectors.Vector(q'range)
        := Zero_Type(mic.nor.all,b,mic.pts.all);
    nbz : constant integer32 := Number_of_Zeroes(ztp);
    vol : natural32;
    ls : Link_to_Solution;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 6,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
   -- put("Type of Cell "); put(k,1); put(" : ");
   -- Standard_Integer_Vectors_io.put(ztp); 
   -- put("  #zeroes : "); put(nbz,1); new_line;
    Mixed_Volume(q'last,mix.all,mic,vol);
    if nbz = q'last then
     -- put("Origin is solution with multiplicity "); put(vol,1); new_line;
      ls := new Solution'(Origin(q'last,integer32(vol)));
      Construct(ls,sols);
    else
     -- put("Solution with "); put(nbz,1); put(" zeroes and volume ");
     -- put(vol,1); new_line;
      Silent_Polyhedral_Continuation
        (q,mix,lif,mic,nbz,vol,ztp,sols,verbose=>verbose-1);
    end if;
  exception
    when others => put_line("exception in main silent polyhedral continuation");
                   raise;
  end Silent_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in Standard_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; k : in integer32;
                sols : in out Standard_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions,Standard_Stable_Homotopies;

    ztp : constant Standard_Integer_Vectors.Vector(q'range)
        := Zero_Type(mic.nor.all,b,mic.pts.all);
    nbz : constant integer32 := Number_of_Zeroes(ztp);
    vol : natural32;
    ls : Link_to_Solution;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 4,");
      put_line("for Laurent systems in double precision ...");
    end if;
    put(file,"Type of Cell "); put(file,k,1); put(file," : ");
    Standard_Integer_Vectors_io.put(file,ztp); 
    put(file,"  #zeroes : "); put(file,nbz,1); new_line(file);
    Mixed_Volume(q'last,mix.all,mic,vol);
    if nbz = q'last then
      put(file,"Origin is solution with multiplicity ");
      put(file,vol,1); new_line(file);
      ls := new Solution'(Origin(q'last,integer32(vol)));
      Construct(ls,sols);
    else
      put(file,"Solution with "); put(file,nbz,1);
      put(file," zeroes and volume "); 
      put(file,vol,1); new_line(file);
      Reporting_Polyhedral_Continuation
        (file,q,mix,lif,mic,nbz,vol,ztp,sols,verbose=>verbose-1);
    end if;
  end Reporting_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; k : in integer32;
                sols : in out DoblDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions,DoblDobl_Stable_Homotopies;

    ztp : constant Standard_Integer_Vectors.Vector(q'range)
        := Zero_Type(mic.nor.all,b,mic.pts.all);
    nbz : constant integer32 := Number_of_Zeroes(ztp);
    vol : natural32;
    ls : Link_to_Solution;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 5,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    put(file,"Type of Cell "); put(file,k,1); put(file," : ");
    Standard_Integer_Vectors_io.put(file,ztp); 
    put(file,"  #zeroes : "); put(file,nbz,1); new_line(file);
    Mixed_Volume(q'last,mix.all,mic,vol);
    if nbz = q'last then
      put(file,"Origin is solution with multiplicity ");
      put(file,vol,1); new_line(file);
      ls := new Solution'(Origin(q'last,integer32(vol)));
      Construct(ls,sols);
    else
      put(file,"Solution with "); put(file,nbz,1);
      put(file," zeroes and volume "); 
      put(file,vol,1); new_line(file);
      Reporting_Polyhedral_Continuation
        (file,q,mix,lif,mic,nbz,vol,ztp,sols,verbose=>verbose-1);
    end if;
  end Reporting_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
              ( file : in file_type;
                q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                b : in double_float; 
                mix : in Standard_Integer_Vectors.Link_to_Vector;
                lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                mic : in out Mixed_Cell; k : in integer32;
                sols : in out QuadDobl_Complex_Solutions.Solution_List;
                verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions,QuadDobl_Stable_Homotopies;

    ztp : constant Standard_Integer_Vectors.Vector(q'range)
        := Zero_Type(mic.nor.all,b,mic.pts.all);
    nbz : constant integer32 := Number_of_Zeroes(ztp);
    vol : natural32;
    ls : Link_to_Solution;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 6,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    put(file,"Type of Cell "); put(file,k,1); put(file," : ");
    Standard_Integer_Vectors_io.put(file,ztp); 
    put(file,"  #zeroes : "); put(file,nbz,1); new_line(file);
    Mixed_Volume(q'last,mix.all,mic,vol);
    if nbz = q'last then
      put(file,"Origin is solution with multiplicity ");
      put(file,vol,1); new_line(file);
      ls := new Solution'(Origin(q'last,integer32(vol)));
      Construct(ls,sols);
    else
      put(file,"Solution with "); put(file,nbz,1);
      put(file," zeroes and volume "); 
      put(file,vol,1); new_line(file);
      Reporting_Polyhedral_Continuation
        (file,q,mix,lif,mic,nbz,vol,ztp,sols,verbose=>verbose-1);
    end if;
  end Reporting_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
               ( q : in Standard_Complex_Laur_Systems.Laur_Sys;
                 b : in double_float; 
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Mixed_Subdivision;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions,Standard_Stable_Homotopies;

    tmp : Mixed_Subdivision := mcc;
    last : Solution_List := sols;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 7,");
      put_line("for Laurent systems in double precision ...");
    end if;
    for i in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        csols : Solution_List;
      begin
       -- Silent_Polyhedral_Continuation(q,b,mix,lif,mic,integer32(i),csols);
        Silent_Polyhedral_Continuation
          (q,b,mix,lif,mic,csols,verbose=>verbose-1);
        Merge_and_Concat(sols,last,csols);
        Clear(csols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Silent_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
               ( q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 b : in double_float; 
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Mixed_Subdivision;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions,DoblDobl_Stable_Homotopies;

    tmp : Mixed_Subdivision := mcc;
    last : Solution_List := sols;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 8,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    for i in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        csols : Solution_List;
      begin
       -- Silent_Polyhedral_Continuation(q,b,mix,lif,mic,integer32(i),csols);
        Silent_Polyhedral_Continuation
          (q,b,mix,lif,mic,csols,verbose=>verbose-1);
        Merge_and_Concat(sols,last,csols);
        Clear(csols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Silent_Polyhedral_Continuation;

  procedure Silent_Polyhedral_Continuation
               ( q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 b : in double_float; 
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Mixed_Subdivision;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions,QuadDobl_Stable_Homotopies;

    tmp : Mixed_Subdivision := mcc;
    last : Solution_List := sols;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Silent_Polyhedral_Continuation 9,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    for i in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        csols : Solution_List;
      begin
       -- Silent_Polyhedral_Continuation(q,b,mix,lif,mic,integer32(i),csols);
        Silent_Polyhedral_Continuation
          (q,b,mix,lif,mic,csols,verbose=>verbose-1);
        Merge_and_Concat(sols,last,csols);
        Clear(csols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Silent_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
               ( file : in file_type;
                 q : in Standard_Complex_Laur_Systems.Laur_Sys;
                 b : in double_float; 
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Mixed_Subdivision;
                 sols : in out Standard_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use Standard_Complex_Solutions,Standard_Stable_Homotopies;

    tmp : Mixed_Subdivision := mcc;
    last : Solution_List := sols;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 7,");
      put_line("for Laurent systems in double precision ...");
    end if;
    for i in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        csols : Solution_List;
      begin
        Reporting_Polyhedral_Continuation
          (file,q,b,mix,lif,mic,integer32(i),csols,verbose=>verbose-1);
        if verbose > 0 then
          put(file,"Length of csols : "); put(file,Length_Of(csols),1);
          new_line(file);
          put(file,"Length of sols before concat : ");
          put(file,Length_Of(sols),1); new_line(file);
        end if;
        Merge_and_Concat(sols,last,csols);
        if verbose > 0 then
          put(file,"Length of sols after concat : ");
          put(file,Length_Of(sols),1); new_line(file);
        end if;
        Clear(csols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Reporting_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
               ( file : in file_type;
                 q : in DoblDobl_Complex_Laur_Systems.Laur_Sys;
                 b : in double_float; 
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Mixed_Subdivision;
                 sols : in out DoblDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use DoblDobl_Complex_Solutions,DoblDobl_Stable_Homotopies;

    tmp : Mixed_Subdivision := mcc;
    last : Solution_List := sols;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 8,");
      put_line("for Laurent systems in double double precision ...");
    end if;
    for i in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        csols : Solution_List;
      begin
        Reporting_Polyhedral_Continuation
          (file,q,b,mix,lif,mic,integer32(i),csols,verbose=>verbose-1);
        Merge_and_Concat(sols,last,csols);
        Clear(csols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Reporting_Polyhedral_Continuation;

  procedure Reporting_Polyhedral_Continuation
               ( file : in file_type;
                 q : in QuadDobl_Complex_Laur_Systems.Laur_Sys;
                 b : in double_float; 
                 mix : in Standard_Integer_Vectors.Link_to_Vector;
                 lif : in Arrays_of_Floating_Vector_Lists.Array_of_Lists;
                 mcc : in Mixed_Subdivision;
                 sols : in out QuadDobl_Complex_Solutions.Solution_List;
                 verbose : in integer32 := 0 ) is

    use QuadDobl_Complex_Solutions,QuadDobl_Stable_Homotopies;

    tmp : Mixed_Subdivision := mcc;
    last : Solution_List := sols;

  begin
    if verbose > 0 then
      put("-> in stable_polyhedral_continuation.");
      put_line("Reporting_Polyhedral_Continuation 9,");
      put_line("for Laurent systems in quad double precision ...");
    end if;
    for i in 1..Length_Of(mcc) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        csols : Solution_List;
      begin
        Reporting_Polyhedral_Continuation
          (file,q,b,mix,lif,mic,integer32(i),csols,verbose=>verbose-1);
        Merge_and_Concat(sols,last,csols);
        Clear(csols);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Reporting_Polyhedral_Continuation;

end Stable_Polyhedral_Continuation;
