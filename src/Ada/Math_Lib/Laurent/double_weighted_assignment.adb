with Ada.Text_IO;                       use Ada.Text_IO;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Boolean_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;

package body Double_Weighted_Assignment is

  function Is_In ( v : Standard_Integer_Vectors.Vector;
                   x : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true if there is some k such that v(k) = x,
  --   returns false otherwise.

  begin
    for k in v'range loop
      if v(k) = x
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  procedure enumerate ( A : in Matrix; k : in integer32;
                        selcols : in out Standard_Integer_Vectors.Vector;
                        selsum : in double_float;
                        mincols : out Standard_Integer_Vectors.Vector;
                        minsum : in out double_float;
                        vrblvl : in integer32 := 0 ) is
  begin
    if vrblvl > 0 then
      put("-> in Double_Weighted_Assignment.enumerate, k = ");
      put(k,1); put_line(" ...");
    end if;
    if k > A'last(1) then
      if selsum < minsum then
        minsum := selsum;
        mincols := selcols;
        if vrblvl > 0 then
          put("Minimum : "); put(minsum,1,3,3);
          put(" at"); put(selcols); new_line;
        end if;
      elsif selsum = minsum then
        if vrblvl > 0 then
          put("Minimum : "); put(minsum,1,3,3);
          put(" at"); put(selcols); new_line;
        end if;
      end if;
    else
      for i in A'range(2) loop
        if not Is_In(selcols,i) then
          selcols(k) := i;
          enumerate(A,k+1,selcols,selsum+A(k,i),mincols,minsum,vrblvl);
          selcols(k) := 0;
        end if;
      end loop;
    end if;
  end enumerate;

  function Value_Selection
             ( A : Matrix; cols : Standard_Integer_Vectors.Vector ) 
             return double_float is

    res : double_float := 0.0;

  begin
    for i in cols'range loop
      res := res + A(i,cols(i));
    end loop;
    return res;
  end Value_Selection;

  procedure Hungarian ( A : in Matrix;
                        u,v : out Standard_Floating_Vectors.Vector;
                        p,way : out Standard_Integer_Vectors.Vector;
                        vrblvl : in integer32 := 0 ) is

    maxfloat : constant double_float := 1.0E+99;
    minv : Standard_Floating_Vectors.Vector(0..A'last(2));
    used : Boolean_Vectors.Vector(0..A'last(2));
    i0,j0,j1,innercnt : integer32;
    vdelta,current : double_float;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Weighted_Assignment.hungarian ...");
    end if;
    u := (u'range => 0.0);
    v := (v'range => 0.0);
    p := (p'range => 0);
    way := (way'range => 0);
    for i in A'range(1) loop
      if vrblvl > 0
       then put("step "); put(i,1); put_line(" ...");
      end if;
      p(0) := i; j0 := 0;
      minv := (minv'range => maxfloat);
      used := (used'range => false);
      innercnt := 0;
      loop
        innercnt := innercnt + 1;
        used(j0) := true;
        i0 := p(j0); j1 := 0;
        vdelta := maxfloat;
        for j in A'range(2) loop
          if not used(j) then
            current := A(i0,j) - u(i0) - v(j);
            if current < minv(j)
             then minv(j) := current; way(j) := j0;
            end if;
            if minv(j) < vdelta
             then vdelta := minv(j); j1 := j;
            end if;
          end if;
        end loop;
        if vrblvl > 100 then -- extra debugging information
          put("used : ");
          for k in used'range loop
            if used(k)
             then put(" T");
             else put(" F");
            end if;
          end loop;
          new_line;
          put("minv : "); put(minv); new_line;
          put("way :"); put(way);
          put(", p :"); put(p); new_line;
        end if;
        for j in 0..A'last(2) loop
          if used(j) then
            u(p(j)) := u(p(j)) + vdelta;
            v(j) := v(j) - vdelta;
          else
            minv(j) := minv(j) - vdelta;
          end if;
        end loop;
        j0 := j1;
        exit when (p(j0) = 0);
      end loop;
      if vrblvl > 0 then
        put("u :"); put(u,3); new_line;
        put("v :"); put(v,3); new_line;
        put("count : "); put(innercnt,1); new_line;
      end if;
      loop
        j1 := way(j0);
        p(j0) := p(j1);
        j0 := j1;
        exit when (j0 <= 0);
      end loop;
      if vrblvl > 0 then
        put("way :"); put(way);
        put(", p :"); put(p); new_line;
      end if;
    end loop;
  end Hungarian;

  function Row_Matching
             ( cols : Standard_Integer_Vectors.Vector;
               nbrows : integer32 )
             return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..nbrows) := (1..nbrows => 0);

  begin
    for i in cols'range loop
      res(cols(i)) := i;
    end loop;
    return res;
  end Row_Matching;

  procedure Cramer_Vector
              ( A : in Matrix;
                b : in Standard_Floating_Vectors.Vector;
                c : out Standard_Floating_Vectors.Vector;
                m : in Standard_Integer_VecVecs.VecVec;
                vrblvl : in integer32 := 0 ) is

    u : Standard_Floating_Vectors.Vector(0..A'last(1));
    v : Standard_Floating_Vectors.Vector(0..A'last(2));
    p : Standard_Integer_Vectors.Vector(0..A'last(2));
    w : Standard_Integer_Vectors.Vector(1..A'last(2));
    match : Standard_Integer_Vectors.Vector(A'range(1));
    mincost : double_float := 0.0;
    wrk : Matrix(A'range(1),A'range(2));

  begin
    if vrblvl > 0
     then put_line("-> in Double_Weighted_Assignment.cramer_vector ...");
    end if;
    Hungarian(A,u,v,p,w,vrblvl-1);
    for i in 1..u'last loop
      mincost := mincost + u(i);
    end loop;
    for i in 1..v'last loop
      mincost := mincost + v(i);
    end loop;
    match := Double_Weighted_Assignment.Row_Matching(p,A'last(1));
    if vrblvl > 0 then
      put("0 : "); 
      put("minimum : "); put(mincost,1,3,3);
      put(" at"); put(match); new_line;
    end if;
    c(0) := mincost;
    for i in match'range loop
      m(0)(i) := match(i);
    end loop;
    for k in A'range(2) loop
      wrk := A;
      for i in A'range(1) loop
        wrk(i,k) := b(i);
      end loop;
      Hungarian(wrk,u,v,p,w,vrblvl-1);
      mincost := 0.0;
      for i in 1..u'last loop
        mincost := mincost + u(i);
      end loop;
      for i in 1..v'last loop
        mincost := mincost + v(i);
      end loop;
      match := Double_Weighted_Assignment.Row_Matching(p,A'last(1));
      if vrblvl > 0 then
        put(k,1); put(" : ");
        put("minimum : "); put(mincost,1,3,3);
        put(" at"); put(match); new_line;
      end if;
      c(k) := mincost;
      for i in match'range loop
        m(k)(i) := match(i);
      end loop;
    end loop;
  end Cramer_Vector;

  function Second_Index
              ( m : Standard_Integer_VecVecs.VecVec )
              return Standard_Integer_Vectors.Vector is

    res : Standard_Integer_Vectors.Vector(1..m'last);
    first : boolean;

  begin
    for k in res'range loop
      res(k) := m(0)(k);
      first := true;
      for i in 1..m'last loop
        if m(i)(k) /= res(k) then   -- different index
          if first then
            res(k) := m(i)(k);
            first := false;
          else
            if m(i)(k) /= m(0)(k)     -- different from first and second
             then res(k) := m'last+1; -- must use b-index
            end if;
          end if;
        end if;
      end loop;
      if res(k) = m(0)(k)
       then res(k) := m'last+1;
      end if;
    end loop;
    return res;
  end Second_Index;

  function Abc_Matrix ( A : Matrix;
                        b,c : Standard_Floating_Vectors.Vector )
                      return Matrix is

    res : Matrix(A'range(1),A'first(2)..A'last(2)+1);

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        res(i,j) := A(i,j) + c(j);
      end loop;
      res(i,res'last(2)) := b(i) + c(0);
    end loop;
    return res;
  end Abc_Matrix;

  procedure Abc_Argmin
              ( abc : in Matrix; tol : in double_float;
                minvals : out Standard_Floating_Vectors.Vector;
                idxmins : out Standard_Integer_Vectors.Vector;
                fail : out boolean; vrblvl : in integer32 := 0 ) is

    minval : double_float;
    cnt : integer32;
    idx : integer32 := 0;

  begin
    if vrblvl > 0
     then put_line("-> in Double_Weighted_Assignment.abc_argmin ...");
    end if;
    fail := false;
    for i in abc'range(1) loop
      minval := abc(i,abc'first(2));
      for j in abc'first(2)..abc'last(2) loop
        if abc(i,j) < minval
         then minval := abc(i,j);
        end if;
      end loop;
      if vrblvl > 0
       then put("min : "); put(minval,2,1,0); put(" at");
      end if;
      cnt := 0;
      for j in abc'range(2) loop
        if abs(abc(i,j) - minval) < tol then
          if vrblvl > 0
           then put(" "); put(j,1);
          end if;
          cnt := cnt + 1;
          if cnt <= 2 then
            idx := idx + 1;
            idxmins(idx) := j;
          end if;
        end if;
      end loop;
      if vrblvl > 0 then
        if cnt >= 2
         then put_line(" okay");
         else put_line(" minimum not attained at least twice, bug!");
        end if;
      end if;
      fail := fail or (cnt /= 2);
      minvals(i) := minval;
    end loop;
  end Abc_Argmin;

  procedure sort ( idx1,idx2 : in out Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Sorts the indices such that idx1(k) < idx2(k), for all k.

    val : integer32;

  begin
    for k in idx1'range loop
      if idx1(k) > idx2(k) then
        val := idx2(k);
        idx2(k) := idx1(k);
        idx1(k) := val;
      end if;
    end loop;
  end sort;

  procedure Abc_Indices
              ( abc : in Matrix; tol : in double_float;
                m : in Standard_Integer_VecVecs.VecVec;
                minvals : in Standard_Floating_Vectors.Vector;
                idx1,idx2 : out Standard_Integer_Vectors.Vector;
                vrblvl : in integer32 := 0 ) is

    corrected : boolean;

  begin
    if vrblvl > 0 then
      put_line("-> in Double_Weighted_Assignment.abc_indices ...");
    end if;
    idx1 := m(0).all;
    idx2 := Double_Weighted_Assignment.second_index(m);
    if vrblvl > 0 then
      put("1st index :"); put(idx1); new_line;
      put("2nd index :"); put(idx2); new_line;
    end if;
    corrected := false;
    for i in minvals'range loop
      if idx2(i) /= abc'last(2) then
        if abs(minvals(i) - abc(i,abc'last(2))) < tol then
          idx2(i) := abc'last(2);
          corrected := true;
        end if;
      end if;
    end loop;
    if vrblvl > 0 then
      if corrected then
        put("corrected :"); put(idx2); new_line;
      end if;
    end if;
    sort(idx1,idx2); -- only sort at end, as idx1 is always correct
  end Abc_Indices;

end Double_Weighted_Assignment;
