with Standard_Common_Divisors;           use Standard_Common_Divisors;
--with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;

package body Standard_Integer_Linear_Equalities is

  procedure Triangulate ( L : in integer32; m : in matrix;
                          index : in integer32; ineq : in out matrix ) is

    a,b,lcmab,faca,facb : integer32;

  begin
    if ineq(index,L) /= 0 then
      a := m(L,L); b := ineq(index,L);
      lcmab := lcm(a,b);
      if lcmab < 0 then lcmab := -lcmab; end if;
      faca := lcmab/a;  facb := lcmab/b;
      if facb < 0
       then facb := -facb; faca := -faca;
      end if;
      for j in L..ineq'last(2) loop
        ineq(index,j) := facb*ineq(index,j) - faca*m(L,j);
      end loop;
    end if;
  end Triangulate;

  procedure Triangulate ( L : in integer32; m : in matrix;
                          first,last : in integer32; ineq : in out matrix ) is
  begin
    for i in first..last loop
      if ineq(i,l) /= 0
       then Triangulate(l,m,i,ineq);
      end if;
    end loop;
  end Triangulate;

  procedure Triangulate ( index,start : in integer32; ineq : in out matrix ) is
  
    column : integer32 := start;  -- current column counter
    firstineq : integer32 := ineq'first(1);
    found : boolean;
    ind2,a,b,lcmab,faca,facb : integer32;

  begin
    --put_line("INEQUALTITY");
    --put(" BEFORE UPDATE : ");
    --for l in ineq'range(2) loop
    --  put(ineq(index,l),1); put(' ');
    --end loop;
    --new_line;
    loop
     -- SEARCH FOR FIRST NONZERO ENTRY IN CURRENT INEQUALITY :
      while column < ineq'last(2) and then ineq(index,column) = 0 loop
        column := column + 1;
      end loop;
      exit when ((ineq(index,column) = 0) or else (column = ineq'last(2)));
                                        -- nothing to eliminate
     -- SEARCH FOR INEQUALITY,
     --    WITH SAME FIRST NONZERO COLUMN, BUT WITH OPPOSITE SIGN :
      found := false;
      for k in firstineq..(index-1) loop
        if ineq(index,column)*ineq(k,column) < 0 then -- check for sign
          found := true;
          for l in start..column-1 loop      -- check if same zero pattern
            -- if ineq(index,l) = 0
            --  then 
            found := (ineq(k,l) = 0);
            -- end if;
            exit when not found;
          end loop;
          if found then ind2 := k; end if;
        end if;
        exit when found;
      end loop;
      exit when not found;  -- no possibility for elimination
     -- if found
     --  then
        -- ELIMINATE BY MAKING A POSITIVE LINEAR COMBINATION :
         a := ineq(index,column);  b := ineq(ind2,column);
         if a < 0
          then lcmab := lcm(-a,b);
               faca := lcmab/(-a); facb := lcmab/b;
          else lcmab := lcm(a,-b);
               facb := lcmab/(-b); faca := lcmab/a;
         end if;
         if ineq(index,ineq'last(2)) >= 0
           or else  -- PRESERVE SIGN OF ineq(index,ineq'last(2)) !!!!
                (faca*ineq(index,ineq'last(2)) 
                            + facb*ineq(ind2,ineq'last(2)) < 0)
          then
            for l in start..ineq'last(2) loop
              ineq(index,l) := faca*ineq(index,l) + facb*ineq(ind2,l);
            end loop;
           -- PROCEED :
            column := column + 1;
            firstineq := ineq'first(1);
          else
           -- TRY TO ELIMINATE WITH OTHER INEQUALITIES :
            firstineq := ind2 + 1;
         end if;
         if (firstineq >= index)
           -- impossible to eliminate with sign preservation
          then firstineq := ineq'first(2);  
               column := column + 1;
         end if;
     --  else
     --    column := column + 1;
     --    firstineq := ineq'first(2);
     -- end if;
      exit when (column >= ineq'last(2)-1);
    end loop;
    --put(" AFTER UPDATE : ");
    --for l in ineq'range(2) loop
    --  put(ineq(index,l),1); put(' ');
    --end loop;
    --new_line;
  end Triangulate;

end Standard_Integer_Linear_Equalities;
