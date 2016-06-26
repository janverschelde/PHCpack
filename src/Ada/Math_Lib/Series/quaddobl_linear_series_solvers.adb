with Standard_Floating_Numbers;         use Standard_Floating_Numbers;

package body QuadDobl_Linear_Series_Solvers is

  function cabs ( c : Complex_Number ) return quad_double is
  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

  function cabs ( s : Series ) return quad_double is
  begin
    return cabs(s.cff(0));
  end cabs;

  procedure LUfac ( A : in out QuadDobl_Dense_Series_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

    kp1,L,nm1 : integer32;
    smax,ikabs : quad_double;
    minone : constant quad_double := create(-1.0);
    nbr : Complex_Number;
    temp : Series;

  begin
    info := 0;
    nm1 := n - 1;
    if nm1 >= 1 then
      for k in 1..nm1 loop
        kp1 := k + 1;                              -- find the pivot index L
        L := k;
        smax := cabs(A(k,k));
        for i in kp1..n loop
          ikabs := cabs(A(i,k));
          if ikabs > smax then
            L := i;
            smax := ikabs;
          end if;
        end loop;
        ipvt(k) := L;
        if is_zero(smax) then       -- this column is already triangularized
          info := k;
        else
          if L /= k then                         -- interchange if necessary
            temp := A(L,k);
            A(L,k) := A(k,k);
            A(k,k) := temp;
          end if;
          nbr := Create(minone);
          temp := nbr/A(k,k);                         -- compute multipliers
          for i in kp1..n loop
            A(i,k) := temp*A(i,k);
          end loop;
          for j in kp1..n loop       -- row elimination with column indexing
            temp := A(L,j);
            if L /= k then
              A(L,j) := A(k,j);
              A(k,j) := temp;
            end if;
            for i in kp1..n loop
              A(i,j) := A(i,j) + temp*A(i,k);
            end loop;
          end loop;
        end if;
      end loop;
    end if;
    ipvt(n) := n;
    if is_zero(cabs(A(n,n)))
     then info := n;
    end if;
  end LUfac;

  procedure LUsolve ( A : in QuadDobl_Dense_Series_Matrices.Matrix;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out QuadDobl_Dense_Series_Vectors.Vector ) is

    ell,nm1,kb : integer32;
    temp : Series;
 
  begin
    nm1 := n-1;
    if nm1 >= 1 then                                       -- solve L*y = b
      for k in 1..nm1 loop
        ell := ipvt(k);
        temp := b(ell);
        if ell /= k then
          b(ell) := b(k);
          b(k) := temp;
        end if;
        for i in (k+1)..n loop
          b(i) := b(i) + temp*A(i,k);
        end loop;
      end loop;
    end if;
    for k in 1..n loop                                     -- solve U*x = y
      kb := n+1-k;
      b(kb) := b(kb)/A(kb,kb);
      temp := -b(kb);
      for j in 1..(kb-1) loop
        b(j) := b(j) + temp*A(j,kb);
      end loop;
    end loop;
  end LUsolve;

end QuadDobl_Linear_Series_Solvers;
