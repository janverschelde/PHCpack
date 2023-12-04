with Standard_Floating_Numbers;           use Standard_Floating_Numbers;

package body HexaDobl_Series_Linear_Solvers is

  function cabs ( c : Complex_Number ) return hexa_double is
  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

  function cabs ( s : Series ) return hexa_double is
  begin
    return cabs(s.cff(0));
  end cabs;

  procedure LUfac ( A : in out HexaDobl_Complex_Series_Matrices.Matrix;
                    n : in integer32;
                    ipvt : out Standard_Integer_Vectors.Vector;
                    info : out integer32 ) is

    zero : constant hexa_double := create(0.0);
    minone : constant hexa_double := create(-1.0);
    kp1,L,nm1 : integer32;
    smax,ikabs : hexa_double;
    nbr : Complex_Number;
    temp : Link_to_Series;
    deg : constant integer32 := A(A'first(1),A'first(2)).deg;
    fac,work : Series(deg);

  begin
    info := 0;
    nm1 := n - 1;
    if nm1 >= 1 then
      for k in 1..nm1 loop
        kp1 := k + 1;                              -- find the pivot index L
        L := k;
        smax := cabs(A(k,k).all);
        for i in kp1..n loop
          ikabs := cabs(A(i,k).all);
          if ikabs > smax then
            L := i;
            smax := ikabs;
          end if;
        end loop;
        ipvt(k) := L;
        if smax = zero then         -- this column is already triangularized
          info := k;
        else
          if L /= k then                         -- interchange if necessary
            temp := A(L,k);
            A(L,k) := A(k,k);
            A(k,k) := temp;
          end if;
          nbr := Create(minone);
          fac := nbr/A(k,k).all;                      -- compute multipliers
          for i in kp1..n loop
            Mul(A(i,k).all,fac);       -- A(i,k) := temp*A(i,k);
          end loop;
          for j in kp1..n loop       -- row elimination with column indexing
            temp := A(L,j);
            if L /= k then
              A(L,j) := A(k,j);
              A(k,j) := temp;
            end if;
            for i in kp1..n loop     -- A(i,j) := A(i,j) + temp*A(i,k);
              work := temp.all*A(i,k).all;
              Add(A(i,j).all,work);
            end loop;
          end loop;
        end if;
      end loop;
    end if;
    ipvt(n) := n;
    if cabs(A(n,n).all) = zero 
     then info := n;
    end if;
  end LUfac;

  procedure LUsolve ( A : in HexaDobl_Complex_Series_Matrices.Matrix;
                      n : in integer32;
                      ipvt : in Standard_Integer_Vectors.Vector;
                      b : in out HexaDobl_Complex_Series_Vectors.Vector ) is

    ell,nm1,kb : integer32;
    temp : Link_to_Series;
    deg : constant integer32 := b(b'first).deg;
    work,fac : Series(deg);
 
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
        for i in (k+1)..n loop          -- b(i) := b(i) + temp*A(i,k);
          work := temp.all*A(i,k).all;
          Add(b(i).all,work);
        end loop;
      end loop;
    end if;
    for k in 1..n loop                                     -- solve U*x = y
      kb := n+1-k;
      Div(b(kb),A(kb,kb));              -- b(kb) := b(kb)/A(kb,kb);
      fac := -b(kb).all;
      for j in 1..(kb-1) loop           -- b(j) := b(j) + temp*A(j,kb);
        work := fac*A(j,kb).all;
        Add(b(j).all,work);
      end loop;
    end loop;
  end LUsolve;

end HexaDobl_Series_Linear_Solvers;
