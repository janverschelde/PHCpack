with Standard_Common_Divisors;           use Standard_Common_Divisors;
with Standard_Integer_Matrices;          use Standard_Integer_Matrices;
with Standard_Integer_Linear_Solvers;    use Standard_Integer_Linear_Solvers;

package body Face_Enumerators_Utilities is

  function Is_Zero ( v : Vector ) return boolean is
  begin
    for i in v'range loop
      if v(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Is_Zero;

  function gcd ( v : Vector ) return integer32 is

    g : integer32 := v(v'first);

  begin
    if g < 0
     then g := -g;
    end if;
    if g /= 1 then
      for i in v'first+1..v'last loop
        if v(i) /= 0 then
          if v(i) < 0
           then g := gcd(g,-v(i));
           else g := gcd(g,v(i));
          end if;
        end if;
        exit when (g = 1);
      end loop;
    end if;
    return g;
  end gcd;

  procedure Scale ( v : in out Vector ) is

    g : constant integer32 := gcd(v);

  begin
    if (g /= 0) and then (g /= 1) then
      for i in v'range loop
        v(i) := v(i)/g;
      end loop;
    end if;
  end Scale;

  function Is_In ( x : integer32; v : Vector ) return boolean is
  begin
    for k in v'range loop
      if v(k) = x
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function In_Edge ( x,a,b : Vector ) return boolean is

    ba,xa,piv : integer32;

  begin
    for i in b'range loop     -- look for first i: b(i) - a(i) /= 0
      ba := b(i) - a(i);
      if ba /= 0
       then piv := i; exit;
      end if;
    end loop;
    if ba = 0 then
      return Equal(x,a);      -- in this case b = a, so test x = a
    else
      if ba < 0 then
        ba := -ba;
        xa := x(piv) - b(piv);
      else
        xa := x(piv) - a(piv);
      end if;
      if xa*ba >= 0 and then xa <= ba
       then return true;
       else return false;
      end if;
    end if;
  end In_Edge;

  function Last_Rows_Zero ( mat : matrix; frst : integer32 ) return boolean is

  -- DESCRIPTION :
  --   Returns true when the last rows of mat, starting at the index frst,
  --   are completely zero.

    function Zero_Row ( mat : matrix; i : integer32 ) return boolean is
    begin
      for j in mat'range(2) loop
        if mat(i,j) /= 0
         then return false;
        end if;
      end loop;
      return true;
    end Zero_Row;

  begin
    for i in frst..mat'last(1) loop
      if not Zero_Row(mat,i)
       then return false;
      end if;
    end loop;
    return true;
  end Last_Rows_Zero;

  function In_Face ( k : in integer32; face,x : Vector; pts : VecVec )
                   return boolean is

    mat : Matrix(x'range,1..k+1);
    sol : Vector(1..k+1) := (1..k+1 => 0);
    sum : integer32 := 0;

  begin
    if k = 0 then
      return (pts(face(face'first)).all = x);
    elsif k = 1 then
      return In_Edge(x,pts(face(face'first)).all,
                       pts(face(face'first+1)).all);
    end if;
   -- put_line("Calling In_Face for ");
   -- put("x : "); put(x); new_line;
   -- for i in face'range loop
   --   put("pts("); put(face(i),1); put(") : "); put(pts(face(i))); new_line;
   -- end loop;
   -- compute the decomposition of x w.r.t. the face :
    for j in 1..k loop
      for i in mat'range(1) loop
        mat(i,j) := pts(face(j))(i) - pts(face(0))(i);
      end loop;
    end loop;
    for i in mat'range(1) loop
      mat(i,k+1) := x(i) - pts(face(0))(i);
    end loop;
   -- put_line("The matrix which defines the decomposition :");
    Upper_Triangulate(mat);
   -- put(mat);
    if not Last_Rows_Zero(mat,k+2) then
      return false;
    else
      declare
        matk1 : matrix(1..k+1,1..k+1);
      begin
        for i in matk1'range(1) loop
          for j in matk1'range(2) loop
            matk1(i,j) := mat(i,j);
          end loop;
        end loop;
        Scale(matk1); Solve0(matk1,sol);
        -- put("The computed solution : "); put(sol); new_line;
        for i in 1..k loop
          if sol(i)*sol(k+1) > 0
           then --put_line("x lies not in the face");
                return false;
           else sum := sum + sol(i);
          end if;
        end loop;
        -- put("The sum : "); put(sum,1); new_line;
        if sum = 0
         then --put_line("x lies not in the face");
              return false;
        end if;
        if sol(k+1) < 0
         then sol(k+1) := -sol(k+1);
         else sum := -sum;
        end if;
        if sum > sol(k+1)
         then --put_line("x lies not in the face");
              return false;
         else --put_line("x lies in the face");
              return true;
        end if;
      end;
    end if;
  end In_Face;

end Face_Enumerators_Utilities;
