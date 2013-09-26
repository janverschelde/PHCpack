with text_io,integer_io;                use text_io,integer_io;
with Communications_with_User;          use Communications_with_User;

procedure ts_nbterms is

-- DESCRIPTION :
--   This procedure counts number of monomials and number of samples
--   needed to interpolate multivariate polynomials of any degree.

  function Polynomial_Size ( n,d : natural ) return natural is

    sum : natural;

  begin
    if n = 2
     then if d = 0
           then return 1;
           else return ((d+1)*(d+2)/2);
          end if;
     else sum := 0;
          for i in 0..d loop
            sum := sum + Polynomial_Size(n-1,i);
          end loop;
          return sum;
    end if;
  end Polynomial_Size;

  function Grid_Size ( n,d : in natural ) return natural is

  -- DESCRIPTION :
  --   Returns the grid size of samples to interpolate a polynomial in
  --   n variables of degree d for use with the bootstrapping Newton.

  -- REQUIRED : n > 0.

    sum : natural;

  begin
    if n = 2
     then return (d+1)*d;
     else sum := 1;              -- extra point
          for i in 1..d loop
            sum := sum + Grid_Size(n-1,d);
          end loop;
          return sum;
    end if;
  end Grid_Size;

  function Full_Grid_Size ( n,d : in natural ) return natural is

  -- DESCRIPTION :
  --   Returns d*(d+1)^(n-1), which is the size of a full grid
  --   needed for the trace form in n variables of degree d.

  -- REQUIRED : n > 0.

    prod : natural := d;

  begin
    for i in 1..n-1 loop
      prod := prod*(d+1);
    end loop;
    return prod;
  end Full_Grid_Size;

  function Linear_Grid_Size ( n,d : in natural ) return natural is

  -- DESCRIPTION :
  --   Returns n*d, which is the size of the grid needed
  --   for the linear trace in n variables for degree d surface.

  -- REQUIRED : n > 0.

  begin
    return n*d;
  end Linear_Grid_Size;

  procedure Display_Table ( kind : in natural ) is

  -- DESCRIPTION :
  --   Determines the size of an interpolation grid of samples needed
  --   to interpolation multivariate polynomial of any degree.
  --   The parameter kind has the following meaning :
  --     0 : number of monomials
  --     1 : number of samples with bootstrapping Newton
  --     2 : number of samples with traces on full grid
  --     3 : number of samples for linear trace

    maxdeg,maxdim : natural;
    colwidth : constant natural := 10;
    nb : natural;

  begin
    new_line;
    put_line("Number of samples for increasing degrees and dimensions.");
    new_line;
    put("Give the maximal degree : "); get(maxdeg);
    put("Give the maximal dimension : "); get(maxdim);
    new_line;
    put(" d\n  ");
    for j in 2..maxdim loop
      put(j,colwidth);
    end loop;
    new_line;
    for i in 1..maxdeg loop
      put(i,3); put(" : ");
      for j in 2..maxdim loop
        case kind is
          when 0 => nb := Polynomial_Size(j,i);
          when 1 => nb := Grid_Size(j,i);
          when 2 => nb := Full_Grid_Size(j,i);
          when 3 => nb := Linear_Grid_Size(j,i);
          when others => nb := 0;
        end case;
        put(nb,colwidth);
      end loop;
      new_line;
    end loop;
  end Display_Table;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Counting monomials and samples needed for interpolation.");
    new_line;
    put_line("MENU for testing enumerators :");
    put_line("  0. Generate table with #terms for degrees and dimensions.");
    put_line("  1. Number of samples with bootstrapping Newton.");
    put_line("  2. Number of samples with trace form on full grid.");
    put_line("  3. Number of samples for the linear trace.");
    put("Type 0, 1, 2, or 3 to select : "); Ask_Alternative(ans,"0123");
    case ans is
      when '0' => Display_Table(0);
      when '1' => Display_Table(1);
      when '2' => Display_Table(2);
      when '3' => Display_Table(3);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_nbterms;
