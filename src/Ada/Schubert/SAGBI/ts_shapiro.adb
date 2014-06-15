with text_io;                            use text_io;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;        use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;       use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;            use Standard_Random_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_Matrices;         use Standard_Floating_Matrices;
with Standard_Floating_Matrices_io;      use Standard_Floating_Matrices_io;
with Standard_Integer_Vectors; 
with Standard_Floating_Linear_Solvers;   use Standard_Floating_Linear_Solvers;
with Osculating_Planes;                  use Osculating_Planes;

procedure ts_shapiro is

-- DESCRIPTION :
--   Tests the generation of special matrices to test Shapiro's conjecture.

  function Determinant
              ( mat : Matrix; rows : Standard_Integer_Vectors.Vector )
              return double_float is

  -- DESCRIPTION :
  --   Computes the determinant of the matrix obtained by selecting rows.

    res : double_float := 1.0;
    sqm : Matrix(rows'range,rows'range);
    piv : Standard_Integer_Vectors.Vector(rows'range);
    inf : integer32;

  begin
    for i in rows'range loop
      piv(i) := i;
      for j in rows'range loop
        sqm(i,j) := mat(rows(i),j);
      end loop;
    end loop;
    lufac(sqm,rows'last,piv,inf);
    for i in rows'range loop
      res := res*sqm(i,i);
    end loop;
    for i in piv'range loop
      if piv(i) > i
       then res := -res;
      end if;
    end loop;
    return res;
  end Determinant;

  procedure Maximal_Minors ( n,d : in natural32; mat : in Matrix;
                             min,max : out double_float ) is

  -- DESCRIPTION :
  --   Computes all maximal minors of a (nxd)-matrix mat, d < n.

    rows : Standard_Integer_Vectors.Vector(1..integer32(d));
    first : boolean := true;
    mindet,maxdet : double_float;

    procedure Select_Rows ( k,start : in integer32 ) is

      det : double_float;

    begin
      if k > integer32(d) then
        det := Determinant(mat,rows);
        -- put("Minor "); put(rows); put(" equals "); put(det); new_line;
        det := abs(det);
        if first then
          mindet := det; maxdet := det; first := false;
        else
          if det > maxdet then
            maxdet := det;
          elsif det < mindet then
            mindet := det;
          end if;
        end if;
      else
        for j in start..integer32(n) loop
          rows(k) := j;
          Select_Rows(k+1,j+1);
        end loop;
      end if;
    end Select_Rows;

  begin
    Select_Rows(1,1);
    put("Min : "); put(mindet,3,3,3);
    put("  Max : "); put(maxdet,3,3,3);
    put("  Max/Min : "); put(maxdet/mindet,3,3,3); new_line;
    min := mindet; max := maxdet;
  end Maximal_Minors;

  procedure Test_Sample ( n,d : natural32; s : double_float ) is

    cheb_mat : Matrix(1..integer32(n),1..integer32(d))
             := Chebychev_Basis(n,d,s);
    orto_mat : Matrix(1..integer32(n),1..integer32(d))
             := Orthogonal_Basis(n,d,s);
    min,max : double_float;

  begin
    put("Osculating "); put(d,1); put("-plane in ");
    put(n,1); put_line("-space : "); put(cheb_mat);
    put_line("All maximal minors : ");
    Maximal_Minors(n,d,cheb_mat,min,max);
    put("Orthogonal respresentation of osculating ");
    put(d,1); put("-plane in "); put(n,1); put_line("-space : ");
    put(orto_mat);
    put_line("All maximal minors : ");
    Maximal_Minors(n,d,orto_mat,min,max);
  end Test_Sample;

  procedure Sampled_Generation is

    n,d : natural32 := 0;
    s : double_float;
    ans : character;

  begin
    new_line;
    put("Give n : "); get(n);
    put("Give d : "); get(d);
    loop
      s := Random;
      put("The s-value : "); put(s); new_line;
      Test_Sample(n,d,s);
      put("Do you want to test another sample ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Sampled_Generation;

  procedure Best_Sampled_Generation is

    n,d,m : natural32 := 0;
    ans : character;
    s,bestratio,bests : double_float;
    first : boolean;

  begin
    new_line;
    put("Give n : "); get(n);
    put("Give d : "); get(d);
    loop
      put("Give the number of samples : "); get(m);
      first := true;
      for i in 1..m loop
        s := Random;
        put("The s-value : "); put(s); new_line;
        declare
          mat : Matrix(1..integer32(n),1..integer32(d))
              := Chebychev_Basis(n,d,s); 
          min,max,ratio : double_float;
        begin
          Maximal_Minors(n,d,mat,min,max);
          ratio := max/min; 
          if first then
            bestratio := ratio; bests := s; first := false;
          else
            if ratio < bestratio
             then bestratio := ratio; bests := s;
            end if;
          end if;
        end;
      end loop;
      put("Best ratio : "); put(bestratio); 
      put("  for s : "); put(bests); new_line;
      put("Do you want more ratio's to test ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Best_Sampled_Generation;

  procedure Update ( s : in out double_float; inc : in double_float ) is
  begin
    s := s + inc;
    if s >= 1.0
     then s := s - 2.0;
    end if;
  end Update;

  procedure Equidistant_Sampling
              ( n,d,nb : in natural32; inits : in double_float;
                rat : out double_float ) is

  -- DESCRIPTION :
  --   Generates nb equidistant s-values and computes the maximum
  --   of all ratios max/min minors.

    inc : constant double_float := 2.0/double_float(nb);
    mat : Matrix(1..integer32(n),1..integer32(d));
    s : double_float := inits;
    first : boolean := true;
    min,max,ratio,maxratio : double_float;

  begin
    for i in 1..nb loop
      mat := Chebychev_Basis(n,d,s);
      Maximal_Minors(n,d,mat,min,max);
      ratio := max/min;
      if first then
        maxratio := ratio; first := false;
      else
        if ratio > maxratio
         then maxratio := ratio;
        end if;
      end if;
      Update(s,inc);
    end loop;
    rat := maxratio;
  end Equidistant_Sampling;

  procedure Init ( mat : in out Matrix ) is
  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        mat(i,j) := 0.0;
      end loop;
    end loop;
  end Init;

  procedure Div ( mat : in out Matrix; d : in double_float ) is
  begin
    for i in mat'range(1) loop
      for j in mat'range(2) loop
        mat(i,j) := mat(i,j)/d;
      end loop;
    end loop;
  end Div;

  procedure Averaged_Equidistant_Sampling
              ( n,d,nb,nav : in natural32; inits : in double_float;
                rat : out double_float ) is

  -- DESCRIPTION :
  --   Generates nb equidistant s-values and computes the maximum
  --   of all ratios max/min minors.
  --   Averages over every interval.

    inc : constant double_float := 2.0/double_float(nb*nav);
    mat : Matrix(1..integer32(n),1..integer32(d));
    s : double_float := inits;
    first : boolean := true;
    min,max,ratio,maxratio : double_float;

  begin
    for i in 1..nb loop
      Init(mat);
      for j in 1..nav loop
        mat := mat + Chebychev_Basis(n,d,s);
        Update(s,inc);
      end loop;
      Div(mat,double_float(nav));
      Maximal_Minors(n,d,mat,min,max);
      ratio := max/min;
      if first then
        maxratio := ratio; first := false;
      else
        if ratio > maxratio
         then maxratio := ratio;
        end if;
      end if;
      Update(s,inc);
    end loop;
    rat := maxratio;
  end Averaged_Equidistant_Sampling;

  procedure Best_Equidistant_Sampling is

    n,d,m,nb : natural32 := 0;
    ans : character;
    s,ratio,bestratio,bests : double_float;
    first : boolean;

  begin
    new_line;
    put("Give n : "); get(n);
    put("Give d : "); get(d);
    loop
      put("Give #samples : "); get(m);
      put("Give #equidistant points : "); get(nb);
      first := true;
      for i in 1..m loop
        s := Random;
        Equidistant_Sampling(n,d,nb,s,ratio);
        if first then
          bestratio := ratio; bests := s; first := false;
        else
          if ratio < bestratio
           then bestratio := ratio; bests := s;
          end if;
        end if;
      end loop;
      put("Best ratio : "); put(bestratio); 
      put("  for s : "); put(bests); new_line;
      put("Do you want more ratio's to test ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Best_Equidistant_Sampling;

  procedure Best_Averaged_Equidistant_Sampling is

    n,d,m,nb,nav : natural32 := 0;
    ans : character;
    s,ratio,bestratio,bests : double_float;
    first : boolean;

  begin
    new_line;
    put("Give n : "); get(n);
    put("Give d : "); get(d);
    put("Give avg #times : "); get(nav);
    loop
      put("Give #samples : "); get(m);
      put("Give #equidistant points : "); get(nb);
      first := true;
      for i in 1..m loop
        s := Random;
        Averaged_Equidistant_Sampling(n,d,nb,nav,s,ratio);
        if first  then
          bestratio := ratio; bests := s; first := false;
        else
          if ratio < bestratio
           then bestratio := ratio; bests := s;
          end if;
        end if;
      end loop;
      put("Best ratio : "); put(bestratio); 
      put("  for s : "); put(bests); new_line;
      put("Do you want more ratio's to test ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Best_Averaged_Equidistant_Sampling;

  function Distance ( v : Standard_Floating_Vectors.Vector;
                      k : integer32 ) return double_float is

    min : double_float := 2.0;

  begin
    for i in v'first..(k-1) loop
      if abs(v(i)-v(k)) < min
       then min := abs(v(i)-v(k));
      end if;
    end loop;
    return min;
  end Distance;

  procedure Spaced_Sampling
              ( n,d,nb : in natural32; mindist : in double_float;
                rat : out double_float ) is

  -- DESCRIPTION :
  --   Generates nb distinct s-values, not closer to each other than
  --   mindist and computes the maximum of all ratios max/min minors.

    mat : Matrix(1..integer32(n),1..integer32(d));
    sva : Standard_Floating_Vectors.Vector(1..integer32(nb));
    first : boolean := true;
    min,max,ratio,maxratio : double_float;

  begin
    for i in 1..integer32(nb) loop
      loop
        sva(i) := Random;
        exit when (Distance(sva,i) > mindist);
      end loop;
      mat := Chebychev_Basis(n,d,sva(i));
      Maximal_Minors(n,d,mat,min,max);
      ratio := max/min;
      if first then
        maxratio := ratio; first := false;
      else
        if ratio > maxratio
         then maxratio := ratio;
        end if;
      end if;
	end loop;
    rat := maxratio;
  end Spaced_Sampling;

  procedure Best_Spaced_Sampling is

    n,d,m,nb : natural32 := 0;
    ans : character;
    inc,mindist,ratio,bestratio : double_float := 0.0;
    first : boolean;

  begin
    new_line;
    put("Give n : "); get(n);
    put("Give d : "); get(d);
    loop
      put("Give #samples : "); get(m);
      put("Give #equidistant points : "); get(nb);
      inc := 2.0/double_float(nb);
      put("The increment : "); put(inc); new_line;
      put("Give Minimal distance : "); get(mindist);
      first := true;
      for i in 1..m loop
        Spaced_Sampling(n,d,nb,mindist,ratio);
        if first then
          bestratio := ratio; first := false;
        else
          if ratio < bestratio
           then bestratio := ratio;
          end if;
        end if;
      end loop;
      put("Best ratio : "); put(bestratio);  new_line;
      put("Do you want more ratio's to test ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Best_Spaced_Sampling;

  procedure Interactive_Generation is

    n,d : natural32 := 0;
    s : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put("Give n : "); get(n);
    put("Give d : "); get(d);
    loop
      put("Give s-value : "); get(s);
      Test_Sample(n,d,s);
      put("Do you want to test another s-value ? (y/n) "); get(ans);
      exit when (ans /= 'y');
    end loop;
  end Interactive_Generation;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Generation of osculating d-planes in n-space.");
    loop
      new_line;
      put_line("Choose one of the following : ");
      put_line("  0. Exit this program.");
      put_line("  1. Interactively input of s-values.");
      put_line("  2. Sampling s-values one after the other.");
      put_line("  3. Sample many times for best ratio.");
      put_line("  4. Equidistant Sampling s-values and select min ratio.");
      put_line("  5. Averaged Equidistant Sampling s-values + min ratio.");
      put_line("  6. Spaced Sampling s-values for min ratio.");
      put("Make your choice (0,1,2,3,4,5 or 6) : "); get(ans);
      exit when (ans = '0');
      case ans is
        when '1' => Interactive_Generation;
        when '2' => Sampled_Generation;
        when '3' => Best_Sampled_Generation;
        when '4' => Best_Equidistant_Sampling;
        when '5' => Best_Averaged_Equidistant_Sampling;
        when '6' => Best_Spaced_Sampling;
        when others => null;
      end case;
    end loop;
  end Main;

begin
  Main;
end ts_shapiro;
