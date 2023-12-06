package body Binomial_Coefficients is

  function binomial ( n,k : integer32 ) return integer32 is

    quot,prod : integer32 := 1;

  begin
    for i in 1..n-k loop
      quot := i*quot;
    end loop;
    for i in k+1..n loop
      prod := i*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return double_float is

    quot,prod : double_float := 1.0;

  begin
    for i in 1..n-k loop
      quot := double_float(i)*quot;
    end loop;
    for i in k+1..n loop
      prod := double_float(i)*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return double_double is

    quot,prod : double_double := create(integer(1));
    i_dd : double_double;

  begin
    for i in 1..n-k loop
      i_dd := create(integer(i));
      quot := i_dd*quot;
    end loop;
    for i in k+1..n loop
      i_dd := create(integer(i));
      prod := i_dd*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return triple_double is

    quot,prod : triple_double := create(integer(1));
    i_dd : triple_double;

  begin
    for i in 1..n-k loop
      i_dd := create(integer(i));
      quot := i_dd*quot;
    end loop;
    for i in k+1..n loop
      i_dd := create(integer(i));
      prod := i_dd*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return quad_double is

    quot,prod : quad_double := create(integer(1));
    i_qd : quad_double;

  begin
    for i in 1..n-k loop
      i_qd := create(integer(i));
      quot := i_qd*quot;
    end loop;
    for i in k+1..n loop
      i_qd := create(integer(i));
      prod := i_qd*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return penta_double is

    quot,prod : penta_double := create(integer(1));
    i_pd : penta_double;

  begin
    for i in 1..n-k loop
      i_pd := create(integer(i));
      quot := i_pd*quot;
    end loop;
    for i in k+1..n loop
      i_pd := create(integer(i));
      prod := i_pd*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return octo_double is

    quot,prod : octo_double := create(integer(1));
    i_od : octo_double;

  begin
    for i in 1..n-k loop
      i_od := create(integer(i));
      quot := i_od*quot;
    end loop;
    for i in k+1..n loop
      i_od := create(integer(i));
      prod := i_od*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return deca_double is

    quot,prod : deca_double := create(integer(1));
    i_da : deca_double;

  begin
    for i in 1..n-k loop
      i_da := create(integer(i));
      quot := i_da*quot;
    end loop;
    for i in k+1..n loop
      i_da := create(integer(i));
      prod := i_da*prod;
    end loop;
    return prod/quot;
  end binomial;

  function binomial ( n,k : integer32 ) return hexa_double is

    quot,prod : hexa_double := create(integer(1));
    i_da : hexa_double;

  begin
    for i in 1..n-k loop
      i_da := create(integer(i));
      quot := i_da*quot;
    end loop;
    for i in k+1..n loop
      i_da := create(integer(i));
      prod := i_da*prod;
    end loop;
    return prod/quot;
  end binomial;

end Binomial_Coefficients;
