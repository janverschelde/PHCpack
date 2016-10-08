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

end Binomial_Coefficients;
