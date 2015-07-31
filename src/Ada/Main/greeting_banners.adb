package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.99 released 2015-07-31";

  begin
    return res;
  end Version;

end Greeting_Banners;
