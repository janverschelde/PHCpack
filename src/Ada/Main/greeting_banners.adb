package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.00 released 2015-08-31";

  begin
    return res;
  end Version;

end Greeting_Banners;
