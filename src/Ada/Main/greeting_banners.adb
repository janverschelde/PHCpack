package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.96 released 2014-11-26";

  begin
    return res;
  end Version;

end Greeting_Banners;
