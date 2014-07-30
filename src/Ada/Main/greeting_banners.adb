package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.90 released 2014-07-30";

  begin
    return res;
  end Version;

end Greeting_Banners;
