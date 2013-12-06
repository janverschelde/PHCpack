package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.85 released 2013-12-06";

  begin
    return res;
  end Version;

end Greeting_Banners;
