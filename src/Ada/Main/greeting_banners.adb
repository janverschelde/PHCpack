package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.08 released 2015-12-29";

  begin
    return res;
  end Version;

end Greeting_Banners;
