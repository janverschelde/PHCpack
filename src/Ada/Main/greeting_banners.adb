package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.03 released 2015-10-22";

  begin
    return res;
  end Version;

end Greeting_Banners;
