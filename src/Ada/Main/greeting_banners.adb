package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.4.01 released 2015-10-07";

  begin
    return res;
  end Version;

end Greeting_Banners;
