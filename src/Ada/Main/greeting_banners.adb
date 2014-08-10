package body Greeting_Banners is

  function Version return string is

    res : constant string := "PHCv2.3.91 released 2014-08-10";

  begin
    return res;
  end Version;

end Greeting_Banners;
