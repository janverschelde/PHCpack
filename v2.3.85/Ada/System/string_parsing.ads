package String_Parsing is

-- DESCRIPTION :
--   This package provides functions for parsing a string,
--   in particular looking for a banner.

  function Scan ( s : string; banner : string ) return integer;

  -- DESCRIPTION :
  --   Returns -1 if the banner does not occur in s.
  --   Otherwise returns the index in the string s.  This index marks
  --   the end of the first occurrence of banner in the string s.
  --   This value of s(index) equals the last character of banner.

end String_Parsing;
