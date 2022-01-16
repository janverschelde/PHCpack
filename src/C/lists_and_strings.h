/* prototypes of functions to write arrays of numbers into strings and
   to convert string representations of lists into arrays of numbers */

int intlist2str ( int n, int *d, char *s );
/*
 * DESCRIPTION :
 *   Given in d is an array of n integers,
 *   on return in s is a list representation.
 *
 * REQUIRED :
 *   The string s is in the Python-style format of a list
 *   of integer nubmers, e.g.: [2, 0, 42, -29].
 *
 * ON ENTRY :
 *   n         number of integers in d;
 *   d         d integer numbers.
 *
 * ON RETURN :
 *   s         string with n integers separated by commas
 *             and enclosed by square brackets [ and ].  */

int dbllist2str ( int n, double *d, char *s );
/*
 * DESCRIPTION :
 *   Given in d is an array of n doubles,
 *   on return in s is a list representation.
 *
 * REQUIRED :
 *   The string s is in the Python-style format of a list
 *   of doubles, e.g.: [1.0231e-9, 0.981]. 
 *
 * ON ENTRY :
 *   n         number of doubles in d;
 *   d         d double numbers.
 *
 * ON RETURN :
 *   s         string with n doubles separated by commas
 *             and enclosed by square brackets [ and ].  */

int itemcount ( char *s );
/*
 * DESCRIPTION :
 *   Given in s a string representation of a list,
 *   counts the number of commas and returns the
 *   number of items in the list the string represents. */

void str2intlist ( int n, char *s, int *d );
/*
 * DESCRIPTION :
 *   Scans the string s of n numbers and puts them in d.
 *
 * REQUIRED :
 *   There is sufficient memory allocated in d to hold
 *   n integer numbers.
 *
 * ON ENTRY :
 *   n         number of integers, equals itemcount(s);
 *   s         string representation of a list of integers.
 *
 * ON RETURN :
 *   d         contains the integers read from s. */

void str2dbllist ( int n, char *s, double *d );
/*
 * DESCRIPTION :
 *   Scans the string s of n numbers and puts them in d.
 *
 * REQUIRED :
 *   There is sufficient memory allocated in d to hold
 *   n double floating-point numbers.
 *
 * ON ENTRY :
 *   n         number of doubles, equals itemcount(s);
 *   s         string representation of a list of doubles.
 *
 * ON RETURN :
 *   d         contains the integers read from s. */
