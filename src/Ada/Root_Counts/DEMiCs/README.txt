
-------------------------------------------------------------
READ file in DEMiCs (Version 0.95)   14 August, 2007
-------------------------------------------------------------

DEMiCs is a software package written in C++
for computing the mixed volume of the support of a polynomial system  
via dynamic enumeration of all mixed cells.


This package provides 

(1) the value of the mixed volume, and 
(2) information about all mixed cells 

of the support of a polynomial system.



This package has been developed by 
Tomohiko Mizutani, Akiko Takeda and Masakazu Kojima 
according to the papers:

(1) T. Mizutani, A. Takeda and M. Kojima,
   ``Dynamic Enumeration of All Mixed Cells,'' 
   Discrete and Computational Geometry 37(3), pp. 351--367 (2007).

(2) T. Mizutani and A. Takeda
   ``DEMiCs: A software package for computing the mixed volume
   via dynamic enumeration of all mixed cells,''
   To appear in IMA Volume on Software for Algebraic Geometry.

You can find the research reports in the directory ``Doc''.

   

E-mail address:

 mizutan8@is.titech.ac.jp
 takeda@is.titech.ac.jp
 kojima@is.titech.ac.jp
  



In this README file, we assume that your operating system is Linux / Unix. 


++++++++++++++++++++++++++++++++++++++++++++++++++++++
(1) INSTALL
++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the directory ``SRC'', we will find ``makefile''. 
When we carry out the command 

> make all

in the directory,
then, the executable file "demics" will be generated in the directory.


++++++++++++++++++++++++++++++++++++++++++++++++++++++
(2) HOW TO USE
++++++++++++++++++++++++++++++++++++++++++++++++++++++

To run the executable file ``demics'' for the input file ``poly.dat''
which is placed in ``SRC'',
we just execute the following command

> demics poly.dat

in ``SRC''. 
Then, this package presents the total number of all mixed cells, 
the value of the mixed volume and cpu time in a terminal.



Furthermore, this package offers 
three options ``-c'', ``-s'' and ``-cs''.

--------------
Option  ``-c'' 
--------------
When we require detailed information about each mixed cell,
this option is used.
The execution of the following command 

> demics -c poly.dat

displays detailed information about all mixed cells in a terminal.

--------------
Option  ``-s'' 
--------------

The option ``-s'' is applied to change a seed number.
The seed number is used to generate a random number for each lifting value.
When we choose a seed number as ``6'' for the input file ``poly.dat'',
the following command is executed

> demics -s poly.dat 6

Note that if no option is selected, the seed number is set as ``1''
automatically.

--------------
Option  ``-cs'' 
--------------
This option is chosen in order to change a seed number, 
and output information about each mixed cell in a terminal.

When a seed number is set to ``2'' for the input file ``poly.dat''
and each mixed cell needs to be displayed,
we carry out the command

> demics -cs poly.dat 2


++++++++++++++++++++++++++++++++++++++++++++++++++++++
(3) INPUT FILE
++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can find some sample files In the directory ``polySys''.
The generators for the sample files can be found in ``polyGen''.
if you type in ``polyGen'',

> make all

you can get the executable files.

In particular, the subdirectory ``polyList_MixedVol'' contains test
problems  from the software package MixedVol, which was developed by 
T. Gao and T. Y. Li based on 

Algorithm 846: MixedVol: A Software Package for Mixed Volume Computation, 
ACM Trans. Math. Software, 31(4), pp. 555 -- 560, (2005).




The input file requires the following information:

(1) the dimension of the polynomial system,
(2) the number of the distinct supports sets,
(3) the number of elements in each support set,
(4) the multiplicity of each support set, and 
(5) all elements of each support set.

As an example, we consider the support sets of
a semi-mixed system of type (2, 1, 1) as follows:

S1 := A1 = A2 
 = {(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1), (0,0,0,0), (1,1,1,1)}

S2 := A3
 = {(2,0,0,0), (0,2,0,0), (0,0,2,0), (0,0,0,2), (0,0,0,0), (2,2,2,2)}

S3 := A4
 = {(3,0,0,0), (0,3,0,0), (0,0,3,0), (0,0,0,3), (0,0,0,0), (3,3,3,3)}.


Then, the input file for S = (S1, S2, S3) is 
written in the following format:


# The dimension or the number of variables
Dim = 4

# The number of the distinct support sets
Support = 3

# The number of elements in each support set
Elem = 6 6 6 

# The multiplicity of each support set
Type = 2 1 1


# The elements of the 1st support set
1 0 0 0 
0 1 0 0 
0 0 1 0 
0 0 0 1 
0 0 0 0 
1 1 1 1 

# The elements of the 2nd support set
2 0 0 0 
0 2 0 0 
0 0 2 0 
0 0 0 2 
0 0 0 0 
2 2 2 2 

# The elements of the 3rd support set
3 0 0 0 
0 3 0 0 
0 0 3 0 
0 0 0 3 
0 0 0 0 
3 3 3 3.





++++++++++++++++++++++++++++++++++++++++++++++++++++++
(4) BUGS OR COMMENTS
++++++++++++++++++++++++++++++++++++++++++++++++++++++

The authors are grateful for any comments from users.
For the bugs, questions and suggestions,
contact with:  

 Tomohiko Mizutani (mizutan8@is.titech.ac.jp),
 Akiko Takeda (takeda@is.titch.ac.jp) or
 Masakazu Kojima (kojima@is.titech.ac.jp).

