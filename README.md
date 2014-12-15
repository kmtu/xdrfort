XDR Fortran Interface with Wrappers
==========

This Fortran interface for reading xdr files (both xtc and trr formats)
is based on the work of [James 'Wes' Barnett](https://github.com/wesbarnett/)
and depends on the [XTC Library](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library).

Barnett's original interface was made only for xtc files. His code and the explanatory article are available on
[James 'Wes' Barnett's blog](http://statthermo.blogspot.jp/2014/03/read-in-gromacs-xtc-files-with-fortran_27.html). 

Inspired by Barnett's well-designed interface, I modified it to incorporate the trr format.
So far only the reading function is available, but adding the writing funciton should be quite straightforward
since it is already implemented in the XTC library.

