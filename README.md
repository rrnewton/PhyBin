% PhyBin: Binning Trees by Topology


PhyBin is a simple command line tool that classifies (bins) a set of
  [Newick tree files](http://en.wikipedia.org/wiki/Newick_format) 
by their topology.  The purpose of it is to take a large set of tree
files and browse through the most common tree topologies.

![(Above figure) Trees corresponding to the three largest bins resulting from a
  phybin run.  The file `binXX_YYY`, where `XX` is the rank of the bin and
  `YYY` is the number of trees having that topology.](trees.jpg)

Invoking PhyBin 
===============

PhyBin is a command-line program that produces output in the form of
text files and pdfs, but to produce pdfs (to visualize trees) the
  [GraphViz program](http://www.graphviz.org/),
including the `dot` command, must be installed on the machine.

The following is a typical invocation of PhyBin:

    phybin *.tree -o output_dir/

The input trees can be specified directly on the command-line, or, if the
name of a directory is provided instead, all contained files are
assumed to be trees in Newick format.

PhyBin, at minimum, produces files of the form
`output_dir/binXX_YY.tr`, one for each bin.  If
requested, it will also produce visual representations of each bin in
the form `output_dir/binXX_YY.pdf`.

Downloading and Installing PhyBin
=================================

The source code to PhyBin can be downloaded here:

  * [Download Source Tarball](phybin-0.1.2.tar.gz)

PhyBin is written in Haskell and if you have 
  [Haskell Platform](http://hackage.haskell.org/platform/).
installed you can install phybin with this one-liner:

    cabal install phybin

PhyBin is also available for download as a statically-linked
executable for Mac-OS and Linux:

  * [Download Mac-OS Binary](phybin-0.1.2.mac) 
  * [Download Linux Binary](phybin-0.1.2.x86_64_linux)
  
It should be possible to build it for Windows, but I haven't done so
yet.



Command-line Options
====================

In addition to input files and directories, `phybin` supports a number
of command-line options.  Run "phybin --help" to see these options.

- - - - - - - - - - - - - - -
Authors: Irene and Ryan Newton

Contact email: `irnewton` `indiana` `edu` (with "at" and "dot" inserted).

.


