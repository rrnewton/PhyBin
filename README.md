% PhyBin (0.3): Binning/Clustering Newick Trees by Topology


PhyBin is a simple command line tool that classifies a set of
  [Newick tree files](http://en.wikipedia.org/wiki/Newick_format) 
by their topology.  The purpose of it is to take a large set of tree
files and browse through the most common tree topologies.

![(Above figure) Trees corresponding to the three largest bins resulting from a
  phybin run.  The file `binXX_YYY`, where `XX` is the rank of the bin and
  `YYY` is the number of trees having that topology.](trees.jpg)

Change Log
==========

In version 0.2, PhyBin was extended to do clustering as well as binning:

 * Computee full all-to-all Robinson-Foulds distance matrices (quickly)
 * Hierarchical clustering of all trees into a tree-of-trees dendrogram based on 
   Robinson Foulds symmetric (tree edit) distance.

In version 0.3, PhyBin gained a number of new features

 * A `--tolerant` mode for computing RF distance matrices even for trees missing taxa. 
 * A `--prune` option for "zooming in" on a specific set of taxa.
 * The `--minboostrap` option was added.


Invoking PhyBin 
===============

PhyBin is a command-line program that produces output in the form of
text files and pdfs, but to produce pdfs (to visualize trees) the
  [GraphViz program](http://www.graphviz.org/),
including the `dot` command, must be installed on the machine.

The following is a simple invocation of PhyBin:

    phybin --bin *.tree -o output_dir/

The input trees can be specified directly on the command-line, or, if the
name of a directory is provided instead, all contained files are
assumed to be trees in Newick format.

PhyBin, at minimum, produces files of the form
`output_dir/clusterXX_YY.tr`, one for each bin.  If
requested, it will also produce visual representations of each bin in
the form `output_dir/clusterXX_YY.pdf`.

Downloading and Installing PhyBin
=================================

The source code to PhyBin can be downloaded here:

  * [Download Source from Hackage](http://hackage.haskell.org/package/phybin)

PhyBin is written in Haskell and if you have 
  [Haskell Platform](http://hackage.haskell.org/platform/).
installed you can install phybin with this one-liner:

    cabal install phybin

Otherwise PhyBin is also available for download as a statically-linked
executable for Mac-OS, Linux, and Windows:

  * [Download Mac-OS Binary](phybin-0.2.11.mac) 
  * [Download Linux Binary](phybin-0.2.11.x86_64_linux)
  * [Download Windows Binary](phybin-0.2.11_windows.exe)
  


Command-line Options
====================

In addition to input files and directories, `phybin` supports a number
of command-line options.  Run "phybin --help" to see these options.
Here is a snapshot of the current help output (version 0.2.11):

    Usage: phybin [OPTION...] files or directories...

    PhyBin takes Newick tree files as input.  Paths of Newick files can
    be passed directly on the command line.  Or, if directories are provided,
    all files in those directories will be read.  Taxa are named based on the files
    containing them.  If a file contains multiple trees, all are read by phybin, and
    the taxa name then includes a suffix indicating the position in the file:
     e.g. FILENAME_0, FILENAME_1, etc.

    When clustering trees, Phybin computes a complete all-to-all Robinson-Foulds distance matrix.
    If a threshold distance (tree edit distance) is given, then a flat set of clusters
    will be produced in files clusterXX_YY.tr.  Otherwise it produces a full dendogram (UNFINISHED).

    Binning mode provides an especially quick-and-dirty form of clustering.
    When running with the --bin option, only exactly equal trees are put in the same cluster.
    Tree pre-processing still applies, however: for example collapsing short branches.

    USAGE NOTES:
     * Currently phybin ignores input trees with the wrong number of taxa.
     * If given a directory as input phybin will assume all contained files are Newick trees.


    Options include:

      -v       --verbose         print WARNINGS and other information (recommended at first)
      -V       --version         show version number
      -o DIR   --output=DIR      set directory to contain all output files (default "./phybin_out/")
	       --selftest        run internal unit tests

				 ----------------------------- Clustering Options ------------------------------
	       --bin             Use simple binning, the cheapest form of 'clustering'
	       --single          Use single-linkage clustering (nearest neighbor)
	       --complete        Use complete-linkage clustering (furthest neighbor)
	       --UPGMA           Use Unweighted Pair Group Method (average linkage) - DEFAULT mode
	       --editdist=DIST   Combine all clusters separated by DIST or less.  Report a flat list of clusters.
				 Irrespective of whether this is activated, a hierarchical clustering (dendogram.pdf) is produced.
				   Select Robinson-Foulds (symmetric difference) distance algorithm:
	       --simple          use direct all-to-all comparison
	       --hashrf          (default) use a variant of the HashRF algorithm for the distance matrix

				 ----------------------------- Visualization --------------------------------
      -g       --graphbins       use graphviz to produce .dot and .pdf output files
      -d       --drawbins        like -g, but open GUI windows to show each bin's tree
      -w       --view            for convenience, "view mode" simply displays input Newick files without binning
	       --showtrees       Print (textual) tree topology inside the nodes of the dendrogram
	       --highlight=FILE  Highlight nodes in the tree-of-trees (dendrogram) consistent with the.
				 given tree file.  Multiple highlights are permitted and use different colors.
	       --interior        Show the consensus trees for interior nodes in the dendogram, rather than just points.

				 ---------------------------- Tree pre-processing -----------------------------
      -n NUM   --numtaxa=NUM     expect NUM taxa for this dataset
      -b LEN   --branchcut=LEN   collapse branches less than LEN

				 --------------------------- Extracting taxa names ----------------------------

      -p NUM   --nameprefix=NUM  Leaf names in the input Newick trees can be gene names, not taxa.
				 Then it is typical to extract taxa names from genes.  This option extracts
				 a prefix of NUM characters to serve as the taxa name.

      -s STR   --namesep=STR     An alternative to --nameprefix, STR provides a set of delimeter characters,
				 for example '-' or '0123456789'.  The taxa name is then a variable-length
				 prefix of each gene name up to but not including any character in STR.

      -m FILE  --namemap=FILE    Even once prefixes are extracted it may be necessary to use a lookup table
				 to compute taxa names, e.g. if multiple genes/plasmids map onto one taxa.
				 This option specifies a text file with find/replace entries of the form
				 "<string> <taxaname>", which are applied AFTER -s and -p.

				 --------------------------- Utility Modes ----------------------------
	       --rfdist          print a Robinson Foulds distance matrix for the input trees
	       --setdiff         for convenience, print the set difference between cluster*.txt files
	       --print           simply print out a concise form of each input tree
	       --printnorms      simply print out a concise and NORMALIZED form of each input tree
	       --consensus       print a strict consensus tree for the inputs, then exit
	       --matching        print a list of tree names that match any --highlight argument


- - - - - - - - - - - - - - -
Authors: Irene and Ryan Newton

Contact email: `irnewton` and `rrnewton` at `indiana` `edu` (with "at" and "dot" inserted).

[Irene's](http://www.bio.indiana.edu/faculty/directory/profile.php?person=irnewton) and 
[Ryan](http://www.cs.indiana.edu/~rrnewton/homepage.html) homepages.

.


