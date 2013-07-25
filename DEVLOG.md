

[2012.11.19]
------------

Right now the rickettsia dataset is causing problems:

    phybin.exe: Rickettsia/phybin_output/bin1_16.txt: openFile: does not exist (No such file or directory)

[2013.07.22] {Beginning optimization of all-to-all RF distance}
---------------------------------------------------------------

Making the labels into Ints and switching to IntSet has sped up the
150-taxa/100-trees example from 3 sec to 1.6 seconds on my macbook
retina.


[2013.07.23] {Debugging}
------------------------

Ok, clustering on RFDistance with --editdist 0 should give the same
answers as phybin --bin, but it's not!

On Irene's current Wolbachia dataset, traditional phybin --bin comes
up with this as its top cluster:

    RAxML_bipartitions.88
    RAxML_bipartitions.506
    RAxML_bipartitions.476
    RAxML_bipartitions.39
    RAxML_bipartitions.386
    RAxML_bipartitions.321
    RAxML_bipartitions.3
    RAxML_bipartitions.276
    RAxML_bipartitions.260
    RAxML_bipartitions.256
    RAxML_bipartitions.233
    RAxML_bipartitions.197
    RAxML_bipartitions.171
    RAxML_bipartitions.165
    RAxML_bipartitions.141
    RAxML_bipartitions.116

By contrast, --UPGMA --editdist0 comes up with this:

    RAxML_bipartitions.120
    RAxML_bipartitions.167
    RAxML_bipartitions.250
    RAxML_bipartitions.38
    RAxML_bipartitions.444
    RAxML_bipartitions.494
    
This is with no branch-collapsing...  Ugh, that's no ovelap at all.
If we look at the second bin for --UPGMA though, we see overlap:

    RAxML_bipartitions.233
    RAxML_bipartitions.3
    RAxML_bipartitions.321
    RAxML_bipartitions.39
    RAxML_bipartitions.476
    RAxML_bipartitions.88

Ok, at least one is a subset of the other this time... and here's the diff:

    RAxML_bipartitions.116
    RAxML_bipartitions.141
    RAxML_bipartitions.165
    RAxML_bipartitions.171
    RAxML_bipartitions.197
    RAxML_bipartitions.256
    RAxML_bipartitions.260
    RAxML_bipartitions.276
    RAxML_bipartitions.386
    RAxML_bipartitions.506

And here's the tree that phybin --bin produced:

   ((14, 3_), (19, (5_, 13)), ((1_, 2_), (7_, (18, 6_))));

I singled out 116 and 233, and then printed their normalized forms
using "phybin -p 2 --printnorms tests/t2/*":

    Tree "116"
    ((1_, 2_), (7_, (18, 6_)), ((14, 3_), (19, (13, 5_))))
    Tree "233"
    ((1_, 2_), (7_, (18, 6_)), ((14, 3_), (19, (13, 5_))))

And UNNORMALIZED:

    Tree "116"
    (3_, ((((6_, 18), 7_), (2_, 1_)), (19, (13, 5_))), 14)
    Tree "233"
    (5_, (19, ((3_, 14), ((2_, 1_), (7_, (6_, 18))))), 13)

Wait.... this makes it look like normalization is WRONG!!!
Tree 233 starts with 13+5 in a three-way split with a big chunk
containing everything else.

In fact, the RFDistance is three, which is probably right.

[2013.07.23] {Updating set reps, timing}
----------------------------------------

Ok, went ahead and switched over to unboxed bitvectors for
DenseLabelSet.  It's 10X slower!

The BITVEC_BIPS option is getting the right answer now though.
Apparently the normalization bug just fixed was the only bug.

For 100 trees, 150 taxa (from the hashrf suite):

  * IntSet takes 2.2s, 7.2G alloc, 114M copied.
  * bit vec takes 21s and 69G alloc, 264M copied.
  
Productivities are high in both cases.

Oops... there was a bug in bipSize.  I'm surprised it was getting the
right answer.  Fixing that brought it down to 39G alloc, 12.3s.
Stil, bit vectors are NOT a win for the direct/simple N^2 version.


[2013.07.25] {Comparison ogainst HashRF}
----------------------------------------

Note that HashRF doesn't have very good error messages.  For example,
if given our messy data with sum trees missing taxa:

    ./hashrf Wolbachia/all_trees.tr 928

    *** Collecting the taxon labels ***
	Number of taxa = 10

    *** Reading tree file and collecting bipartitions ***
    libc++abi.dylib: terminate called throwing an exception
    Abort trap: 6

In fact... even if I use phybin to prune out the 503 trees with 10
taxa, I still get the error.


