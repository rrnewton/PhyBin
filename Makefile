


default: haskell

all: find_tree_matches.exe haskell

# GHC=ghc-7.4.2
GHC=ghc-7.6.3

release: haskell
	upx phybin.exe
	du -h phybin.exe

haskell:
	cabal configure
	cabal build
	cp ./dist/build/phybin/phybin phybin.exe
	du -h phybin.exe
	strip phybin.exe
	du -h phybin.exe

win:
	$(GHC) -DWIN32 --make Main.hs -o phybin.exe

find_tree_matches.exe: find_tree_matches.ss
	mzc --exe find_tree_matches.exe find_tree_matches.ss

clean:
	rm -f phybin *.exe *.hi *.o Bio/Phylogeny/PhyBin/*.o Bio/Phylogeny/PhyBin/*.hi

################################################################################

#launch: clean release doc distro push upload
launch: doc distro push upload

doc:
	pandoc -s -S --toc -c website.css README.md -o website/index.html

release_website:
	cp website/* ~/.hyplan/projects/phybin/

distro:
	cabal sdist
	cp -v dist/*.tar.gz website/

# HASKELLSRV=haskell.intel
HASKELLSRV=community.haskell.org


################################################################################
# Collecting shortcuts for running different datasets:


# LOCALSETS=~/newton_and_newton_local/datasets
LOCALSETS=.

rhizobia:
	rm -rf $(LOCALSETS)/rhizobia/phybin_outputs
	mkdir $(LOCALSETS)/rhizobia/phybin_outputs
	./phybin.exe $(LOCALSETS)/rhizobia/gde/*.dnd -o $(LOCALSETS)/rhizobia/phybin_output/ -p 2 -n 7 -g -v 

#yersinia_test_names:
#	./phybin.exe -w ~/newton_and_newton_local/datasets/yersinia/yersinia_trees/111.dnd -s .  -m ../datasets/yersinia/name_table_hack_yersinia.txt 

yersinia:
	rm -rf $(LOCALSETS)/yersinia/phybin_outputs
	mkdir  $(LOCALSETS)/yersinia/phybin_outputs
	./phybin.exe $(LOCALSETS)/yersinia/yersinia_trees/ -o $(LOCALSETS)/yersinia/phybin_outputs/ -m ../datasets/yersinia/name_table_hack_yersinia.txt -s . -n 9 -g

legionella:
	rm -rf $(LOCALSETS)/legionella/phybin_outputs
	mkdir  $(LOCALSETS)/legionella/phybin_outputs
	./phybin.exe $(LOCALSETS)/legionella/legionella_orthologs_aa/ -o $(LOCALSETS)/legionella/phybin_outputs/ -m ../datasets/legionella/name_table_hack_legionella.txt -s 0123456789 -n 4 -g



# Newer ones [2012.11.19]:
newbatch: rickettsia rickettsiales wolbachia

# BRANCHTHRESH=0.01
BRANCHTHRESH=0.05
PHYBIN= ./phybin.exe -b $(BRANCHTHRESH) -g

#------------------------------------------------------------
rickettsia:    $(LOCALSETS)/Rickettsia/renaming_table.txt
	rm -rf $(LOCALSETS)/Rickettsia/phybin_outputs
	mkdir  $(LOCALSETS)/Rickettsia/phybin_outputs
	$(PHYBIN) -n 15 -m $(LOCALSETS)/Rickettsia/renaming_table.txt -s '_'  -o $(LOCALSETS)/Rickettsia/phybin_output/ $(LOCALSETS)/Rickettsia/final_trees/*BranchLab*.out 

$(LOCALSETS)/Rickettsia/renaming_table.txt: $(LOCALSETS)/Rickettsia/Rickettsia_orthololgs.txt
	runghc stripTable.hs $^ > $@

#------------------------------------------------------------
rickettsiales: $(LOCALSETS)/Rickettsiales/renaming_table.txt
	rm -rf $(LOCALSETS)/Rickettsiales/phybin_outputs
	mkdir  $(LOCALSETS)/Rickettsiales/phybin_outputs
	$(PHYBIN) -n 29 -m $(LOCALSETS)/Rickettsiales/renaming_table.txt -s '_'  -o $(LOCALSETS)/Rickettsiales/phybin_output/ $(LOCALSETS)/Rickettsiales/final_trees/*BranchLab*.out 

$(LOCALSETS)/Rickettsiales/renaming_table.txt: $(LOCALSETS)/Rickettsiales/Rickettsiales_orthologs.txt
	runghc stripTable.hs $^ > $@

#------------------------------------------------------------
wolbachia:     $(LOCALSETS)/Wolbachia/renaming_table.txt
	rm -rf $(LOCALSETS)/Wolbachia/phybin_outputs
	mkdir  $(LOCALSETS)/Wolbachia/phybin_outputs
	$(PHYBIN) -n 4 -m $(LOCALSETS)/Wolbachia/renaming_table.txt -s '_'  -o $(LOCALSETS)/Wolbachia/phybin_output/ $(LOCALSETS)/Wolbachia/final_trees/*BranchLab*.out 

$(LOCALSETS)/Wolbachia/renaming_table.txt: $(LOCALSETS)/Wolbachia/Wolbachia_orthologs.txt
	runghc stripTable.hs $^ > $@



vary_branchlen: $(LOCALSETS)/Rickettsia/renaming_table.txt $(LOCALSETS)/Rickettsiales/renaming_table.txt $(LOCALSETS)/Wolbachia/renaming_table.txt
	./vary_branchlen.sh 15 Rickettsia
	./vary_branchlen.sh 29 Rickettsiales
	./vary_branchlen.sh  4 Wolbachia

temp:
	./phybin.exe ~/newton_and_newton_local/datasets/rhizobia/gde/*.dnd -o temp -p 2 -n 7 

