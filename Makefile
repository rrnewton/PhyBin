


default: haskell

all: find_tree_matches.exe haskell

haskell:
	ghc -O2 --make Main.hs -o phybin.exe
	strip phybin.exe

win:
	ghc -DWIN32 --make Bio/Phylogeny/PhyBin/Main.hs -o phybin.exe

release: haskell
	upx phybin.exe

find_tree_matches.exe: find_tree_matches.ss
	mzc --exe find_tree_matches.exe find_tree_matches.ss

clean:
	rm -f phybin *.exe *.hi *.o Bio/Phylogeny/PhyBin/*.o Bio/Phylogeny/PhyBin/*.hi

################################################################################

#launch: clean release doc distro push upload
launch: doc distro push upload

doc:
	pandoc -s -S --toc -c website.css README.txt -o website/index.html

distro:
	cabal sdist
	cp -v dist/*.tar.gz website/

# HASKELLSRV=haskell.intel
HASKELLSRV=community.haskell.org

upload:
	rsync -vrplt website/ rrnewton@$(HASKELLSRV):/srv/code/phybin/; ssh rrnewton@$(HASKELLSRV) ./setperms.sh
#	ssh rrnewton@$(HASKELLSRV) chmod o+r -R '/srv/code/phybin/*'
#	scp -pr website/* rrnewton@$(HASKELLSRV):/srv/code/phybin/
#	scp -pr website rrnewton@haskell.intel:

push:
#	hg push ssh://rrnewton@community.haskell.org/phybin/repo/ --remotecmd '~/bin/hg'
# Problem: Different versions of darcs right now
#	darcs push -a rrnewton@community.haskell.org:phybin/repo/
	darcs push -a ./website/repo/





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

#------------------------------------------------------------
rickettsia:    $(LOCALSETS)/Rickettsia/renaming_table.txt
	rm -rf $(LOCALSETS)/Rickettsia/phybin_outputs
	mkdir  $(LOCALSETS)/Rickettsia/phybin_outputs
	./phybin.exe -g -n 15 -m $(LOCALSETS)/Rickettsia/renaming_table.txt -s '_'  -o $(LOCALSETS)/Rickettsia/phybin_output/ $(LOCALSETS)/Rickettsia/final_trees/*BranchLab*.out 

$(LOCALSETS)/Rickettsia/renaming_table.txt: $(LOCALSETS)/Rickettsia/Rickettsia_orthololgs.txt
	runghc stripTable.hs $^ > $@

#------------------------------------------------------------
rickettsiales: $(LOCALSETS)/Rickettsiales/renaming_table.txt
	rm -rf $(LOCALSETS)/Rickettsiales/phybin_outputs
	mkdir  $(LOCALSETS)/Rickettsiales/phybin_outputs
	./phybin.exe -g -n 29 -m $(LOCALSETS)/Rickettsiales/renaming_table.txt -s '_'  -o $(LOCALSETS)/Rickettsiales/phybin_output/ $(LOCALSETS)/Rickettsiales/final_trees/*BranchLab*.out 

$(LOCALSETS)/Rickettsiales/renaming_table.txt: $(LOCALSETS)/Rickettsiales/Rickettsiales_orthologs.txt
	runghc stripTable.hs $^ > $@

#------------------------------------------------------------
wolbachia:     $(LOCALSETS)/Wolbachia/renaming_table.txt
	rm -rf $(LOCALSETS)/Wolbachia/phybin_outputs
	mkdir  $(LOCALSETS)/Wolbachia/phybin_outputs
	./phybin.exe -g -n 4 -m $(LOCALSETS)/Wolbachia/renaming_table.txt -s '_'  -o $(LOCALSETS)/Wolbachia/phybin_output/ $(LOCALSETS)/Wolbachia/final_trees/*BranchLab*.out 

$(LOCALSETS)/Wolbachia/renaming_table.txt: $(LOCALSETS)/Wolbachia/Wolbachia_orthologs.txt
	runghc stripTable.hs $^ > $@


temp:
	./phybin.exe ~/newton_and_newton_local/datasets/rhizobia/gde/*.dnd -o temp -p 2 -n 7 

