#!/bin/sh

#unpacks directory with pdb files 
gunzip pdbFiles.tar.gz
tar -xvf pdbFiles.tar
#rm pdbFiles.tar
#creates the csv file that will hold the molprobity data
if [ ! -f "finalrefinementdata.csv" ]; then
		header="Gene Name,Residue Range,Repository,Full Protein Length,% Sequence Identity,% Disordered,% Ramachandran Favored,% Ramachandran Outliers,Poor Rotamers,Clash Score,MolProbity Score,Energy After FFX Minimization 0.8 (Kcal/mol),Energy After FFX Minimization 0.1 (Kcal/mol),FFX % Ramachandran Favored,FFX % Ramachandran Outliers,FFX Poor Rotamers,FFX Clash Score,FFX MolProbity Score,Final Energy (Kcal/mol)"
        echo $header > finalrefinementdata.csv
fi
#enters into directory that holds pdb files
cd pdbFiles
#stores each pdb file directory name into an array called "gene"
shopt -s nullglob
gene=(*/)
i=0
#while loops through each directory
while [  $i -lt ${#gene[@]} ]; do
	#removes the "/" from the gene directory name
	gene[i]=${gene[i]%/}
	cd ${gene[i]}
		#adds each residue range to an array called "residue"
		shopt -s nullglob
		residue=(*/)
		x=0
		while [ $x -lt ${#residue[@]} ]; do
			residue[x]=${residue[x]%/}
			cd ${residue[x]}
			#runs phenix.molprobity on the pdb file
			phenix.molprobity ${gene[i]}_${residue[x]}.pdb
			rm molprobity_coot.py
			#stores the molprobity data (found in molprobity.out) to the file Original_MolProbity.txt
			tee Original_MolProbity.txt <<<"$(tail molprobity.out)"
			rm molprobity.out
			#copies the properties, minimize, and manybody job files to the present directory
			cp ../../../jobFiles/minimize.properties .
			cp ../../../jobFiles/minimize.job .
			cp ../../../jobFiles/secondminimize.job .
			cp ../../../jobFiles/rotamer.job .
			cp ../../../jobFiles/finalminimize.job .
			#renames minimize.properties to GENE_RESIDUERANGE.properties
			mv minimize.properties ${gene[i]}_${residue[x]}.properties
			cp ${gene[i]}_${residue[x]}.properties extraMinimize.properties
			#replaces "fileName" in each job file with "GENE_RESIDUERANGE"
			sed -i "s/fileName/${gene[i]}_${residue[x]}/g" minimize.job
			sed -i "s/fileName/${gene[i]}_${residue[x]}/g" secondminimize.job
			sed -i "s/fileName/${gene[i]}_${residue[x]}/g" rotamer.job
			sed -i "s/fileName/${gene[i]}_${residue[x]}/g" finalminimize.job
			#submits minimize.job
			minimize_jid=$(qsub -terse minimize.job)
			secondminimize_jid=$(qsub -terse -hold_jid $minimize_jid secondminimize.job)
			rotamer_jid=$(qsub -terse -hold_jid $minimize_jid rotamer.job)
			qsub -hold_jid $rotamer_jid finalminimize.job
			#echo "this is where the Minimize job would submit for ${gene[i]}_${residue[x]}.pdb"
			#exits the current directory so the loop can continue
			cd ../
			let x=x+1
		done
		cd ../
	let i=i+1
done