imgtSearch="/work/pi_gblanck/old_blanck_group/shared/imgtSearchIgTcr.php"
searchScript="/work/pi_gblanck/Arpan/NBL_WXS/scripts/t2_SearchReadsFix.sh"

#Set bin paths
tsv2xlsx="/work/pi_gblanck/old_blanck_group/shared/tsv2xlsx.py"
#parallel --citation
phpBin="/work/pi_gblanck/old_blanck_group/shared/bin/php"
pigz="/work/pi_gblanck/old_blanck_group/shared/bin/pigz"

# this is a sort of manager script that takes a bunch of tasks and then executes them.

# no configuration below this line !!!

module purge
module add apps/python/2.7.11
module load apps/samtools/1.3.1

samplesName="${bamFolder##*/}"

#Start the processing
echo "Starting up..."
echo `date +"%a %x %X"`

#Make results folder
resultsDir="${bamFolder}_Results"
echo bamFolder is ${bamFolder}
#resultsDir="${HOME}/${samplesName}_Results"


##############################################################
if [[ "${toDo[*]}" == *quickcheck* ]]; then
	echo "Running quickcheck..."
	echo `date +"%a %x %X"`
	numBams=`find $bamFolder -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam" -printf '.' | wc -c`
	numBais=`find $bamFolder -type f -name "*.bai" -printf '.' | wc -c`
	
	echo "numBams = $numBams"
	echo "numBais = $numBais"

	if [ "$numBams" -ne "$numBais" ]; then
		echo "Warning: Some index files are missing. Indexing files..."
		echo `date +"%a %x %X"`

		fileListM=`find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam"`
		
		for bamFile in $fileListM
		do
			baseName=${bamFile%.bam}
			baiName="${baseName}.bai"

			if [ ! -e "$baiName" ]; then
				echo "Indexing ${bamFile}..."
				samtools index $bamFile
			fi
		done
	fi

	
	fileList=`find $bamFolder -type f -name "*.bam"`
	for bamFile in $fileList
	do
		samtools quickcheck -v $bamFile
	done
fi


##############################################################

#Make the results folders for each receptor
echo "Making folders..."
echo `date +"%a %x %X"`

mkdir -p $resultsDir

for IgTcr in "${IgTcrs[@]}"
do
	resultsDirIgTcr="${resultsDir}/$IgTcr"

	mkdir -p $resultsDirIgTcr
done


##############################################################

echo "Fetching file list..."
echo `date +"%a %x %X"`

#Get list of file paths
#Check if Ig/TCR array requires unmapped reads
if [[ "${toDo[*]}" == *extractUnmapped* ]]; then
	if [[ "${IgTcrs[*]}" == *UM* ]]; then
		numBams=`find $bamFolder -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam" -printf '.' | wc -c`
		numUM=`find $bamFolder -type f -name "*unmapped.bam" -printf '.' | wc -c`
		echo "Found $numBams bam files and $numUM unmapped bam files..."

		# # there isn't an unmapped file for each bam file, so we need to extract the unmapped regions and make these *_unmapped.bam files
		# if [ "$numBams" -ne "$numUM" ]; then
		echo "Extracting unmapped regions..."
		fileListM=$(find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam")
		echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -f4 {} -h -o {.}_unmapped.bam"
		# fi
	fi
fi

##############################################################

if [[ "${toDo[*]}" == *extractReceptors* ]]; then
	echo "Extracting receptor regions..."
	echo `date +"%a %x %X"`
	regionTRA="${chr}14:20000000-24000000"
	regionTRB="${chr}7:140000000-145000000"
	regionTRD="${chr}14:20000000-25000000"
	regionTRG="${chr}7:36000000-41000000"
	regionIGH="${chr}14:103000000-107349000"
	regionIGK="${chr}2:86000000-92000000"
	regionIGL="${chr}22:21000000-25000000"

	fileListM=$(find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam")
	echo "Extracting TRA..."
	echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -h {} $regionTRA -o {.}_TRA.bam"
	echo "Extracting TRB..."
	echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -h {} $regionTRB -o {.}_TRB.bam"
	echo "Extracting TRD..."
	echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -h {} $regionTRD -o {.}_TRD.bam"
	echo "Extracting TRG..."
	echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -h {} $regionTRG -o {.}_TRG.bam"
	echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -h {} $regionIGH -o {.}_IGH.bam"
	echo "Extracting IGK..."
	echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -h {} $regionIGK -o {.}_IGK.bam"
	echo "Extracting IGL..."
	echo "$fileListM" | parallel -j${numCores} --eta  "echo {} ; samtools view -b -h {} $regionIGL -o {.}_IGL.bam"
fi

##############################################################
if [[ "${toDo[*]}" == *runSearches* ]]; then
	echo "Running Searches..."
	echo `date +"%a %x %X"`

	# unmapped reads
	fileListUM=$(find ${bamFolder} -type f -name "*unmapped.bam")

	#mapped reads corresponding to the 7 immune receptors
	fileListM=$(find ${bamFolder} -type f -name "*.bam" ! -name "*unmapped*"  ! -name "*TRA.bam" ! -name "*TRB.bam" ! -name "*TRD.bam" ! -name "*TRG.bam" ! -name "*IGH.bam" ! -name "*IGK.bam" ! -name "*IGL.bam")

	for IgTcr in "${IgTcrs[@]}"
	do
		echo "Processing $IgTcr..."
		echo `date +"%a %x %X"`
		if [[ "$IgTcr" == *UM* ]]; then
			echo "$fileListUM" | parallel -j${numCores} --eta "echo {} ; sh $searchScript samtools {} ${resultsDir}/${IgTcr}/{/}.tsv ${IgTcr} ${chr}"
		else
			echo "$fileListM" | parallel -j${numCores} --eta "echo {} ; sh $searchScript samtools {} ${resultsDir}/${IgTcr}/{/}.tsv ${IgTcr} ${chr}"
		fi
	done
fi
##############################################################

if [[ "${toDo[*]}" == *tsv2xlsx* ]]; then
	echo "Running tsv2xlsx..."
	echo `date +"%a %x %X"`

	for IgTcr in "${IgTcrs[@]}"
	do
		echo "Running tsv2xlsx ${IgTcr}..."
		echo `date +"%a %x %X"`

		mkdir -p "${resultsDir}_xlsx/$IgTcr"

		if [[ "$IgTcr" == *UM* ]]; then
			echo "$fileListUM" | parallel -j${numCores} --eta "python $tsv2xlsx $resultsDir/${IgTcr}/{/}.tsv ${resultsDir}_xlsx/${IgTcr}/{/}.xlsx"
		else
			echo "$fileListM" | parallel -j${numCores} --eta "python $tsv2xlsx $resultsDir/${IgTcr}/{/}.tsv ${resultsDir}_xlsx/${IgTcr}/{/}.xlsx"
		fi
	done
fi

##############################################################
if [[ "${toDo[*]}" == *imgtSearch* ]]; then
	echo "Running IMGT Search..."
	echo `date +"%a %x %X"`

	for IgTcr in "${IgTcrs[@]}"
	do
		echo "Running IMGT for ${IgTcr}..."
		echo `date +"%a %x %X"`
		$phpBin $imgtSearch $resultsDir/${IgTcr} $resultsDir/${samplesName}_${IgTcr}_vjMatchList.tsv ${IgTcr}
		if [[ "${toDo[*]}" == *tsv2xlsx* ]]; then
			python $tsv2xlsx $resultsDir/${samplesName}_${IgTcr}_vjMatchList.tsv ${resultsDir}_xlsx/${samplesName}_${IgTcr}_vjMatchList.xlsx
		fi
	done
fi

##############################################################

if [[ "${toDo[*]}" == *pigz* ]]; then
	tar -I pigz -cf ${resultsDir}.tar.gz ${resultsDir}
	if [[ "${toDo[*]}" == *tsv2xlsx* ]]; then
		tar -I pigz -cf ${resultsDir}_xlsx.tar.gz ${resultsDir}_xlsx
	fi
fi


##############################################################
echo "Finished."
echo `date +"%a %x %X"`