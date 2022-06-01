#!/bin/csh -f

#Submission script suitable for high-throughput calculations

set origLoc=`pwd`
set dirControl="$origLoc/"
set maxDepth=5
set add="*/"
set rnCount=0

#make dirControl in folder to avoid early termination
set z=0
cd $dirControl
mkdir depthControl
cd depthControl
while ($z < $maxDepth)
        mkdir $z
        cd $z
        @ z++
end

set i=0
while ($i < $maxDepth)

	foreach j ($dirControl)
		cd j >~/scratch/uselessFile
	
		if (-f thisPoscarWithVac) then

#do stuff here------------------------------------------------------------------------------
	echo "Working on file $rnCount"
	@ rnCount++

	set wd=`echo "$j" | awk -F "/" '{print $(NF-1)}'`
		
	$origLoc/desc -i thisPoscarWithVac -o thisOut		
	grep "CSV" thisOut >thisOut_
	sed -i 's/[^ ]* //' thisOut_
	set descri=`more thisOut_`	
	
	echo "`pwd`,$wd,$descri" >>$origLoc/allDescriptions.csv

	rm thisOut_ thisOut
#stop doing stuff here----------------------------------------------------------------------
		
		cd ../
		endif
	end

set dirControl="$dirControl$add"
@ i++

end

cd $origLoc
rm -rf depthControl
echo "Script Finished"
exit 0
