#!/bin/csh -f

set origLoc=`pwd`
set dirControl="$origLoc/"
set maxDepth=5
set add="*/"
set rnCount=0

#make dirControl in folder to avoid early termination
set z=0
cd $dirControl
mkdir depthControl__
cd depthControl__
while ($z < $maxDepth)
        mkdir $z
        cd $z
        @ z++
end

set i=0
while ($i < $maxDepth)

	foreach j ($dirControl)
		cd j >$origLoc/.uselessFile 
	
		if (-f POSCAR0) then

#do stuff here------------------------------------------------------------------------------

set wd=`echo "$j" | awk -F "/" '{print $(NF-1)}'`
mkdir $origLoc/move/$wd
cp POSCAR0 $origLoc/move/$wd/.
cp CONTCAR $origLoc/move/$wd/.

#stop doing stuff here----------------------------------------------------------------------
		
		cd ../
		endif
	end

set dirControl="$dirControl$add"
@ i++

end

cd $origLoc
rm -rf depthControl__
echo "Script Finished"
exit 0
