#!/bin/bash
for i in {1..9}
do
	hipo-utils -merge /volatile/clas12/osg/job_2443/simu_"$i"*/dst.hipo -o /volatile/clas12/matthewn/Simulations/Dibaryon/KpKm_RGA_FALL18_10M_240221_"$i".hipo
done
