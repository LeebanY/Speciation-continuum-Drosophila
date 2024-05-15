#!/bin/bash/

for all in *;do
dxy_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "D_xy" | sed 's/  */ /g' | cut -d' ' -f4);
sample_A_WG=$( less $all/calls/Summarystats.wholegenome.info.txt | grep "A = " | cut -d' ' -f7);
sample_B_WG=$( less $all/calls/Summarystats.wholegenome.info.txt | grep "B = " | cut -d' ' -f7);
totalblocks_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "total blocks" | sed 's/  */ /g' | cut -d' ' -f5);
interval_cov_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "interval coverage" | sed 's/  */ /g' | cut -d' ' -f5);
invar_blocks_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "invariant blocks" | sed 's/  */ /g' | cut -d' ' -f5);
four_gam_vio_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "four-gamete-violation blocks" | sed 's/  */ /g' | cut -d' ' -f5);
hetA_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "heterozygosity (A)" | sed 's/  */ /g' | cut -d' ' -f5);
hetB_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "heterozygosity (B)" | sed 's/  */ /g' | cut -d' ' -f5);
interval_length_WG=$(less $all/calls/Summarystats.wholegenome.info.txt | grep "intervals" | cut -d' ' -f5);
dxy_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "D_xy" | sed 's/  */ /g' | cut -d' ' -f4);
sample_A_IN=$( less $all/calls/Summarystats.introns.info.txt | grep "A = " | cut -d' ' -f7);
sample_B_IN=$( less $all/calls/Summarystats.introns.info.txt | grep "B = " | cut -d' ' -f7);
totalblocks_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "total blocks" | sed 's/  */ /g' | cut -d' ' -f5);
interval_cov_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "interval coverage" | sed 's/  */ /g' | cut -d' ' -f5);
invar_blocks_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "invariant blocks" | sed 's/  */ /g' | cut -d' ' -f5);
four_gam_vio_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "four-gamete-violation blocks" | sed 's/  */ /g' | cut -d' ' -f5);
hetA_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "heterozygosity (A)" | sed 's/  */ /g' | cut -d' ' -f5);
hetB_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "heterozygosity (B)" | sed 's/  */ /g' | cut -d' ' -f5);
interval_length_IN=$(less $all/calls/Summarystats.introns.info.txt | grep "intervals" | cut -d' ' -f5);
echo $all $dxy_WG $sample_A_WG $sample_B_WG $hetA_WG $hetB_WG $totalblocks_WG $interval_cov_WG $interval_length_WG $four_gam_vio_WG $invar_blocks_WG \
$dxy_IN $sample_A_IN $sample_B_IN $hetA_IN $hetB_IN $totalblocks_IN $interval_cov_IN $interval_length_IN $four_gam_vio_IN $invar_blocks_IN;
done
