#!/bin/sh -f

python3 ./site_fullrun.py \
      --site SC-SUMTER --sitegroup USRgroup --caseidprefix OLMT \
      --nyears_ad_spinup 200 --nyears_final_spinup 600 --nyears_transient 130 --tstep 1 \
      --machine cades --compiler gnu --mpilib mpi-serial \
      --cpl_bypass --gswp3 --daymet4 \
      --model_root $HOME/models/E3SM \
      --caseroot $HOME/cases \
      --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
      --spinup_vars \
      --nopointdata \
      --metdir /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_sc-sumter \
      --domainfile /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_sc-sumter/domain.nc \
      --surffile /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_sc-sumter/surfdata.nc \
      --landusefile /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_sc-sumter/surfdata.pftdyn.nc


