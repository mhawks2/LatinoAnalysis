# Before Starting

If you are here, you probably want to use the `Latinos` framework to analyze ultra-legacy  Run 2 data. Please consider that this version of the framework is obsolete, and we suggest to use instead [mkShapesRDF](https://github.com/latinos/mkShapesRDF/tree/master).

In case you really need to use this version of the framework, be aware that it is based on a `CMSSW` version not running on `el9`. To use the framework on lxplus, follow the instructions [here](https://gitlab.cern.ch/cms-cat/cmssw-lxplus). In particular, the framework can run only inside a singularity, and you need to create a script called `start_el7.sh`, containing:

    #!/bin/bash
    export APPTAINER_BINDPATH=/afs,/cvmfs,/cvmfs/grid.cern.ch/etc/grid-security:/etc/grid-security,/cvmfs/grid.cern.ch/etc/grid-security/vomses:/etc/vomses,/eos,/etc/pki/ca-trust,/etc/tnsnames.ora,/run/user,/tmp,/var/run/user,/etc/sysconfig,/etc:/orig/etc
schedd=`myschedd show -j | jq .currentschedd | tr -d '"'`

    apptainer -s exec /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cat/cmssw-lxplus/cmssw-el7-lxplus:latest/ sh -c "source /app/setupCondor.sh && export _condor_SCHEDD_HOST=$schedd && export _condor_SCHEDD_NAME=$schedd && export _condor_CREDD_HOST=$schedd && /bin/bash  "


Then, make the script executable using `chmod +x start_el7.sh`, and run it:

    ./start_el7.sh


Good luck!

# Install

Download the framework:

    cmsrel CMSSW_10_6_28
    cd CMSSW_10_6_28/src/
    cmsenv
    git clone --branch 13TeV git@github.com:latinos/setup.git LatinosSetup

Before running setup, edit `scripts/bootstrap.sh` and replace:

    git clone git@github.com:latinos/LatinoAnalysis.git LatinoAnalysis
    cd LatinoAnalysis
    git checkout UL_production
    
with the url for this repo `mhawks2/LatinoAnalysis` and branch `azhPostProc`:

    git clone git@github.com:mhawks2/LatinoAnalysis.git LatinoAnalysis
    cd LatinoAnalysis
    git checkout azhPostProc

then run the setup script:

    source LatinosSetup/SetupShapeOnly.sh
    scram b -j 10

Now we are in the `correctionlib/pybind11/` directory. Update correctionlib to `v2.9.2`:

    git checkout v2.9.2
    cd ../..
    scram b -j 10

Edit the following python files to specify your main directories, i.e. the directories in which your job related information and output will be stored:

    LatinoAnalysis/Tools/python/userConfig.py #for baseDir, jobDir, workDir
    NanoGardener/python/framework/Sites_cfg.py #for xrootdPath, treeBaseDir


# Latino trees post-processing

### Postprocessing script
The mkPostProc.py script is provided, that automates the submission of a full postprocessing campaing. The basic idea is that this script creates one python executable similar to the example quoted above (https://github.com/latinos/LatinoAnalysis/blob/master/NanoGardener/test/postproc.py), with automated definition of the input and output files and the list of modules to be run.

This script is based on three master configuration files:

   * `Sites_cfg.py` (https://github.com/latinos/LatinoAnalysis/blob/master/NanoGardener/python/framework/Sites_cfg.py) defines the sites on which one is willing to write the output. By default, if the postprocessing is run from one of these sites, the output will go to that site.
   * `Productions_cfg.py` (https://github.com/latinos/LatinoAnalysis/blob/master/NanoGardener/python/framework/Productions_cfg.py) Defines the path to the list of samples.
   * `Steps_cfg.py` (https://github.com/latinos/LatinoAnalysis/blob/master/NanoGardener/python/framework/Steps_cfg.py) defines the different steps and the chains of steps to be run.
   
Sample names and paths are found in `NanoGardener/python/framework/samples`, e.g. for 2017 the filename is `Summer20UL17_106x_nAODv9.py`

 Examples:
 
    mkPostProc.py -p Summer20UL17_106ix_nAODv9_Full2017v9 -i MCl1loose2017v9 -s MCCorr2017v9NoJERInHorn -T TWZ_thad_Wlep-DR1 -b -Q nextweek 
 
 this will submit the `MCCorr2017v9NoJERInHorn` chain on the `TWZ_thad_Wlep-DR1` sample defined for the production version `Summer20UL17_106x_nAODv9_Full2017v9`.
 
 Options:
     
         -i : step to start from [default is 'Prod' mode] 
         -s : step to run 
         -b : submit to batch [default is interactive execution] 
         -n : dry-run  just produce script in job directory but do not submit  
         -T <sample1>, ... ,< sampleN > : run only on these samples 
         -E <sample1>, ... ,< sampleN > : do not run on these samples 
         -R : redo all jobs even if output file exist 
         -Q < queuename > : specify queue like 8nh [default btw, see  Site_cfg.py ],   
         Not needed by default 
         --sitescfg  <File> : alternative site cfg
         --modcfg <File> : alternative step/module  cfg
         --datacfg <File> : alternative production cfg


Steps for full postprocessing of nominal samples are `MCl1loose2017v9 -> MCCorr2017v9NoJERInHorn -> l2tightOR2017v9`. The commands to run the chain are shown for a single TWZ sample for 2017:

    mkPostProc.py -p Summer20UL17_106ix_nAODv9_Full2017v9 -s MCl1loose2017v9 -T TWZ_thad_Wlep-DR1 -b -Q nextweek 
    mkPostProc.py -p Summer20UL17_106ix_nAODv9_Full2017v9 -i MCl1loose2017v9 -s MCCorr2017v9NoJERInHorn -T TWZ_thad_Wlep-DR1 -b -Q nextweek 
    mkPostProc.py -p Summer20UL17_106ix_nAODv9_Full2017v9 -i MCl1loose2017v9__MCCorr2017v9NoJERInHorn -s l2tightOR2017v9 -T TWZ_thad_Wlep-DR1 -b -Q nextweek 

### Systematics (Up/Down variations)

Systematics are run from the last step of the base chain above. Up and Down variations are run separately. For example:

    mkPostProc.py -p Summer20UL17_106ix_nAODv9_Full2017v9 -i MCl1loose2017v9__MCCorr2017v9NoJERInHorn__l2tightOR2017v9 -s ElepTup_suffix -T TWZ_thad_Wlep-DR1 -b -Q nextweek 

The list of systematics steps for the Up variations are `ElepTup_suffix, MupTup_suffix, METup_suffix, JERup_suffix`. The Down variation steps are identical but swapping `up_suffix` for `do_suffix`.

### NOTE

JES systematics for the Run 2 UL postprocessing campaign were run with mkShapesRDF, using the `jes-production` branch [here](https://github.com/latinos/mkShapesRDF/tree/jes-production)

Installation instructions are found in the link above.

Update with AZH-specific instructions for following the Run 2 UL prescription are coming soon! 



