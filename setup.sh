## based on Latinos setup script 'SetupShapeOnly.sh' 
## from https://github.com/latinos/setup.git on branch '13TeV'

if [ -z $CMSSW_BASE ]; then
    echo "======================================="
    echo "No CMS environment detected; stopping..."
    echo "======================================="
    exit 1
fi

source $CMSSW_BASE/src/LatinoAnalysis/functions.sh

if [[ "$CMSSW_VERSION" == CMSSW_10_*_* ]]; then
    echo "======================================="
    echo "running with $CMSSW_VERSION - this is a 13 TeV setup!"
    echo "Current time:" $(date)
    echo "checking out additional repositories; this could take a while ..."
    echo "======================================="

    cd $CMSSW_BASE/src    

    #echo "++++++++++++ WARNING: Assuming UL setup: using UL_production branch +++++++++++++"
    #
    #echo " - Basic Code"
    #
    #git clone git@github.com:latinos/LatinoAnalysis.git LatinoAnalysis
    #cd LatinoAnalysis
    #git checkout UL_production
    #cd -

    echo " - Nano Tools"

    git clone git@github.com:cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
    
    cp LatinoAnalysis/Tools/data/JECs/*txt PhysicsTools/NanoAODTools/data/jme/
    cp LatinoAnalysis/Tools/data/JECs/*tgz PhysicsTools/NanoAODTools/data/jme/
    cp LatinoAnalysis/Tools/data/JERs/Summer19UL17_JRV2_MC.tgz PhysicsTools/NanoAODTools/data/jme/
    cp LatinoAnalysis/NanoGardener/python/data/Summer16_25nsV1b_MC.tgz PhysicsTools/NanoAODTools/data/jme/ 
    cp LatinoAnalysis/NanoGardener/python/data/Fall17_V3b_MC.tgz PhysicsTools/NanoAODTools/data/jme/
    cp LatinoAnalysis/NanoGardener/python/data/Autumn18_V7b_MC.tgz PhysicsTools/NanoAODTools/data/jme/
    

    echo " - Plotting Tools"

    git clone git@github.com:yiiyama/multidraw.git LatinoAnalysis/MultiDraw
    cd LatinoAnalysis/MultiDraw
    git checkout 2.0.12 2>/dev/null
    ./mkLinkDef.py --cmssw
    cd ../..

    echo " - MELA new version"

    git clone git@github.com:MELALabs/MelaAnalytics.git MelaAnalytics
    cd MelaAnalytics ; git checkout -b from-v22 v2.2 ; cd ..
    git clone https://github.com/JHUGen/JHUGenMela.git JHUGenMELA
    cd JHUGenMELA; git checkout -b from-v235 v2.3.5 ; source setup.sh -j 12 ; cd ..


    scram b -j 8

    echo " - Correction Lib"
    git clone ssh://git@gitlab.cern.ch:7999/cms-crossPOG/jsonpog-integration.git
    git clone --recursive git@github.com:cms-nanoAOD/correctionlib.git
    cd correctionlib
    echo "   - enabling python2 bindings"
    cd pybind11
    git checkout v2.9.2
    cd ..
    make PYTHON=python
    make install
    cp -R correctionlib $CMSSW_BASE/src/LatinoAnalysis/NanoGardener/python/modules
    cd ..

    echo " - extras"
    cp LatinoAnalysis/extras/* $CMSSW_BASE/src/

else
    echo "======================================="
    echo "You are using release $CMSSW_VERSION which is not supported by this script."
    echo "A CMSSW_10_*_* release is required for postprocessing UL samples."
    echo "======================================="

fi
