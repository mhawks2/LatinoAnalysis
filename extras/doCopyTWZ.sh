#! /bin/sh


########################################
##########       2018       ############
########################################

#outDir="/eos/user/m/mihawksw/azh/postprocessing/TWZ_samples_18"
#mkdir $outDir
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_thad_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > twz-thad-wlep-18.txt
#echo $(cat twz-thad-wlep-18.txt | wc -l) files for TWZToLL_thad_Wlep-DR1 2018.....
#while IFS= read -r line; do 
#    fname=$(basename $line)
#    echo $fname
#    env -i  X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_thad_Wlep-DR1/"$fname" &> /dev/null
#done < twz-thad-wlep-18.txt
#
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_tlep_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > twz-tlep-wlep-18.txt
#echo $(cat twz-tlep-wlep-18.txt | wc -l) files for TWZToLL_tlep_Wlep-DR1 2018.....
#while IFS= read -r line; do 
#    fname=$(basename $line)
#    echo $fname
#    env -i X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_tlep_Wlep-DR1/"$fname" &> /dev/null
#done < twz-tlep-wlep-18.txt
#
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_tlep_Whad-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM" > twz-tlep-whad-18.txt
#echo $(cat twz-tlep-whad-18.txt | wc -l) files for TWZToLL_tlep_Whad-DR1 2018.....
#while IFS= read -r line; do 
#    fname=$(basename $line)
#    echo $fname
#    env -i X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_tlep_Whad-DR1/"$fname" &> /dev/null
#done < twz-tlep-whad-18.txt




########################################
##########       2017       ############
########################################

#outDir="/eos/user/m/mihawksw/azh/postprocessing/TWZ_samples_18"
#mkdir $outDir
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_thad_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM" > twz-thad-wlep.txt
#while IFS= read -r line; do 
#    echo $line
#    env -i X509_USER_PROXY=/tmp/x509up_u139558 gfal-copy davs://skynet013.crc.nd.edu:1094/store/mc/RunIISummer20UL17NanoAODv9/TWZToLL_thad_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/30000/"$line" "$outDir"/TWZToLL_thad_Wlep-DR1/"$line" &> /dev/null
#done < twz-thad-wlep.txt
#
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_thad_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM" > twz-tlep-whad.txt
#while IFS= read -r line; do 
#    echo " "
#    echo "$line"
#    #env -i X509_USER_PROXY=/tmp/x509up_u139558 gfal-copy davs://skynet013.crc.nd.edu:1094/store/mc/RunIISummer20UL17NanoAODv9/TWZToLL_tlep_Whad-DR1_TuneCP5_13TeV_amcatnlo-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/30000/"$line" "$outDir"/TWZToLL_tlep_Whad-DR1/"$line"
#done < twz-tlep-whad.txt
#
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_tlep_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v1/NANOAODSIM" > twz-tlep-wlep.txt
#while IFS= read -r j; do 
#    env -i X509_USER_PROXY=/tmp/x509up_u139558 gfal-copy davs://skynet013.crc.nd.edu:1094/store/mc/RunIISummer20UL17NanoAODv9/TWZToLL_tlep_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/30000/"$j" "$outDir"/TWZToLL_tlep_Wlep-DR1/"$j"
#done < <(twz-tlep-wlep.txt)




########################################
##########    2016noHIPM    ############
########################################

outDir="/eos/user/m/mihawksw/azh/postprocessing/TWZ_samples_16noHIPM"
mkdir $outDir

cp /eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/TWZToLL_thad_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/30000/*   $outDir/TWZToLL_thad_Wlep-DR1/
cp /eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/TWZToLL_tlep_Whad-DR1_TuneCP5_13TeV_amcatnlo-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/2820000/* $outDir/TWZToLL_tlep_Whad-DR1/
cp /eos/opendata/cms/mc/RunIISummer20UL16NanoAODv9/TWZToLL_tlep_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/30000/*   $outDir/TWZToLL_tlep_Wlep-DR1/

#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_thad_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM" > twz-thad-wlep-16noHIPM.txt
#echo $(cat twz-thad-wlep-16noHIPM.txt | wc -l) files for TWZToLL_thad_Wlep-DR1 2016noHIPM.....
#while IFS= read -r line; do 
#    fname=$(basename $line)
#    echo $fname
#    env -i  X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_thad_Wlep-DR1/"$fname" &> /dev/null
#done < twz-thad-wlep-16noHIPM.txt
#
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_tlep_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM" > twz-tlep-wlep-16noHIPM.txt
#echo $(cat twz-tlep-wlep-16noHIPM.txt | wc -l) files for TWZToLL_tlep_Wlep-DR1 2016noHIPM.....
#while IFS= read -r line; do 
#    fname=$(basename $line)
#    echo $fname
#    env -i X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_tlep_Wlep-DR1/"$fname" &> /dev/null
#done < twz-tlep-wlep-16noHIPM.txt
#
#
#dasgoclient -query="instance=prod/global file dataset=/TWZToLL_tlep_Whad-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM" > twz-tlep-whad-16noHIPM.txt
#echo $(cat twz-tlep-whad-16noHIPM.txt | wc -l) files for TWZToLL_tlep_Whad-DR1 2016noHIPM.....
#while IFS= read -r line; do 
#    fname=$(basename $line)
#    echo $fname
#    env -i X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_tlep_Whad-DR1/"$fname" &> /dev/null
#done < twz-tlep-whad-16noHIPM.txt




########################################
##########     2016HIPM     ############
########################################

outDir="/eos/user/m/mihawksw/azh/postprocessing/TWZ_samples_16HIPM"
mkdir $outDir

dasgoclient -query="instance=prod/global file dataset=/TWZToLL_thad_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM" > twz-thad-wlep-16HIPM.txt
echo $(cat twz-thad-wlep-16HIPM.txt | wc -l) files for TWZToLL_thad_Wlep-DR1 2016HIPM.....
while IFS= read -r line; do 
    fname=$(basename $line)
    echo $fname
    env -i  X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_thad_Wlep-DR1/"$fname" &> /dev/null
done < twz-thad-wlep-16HIPM.txt


dasgoclient -query="instance=prod/global file dataset=/TWZToLL_tlep_Wlep-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM" > twz-tlep-wlep-16HIPM.txt
echo $(cat twz-tlep-wlep-16HIPM.txt | wc -l) files for TWZToLL_tlep_Wlep-DR1 2016HIPM.....
while IFS= read -r line; do 
    fname=$(basename $line)
    echo $fname
    env -i X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_tlep_Wlep-DR1/"$fname" &> /dev/null
done < twz-tlep-wlep-16HIPM.txt


dasgoclient -query="instance=prod/global file dataset=/TWZToLL_tlep_Whad-DR1_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM" > twz-tlep-whad-16HIPM.txt
echo $(cat twz-tlep-whad-16HIPM.txt | wc -l) files for TWZToLL_tlep_Whad-DR1 2016HIPM.....
while IFS= read -r line; do 
    fname=$(basename $line)
    echo $fname
    env -i X509_USER_PROXY=/tmp/x509up_u"$UID" gfal-copy davs://skynet013.crc.nd.edu:1094"$line" "$outDir"/TWZToLL_tlep_Whad-DR1/"$fname" &> /dev/null
done < twz-tlep-whad-16HIPM.txt

