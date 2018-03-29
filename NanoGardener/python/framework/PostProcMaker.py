#!/usr/bin/env python
import sys, re, os, os.path, math, copy
import string
import subprocess

# configuration auto-loaded where the job directory and the working directory is defined
from LatinoAnalysis.Tools.userConfig  import *

# functions used in everyday life ...
from LatinoAnalysis.Tools.commonTools import *
from LatinoAnalysis.Tools.batchTools  import *


class PostProcMaker():

# ------------- Configuration

   def __init__(self) : 

     self._cmsswBasedir = os.environ["CMSSW_BASE"] 

     self._aaaXrootd = 'root://cms-xrd-global.cern.ch//'

     # root tree prefix
     self._treeFilePrefix= 'nanoLatino_'

     # site
     self._LocalSite    = None
     self._TargetSite   = None

     # Cfg
     self._Sites        = {}
     self._Steps        = {}
     self._Productions  = {}

     # job mode = Interactive / Batch / Crab / DryRun
     self._pretend      = False 
     self._jobMode      = 'Interactive'
     self._batchQueue   = '8nh'

     # What to do
     self._prodList     = []
     self._stepList     = [] 
     self._iniStep      = None 
     self._selTree      = []
     self._excTree      = []
     self._redo         = False
  
     # Samples
     self._Samples     = {}

   def Reset(self) : 

     # Samples
     self._Samples     = {}

   def configSite(self,TargetSite=None):
     osName = os.uname()[1]
     for iSite in self._Sites :
       if iSite in osName : self._LocalSite = iSite
     if self._Sites == None :
       print 'ERROR: Unknown site : ', osName
       exit()
     print '_LocalSite  = ',self._LocalSite
     print '_TargetSite = ',self._TargetSite

   def configBatch(self,queue):
     self._batchQueue = queue

   def readSampleFile(self,iProd):
     prodFile=self._cmsswBasedir+'/src/'+self._Productions[iProd]['samples']
     if os.path.exists(prodFile):
       handle = open(prodFile,'r')
       exec(handle)
       self._Samples     = Samples
       handle.close()  
     keys2del = []
     if len(self._selTree) > 0 :
       for iSample in self._Samples :
         if not iSample in self._selTree : keys2del.append(iSample)
     if len(self._excTree) > 0 :
       for iSample in self._Samples :
         if iSample in self._excTree : keys2del.append(iSample)
     for iSample in keys2del : del self._Samples[iSample]

# -------------- File Handling

   def selectSample(self,iProd,iStep,iSample):
      # From Production
      if     'onlySample' in  self._Productions[iProd]              \
         and len(self._Productions[iProd]['onlySample']) > 0        \
         and not iSample in self._Productions[iProd]['onlySample']  : return False
      if     'excludeSample' in self._Productions[iProd]            \
         and len(self._Productions[iProd]['excludeSample']) > 0     \
         and iSample in self._Productions[iProd]['excludeSample']   : return False
      # From Step
      if     'onlySample' in  self._Steps[iStep]              \
         and len(self._Steps[iStep]['onlySample']) > 0        \
         and not iSample in self._Steps[iStep]['onlySample']  : return False
      if     'excludeSample' in self._Steps[iStep]            \
         and len(self._Steps[iStep]['excludeSample']) > 0     \
         and iSample in self._Steps[iStep]['excludeSample']   : return False
      # ---
      return True

   def getTargetFileDic(self,iProd,iStep,iSample,FileList):
     FileDic = {}

     # fileCmd .... Directory
     fileCmd = self._Sites[self._LocalSite]['lsCmd']+' '+self._targetDir
     # fileCmd .... Files
     if len(FileList) == 1 : fileCmd += self._treeFilePrefix+iSample+'.root'
     else                  : fileCmd += self._treeFilePrefix+iSample+'__part*.root'

     # fileCmd .... Exec
     proc=subprocess.Popen(fileCmd, stderr = subprocess.PIPE,stdout = subprocess.PIPE, shell = True)
     out, err = proc.communicate()
     FileExistList=string.split(out)
     # Now Check 
     toSkip=[]
     if not self._redo : 
       if  len(FileExistList) == len(FileList) : return FileDic
       for iFile in FileExistList: toSkip.append(iFile.replace('.root','').split('__part')[1])

     if not self._iniStep == 'Prod' :
       for iFile in FileList : 
         iPart = iFile.replace('.root','').split('__part')[1]
         if not iPart in toSkip : 
           FileDic[iFile] = self._targetDir+os.path.basename(iFile)
     else :
       # Here I have to assume/fix the ordering !!!!
       iPart = 0
       for iFile in FileList:
         if not str(iPart) in toSkip :
           PartName=''
           if len(FileList)>0 : PartName='__part'+str(iPart)
           fileTargetName = self._targetDir+self._treeFilePrefix+iSample+PartName+'.root'
           FileDic[self._aaaXrootd+iFile] = fileTargetName
         iPart +=1 

     return FileDic

   def getTargetFiles(self,iProd,iStep):

     self._targetDic = {}

     for iSample in self._Samples :
       if self.selectSample(iProd,iStep,iSample) :
         # From central nanoAOD 
         if self._iniStep == 'Prod' : 
           FileDic = self.getTargetFileDic(iProd,iStep,iSample,self.getFilesFromDAS(self._Samples[iSample]['nanoAOD']))
         else :
           FileDic = self.getTargetFileDic(iProd,iStep,iSample,getSampleFiles(self._sourceDir,iSample,True,'nanoLatino_',True))
         if len(FileDic) : self._targetDic[iSample] = FileDic

   def getFilesFromDAS(self,dataset):
     dasCmd='dasgoclient -query="file dataset='+dataset+'"'
     proc=subprocess.Popen(dasCmd, stderr = subprocess.PIPE,stdout = subprocess.PIPE, shell = True)
     out, err = proc.communicate()
     FileList=string.split(out)
     return FileList

   def mkFileDir(self,iProd,iStep):

     self._targetDir = None
     self._sourceDir = None
     if not self._iniStep == 'Prod' :
       self._sourceDir = self._Sites[self._LocalSite]['treeBaseDir']+'/'+iProd+'/'+self._iniStep+'/'

     
     if not iStep == 'UEPS' : 

       self._targetDir = self._Sites[self._LocalSite]['treeBaseDir']+'/'+iProd+'/'
       if not self._iniStep == 'Prod' : self._targetDir += self._iniStep+'__'+iStep+'/'
       else                           : self._targetDir += iStep+'/'

       if self._Sites[self._LocalSite]['mkDir'] : os.system('mkdir -p '+ self._targetDir )

     # UEPS
     else:
       for iUEPS in Steps[iStep]['cpMap'] :
         os.system('mkdir -p '+ self._Sites[self._LocalSite]['treeBaseDir']+'/'+iProd+'/'+self._iniStep+'__'+iUEPS)
  
  


# --------------- Job Jandling

   def prepareJobs(self,iProd,iStep):

     bpostFix=''
     if not self._iniStep == 'Prod' : bpostFix='____'+self._iniStep

     # Make job directories
     jDir = jobDir+'/NanoGardening__'+iProd
     if not os.path.exists(jDir) : os.system('mkdir -p '+jDir)
     wDir = workDir+'/NanoGardening__'+iProd
     if not os.path.exists(wDir) : os.system('mkdir -p '+wDir)
   

     # prepare targetList
     targetList = []
     for iSample in self._targetDic :
       for iFile in self._targetDic[iSample] :
         iTarget = os.path.basename(self._targetDic[iSample][iFile]).replace(self._treeFilePrefix,'').replace('.root','')
         pidFile=jDir+'/NanoGardening__'+iProd+'__'+iStep+'__'+iTarget+bpostFix+'.jid'
         if os.path.isfile(pidFile) :
           print "pidFile", pidFile
           print '--> Job Running already : '+iTarget
         else: targetList.append(iTarget)

     # Dummy stepList for jobs
     stepList=[]
     stepList.append(iStep)

     # batchMode Preparation
     if self._jobMode == 'Batch':
       self._jobs = batchJobs('NanoGardening',iProd,stepList,targetList,'Targets,Steps',bpostFix)
       self._jobs.AddPy2Sh()
 

     for iSample in self._targetDic :
       for iFile in self._targetDic[iSample] :
         iTarget = os.path.basename(self._targetDic[iSample][iFile]).replace(self._treeFilePrefix,'').replace('.root','')
         if iTarget in targetList :
           pyFile=jDir+'/NanoGardening__'+iProd+'__'+iStep+'__'+iTarget+bpostFix+'.py'
           if os.path.isfile(pyFile) : os.system('rm '+pyFile)
           outTree=self._treeFilePrefix+iTarget+'__'+iStep+'.root'
           #GarbageCollector.append(outTree)
           self.mkPyCfg([iFile],iStep,pyFile,outTree,self._Productions[iProd]['isData'])
           command=''
           if self._jobMode == 'Interactive' : command='cd '+wDir+' ; python '+pyFile+' ; '
           if self._jobMode == 'Interactive' : os.system(command)
     
     if self._jobMode == 'Batch': self._jobs.Sub()


   def mkPyCfg(self,inputRootFiles,iStep,fPyName,haddFileName=None,isData=False):


     fPy = open(fPyName,'a') 
     
     # Common Header
     fPy.write('#!/usr/bin/env python \n')
     fPy.write('import os, sys \n')
     fPy.write('import ROOT \n')
     fPy.write('ROOT.PyConfig.IgnoreCommandLineOptions = True \n')
     fPy.write(' \n')
     fPy.write('from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor \n')
     fPy.write(' \n')

     # Import(s) of modules
     if self._Steps[iStep]['isChain'] :
       for iSubStep in  self._Steps[iStep]['subTargets'] :
         if 'import' in self._Steps[iSubStep] :
           fPy.write('from '+self._Steps[iSubStep]['import']+' import *\n')  
     else:
       if 'import' in self._Steps[iStep] :
         fPy.write('from '+self._Steps[iStep]['import']+' import *\n') 
     fPy.write(' \n')

     # Declaration(s) of in-line modules
     if self._Steps[iStep]['isChain'] :
       for iSubStep in  self._Steps[iStep]['subTargets'] :
         if 'declare' in self._Steps[iSubStep] :
           fPy.write(self._Steps[iSubStep]['declare']+'\n')
     else:
       if 'declare' in self._Steps[iStep] :
         fPy.write(self._Steps[iStep]['declare']+'\n') 
     fPy.write(' \n')

     # Files
     fPy.write('files=[')
     for iFile in inputRootFiles : fPy.write('"'+iFile+'",')
     fPy.write(']\n') 
     fPy.write(' \n')
     
     # Configure modules
     fPy.write('p = PostProcessor(  "."   ,          \n')
     fPy.write('                    files ,          \n')
     fPy.write('                    cut=None ,       \n')
     fPy.write('                    branchsel=None , \n')
     fPy.write('                    modules=[        \n')
     if self._Steps[iStep]['isChain'] :
       for iSubStep in  self._Steps[iStep]['subTargets'] :
         fPy.write('                          '+self._Steps[iSubStep]['module']+',\n')
     else:
       fPy.write('                          '+self._Steps[iStep]['module']+'\n') 
     fPy.write('                            ],      \n') 
     fPy.write('                    provenance=True, \n')
     fPy.write('                    fwkJobReport=True, \n')
     if not haddFileName == None :
       fPy.write('                    haddFileName="'+haddFileName+'", \n')
     fPy.write('                 ) \n')
     fPy.write(' \n')

     # Common footer
     fPy.write('p.run() \n')
     fPy.write(' \n')

     # Close file
     fPy.close()

#------------- Main 

   def process(self):

     for iProd in self._prodList:
       print '----------- Running on production: '+iProd
       self.readSampleFile(iProd) 

       for iStep in self._stepList:
         if    ( not self._Productions[iProd]['isData'] and self._Steps[iStep]['do4MC'] ) \
            or (     self._Productions[iProd]['isData'] and self._Steps[iStep]['do4Data'] ) :
           print '---------------- for Step : ',iStep
           self.mkFileDir(iProd,iStep)
           if not iStep == 'hadd' and not iStep == 'UEPS' :
             self.getTargetFiles(iProd,iStep)
             self.prepareJobs(iProd,iStep)
 
       self.Reset()
       

# ---- Testing ----

#Samples = { 
#              'GluGluHToWWTo2L2Nu_M125' : '/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/NANOAODSIM'
#          }
#iSample='GluGluHToWWTo2L2Nu_M125' 

#pp = PostProcMaker(PostProcSteps)
#Files=pp.getFilesFromDAS(Samples[iSample]) 
#FilesXrootd=[]
## Let's take 1 file only for test
#for iFile in Files:
##  FilesXrootd.append('root://cms-xrd-global.cern.ch//'+iFile)
#FilesXrootd.append('root://cms-xrd-global.cern.ch//'+Files[0])
#pp.mkPyCfg (FilesXrootd,'TestChain','aa.py','nanoLatino_'+iSample+'.root')