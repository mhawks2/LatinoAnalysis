import ROOT
import os
import re
import numpy
import math
import time
import copy
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import collectionMerger
from LatinoAnalysis.NanoGardener.framework.BranchMapping import mappedOutputTree, mappedEvent
 
import math
from itertools import combinations, permutations
import itertools
class BTagEventWeight():

    def _openRootFile(self,path, option=''):
        f =  ROOT.TFile.Open(path,option)
        if not f.__nonzero__() or not f.IsOpen():
            raise NameError('File '+path+' not open')
        return f

    def _getRootObj(self,d,name):
        o = d.Get(name)
        if not o.__nonzero__():
            print 'Object '+name+' doesn\'t exist in '+d.GetName(), ' BE CAREFUL!'
        return o

    def __init__ (self,cmssw):
        self.cmssw = cmssw
        cmssw_base = os.getenv('CMSSW_BASE')
  
        #deepjet medium WPs
        WPdic = {'Full2018v9': 0.2738, 'Full2017v9': 0.3040, 'Full2016v9noHIPM': 0.2489, 'Full2016v9HIPM': 0.2598} 
        self.WP = WPdic[self.cmssw]

        # Root File
        self.rootfile = self._openRootFile(cmssw_base+'/src/LatinoAnalysis/NanoGardener/python/data/AZH_btagEff/'+self.cmssw+'/ttZ_btagEff_deepjet_M.root')

        # Root Histos
        self.h_b_Eff    = self._getRootObj(self.rootfile, 'b_btagEff_deepjet_M_ptbin')
        self.h_c_Eff    = self._getRootObj(self.rootfile, 'c_btagEff_deepjet_M_ptbin')
        self.h_udsg_Eff = self._getRootObj(self.rootfile, 'udsg_btagEff_deepjet_M_ptbin')

    def _getEff(self, pt, flavour):

        if   (flavour == 5): hist = self.h_b_Eff
        elif (flavour == 4): hist = self.h_c_Eff
        elif (flavour == 0): hist = self.h_udsg_Eff
        
        nbins = hist.GetNbinsX()
        ptmax = hist.GetXaxis().GetBinCenter(nbins)

        eff_value = hist.GetBinContent(hist.FindBin(min(pt, ptmax)))
        eff_error = hist.GetBinError  (hist.FindBin(min(pt, ptmax)))
        
        return eff_value


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # _getProductMC
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _getProductMC(self, cleanjets_col, jets_col):

        CleanJets = cleanjets_col
        Jets = jets_col
        
        p_tagged     = numpy.ones(1, dtype=numpy.float32)
        p_not_tagged = numpy.ones(1, dtype=numpy.float32)
        
        for cj in CleanJets:
            if (cj.pt > 30.0 and abs(cj.eta) < 2.5):
                pt = cj.pt
                flav = Jets[cj.jetIdx].hadronFlavour
                disc = Jets[cj.jetIdx].btagDeepFlavB

                isTagged = True if (disc > self.WP) else False
                if (isTagged):
                    p_tagged     *=      self._getEff(pt, flav)
                else:
                    p_not_tagged *= (1 - self._getEff(pt, flav))

        btagweight = p_tagged * p_not_tagged

        return btagweight


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # _getProductData
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def _getProductDATA(self, cleanjets_col, jets_col):

        CleanJets = cleanjets_col
        Jets = jets_col
        
        p_tagged     = numpy.ones(1, dtype=numpy.float32)
        p_not_tagged = numpy.ones(1, dtype=numpy.float32)
        
        for cj in CleanJets:
            if (cj.pt > 30.0 and abs(cj.eta) < 2.5):
                pt   = cj.pt
                flav = Jets[cj.jetIdx].hadronFlavour
                disc = Jets[cj.jetIdx].btagDeepFlavB
                sf   = Jets[cj.jetIdx].btagSF_deepjet_M

                isTagged = True if (disc > self.WP) else False
                if (isTagged):
                    p_tagged     *=      (sf * self._getEff(pt, flav))
                else:
                    p_not_tagged *= (1 - (sf * self._getEff(pt, flav)))

        btagweight = p_tagged * p_not_tagged

        return btagweight


class BTagEventWeightProducer_WP(Module):
    '''
    Produce branches with btag event weights calculated with method 1(a)
    https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a_Event_reweighting_using_scale

    only produces weights for events with >=4 jets with pt > 30 and abs(eta) < 2.5
    '''

    def __init__(self, cmssw):
        self.cmssw = cmssw     
        cmssw_base = os.getenv('CMSSW_BASE')
        print " cmssw = ", self.cmssw

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        print 'Creating branch : AZH_btagWeight4j' 
        self.out.branch('AZH_btagWeight4j', 'F')            

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jet_col      = Collection(event, 'Jet')
        CleanJet_col = Collection(event, 'CleanJet') 

        njet = 0
        for j in CleanJet_col:
            if (j.pt > 30.0 and abs(j.eta) < 2.5): njet += 1
 
        weight = 1.0
        if (njet >= 4):
            p_MC   = BTagEventWeight(self.cmssw)._getProductMC(CleanJet_col, jet_col)
            p_DATA = BTagEventWeight(self.cmssw)._getProductDATA(CleanJet_col, jet_col)            
            weight *= (p_DATA / p_MC)

        self.out.fillBranch('AZH_btagWeight4j', weight)

        return True
