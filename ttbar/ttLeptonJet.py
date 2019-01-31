import os,sys
basepath = os.path.abspath(__file__).rsplit('/ttCMSDAS/',1)[0]+'/ttCMSDAS/'
sys.path.append(basepath)

from framework.analysis import analysis
from framework.functions import DeltaPhi, DiPt, InvMass, lepton, jet
from ROOT.TMath import Sqrt as sqrt
from ROOT import *

################ Analysis
class ttLeptonJet(analysis):
  def init(self):
    # Load SF files
    if not self.isData:
      self.LoadHisto('MuonIsoSF', basepath+'./inputs/MuonISO.root', 'NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta') # pt, abseta
      self.LoadHisto('MuonIdSF',  basepath+'./inputs/MuonID.root',  'NUM_TightID_DEN_genTracks_pt_abseta') # pt, abseta
      self.LoadHisto('ElecSF',    basepath+'./inputs/ElecTightCBid94X.root',  'EGamma_SF2D') # eta, pt

    # Objects for the analysis
    self.selLeptons = []
    self.selJets = []
    self.pmet = TLorentzVector()

    # Create output histograms
    self.CreateTH1F("Lep0Pt",   "", 24, 0, 120)
    self.CreateTH1F("Lep0Eta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("BjetPt", "", 40, 0. , 400)
    self.CreateTH1F("BjetEta", "", 50, -2.5, 2.5)
    self.CreateTH1F("JetPt","", 40, 0. , 200)
    self.CreateTH1F("JetEta", "", 50, -2.5, 2.5)
    self.CreateTH1F("JetHT", "" , 40, 0. , 200)

  def resetObjects(self):
    ''' Reset the list where the objects are stored '''
    self.selLeptons = []
    self.selJets = []
    self.selBJets = []
    self.pmet = TLorentzVector()

  def FillHistograms(self, leptons, jets, bjets, pmet):
    ''' Fill all the histograms. Take the inputs from lepton list, jet list, pmet '''
    if not len(leptons) >= 1: return # Just in case
    self.weight = self.EventWeight * self.SFmuon * self.SFelec * self.PUSF

    # Re-calculate the observables
    lep0  = leptons[0]
    l0pt  = lep0.Pt()  
    l0eta = lep0.Eta()
    bjets0 = bjets[0]; bjets1 = bjets[1]
    bjets0pt = bjets0.Pt();bjets1pt = bjets1.Pt()
    bjets0eta = bjets0.Eta();bjets1eta = bjets1.Eta()

    #dphi  = DeltaPhi(lep0, lep1)
    #mll   = InvMass(lep0, lep1)
    #dipt  = DiPt(lep0, lep1)
    
    ### Fill the histograms
    self.obj['Lep0Pt'].Fill(l0pt, self.weight)
    #self.obj['Lep1Pt'].Fill(l1pt, self.weight)
    self.obj['Lep0Eta'].Fill(l0eta, self.weight)
    #self.obj['Lep1Eta'].Fill(l1eta, self.weight)
    #self.obj["InvMass"].Fill(mll, self.weight)
    #self.obj['DilepPt'].Fill(dipt, self.weight)
    #self.obj['DeltaPhi'].Fill(dphi/3.141592, self.weight)

  def insideLoop(self, t):
    self.resetObjects()

    ### Lepton selection
    ###########################################
    if not self.isData: nGenLep = t.nGenDressedLepton 
    
    ##### Jets
    for i in range (t.nJets):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Jets_pt[i], t.Jets_eta[i], t.Jets_phi[i], t.Jets_mass[i])
      
      if p.Pt() < 30 or abs(p.Eta()) > 2.4:continue

      ## medium working point https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
      if t.Jet_DeepC > 0.4941 : 
        self.selBJets.append(p)
      else:
        self.selJets.append(p)

    ##### Muons
    for i in range(t.nMuon):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Muon_pt[i], t.Muon_eta[i], t.Muon_phi[i], t.Muon_mass[i])
      charge = t.Muon_charge[i]

      # Tight ID, tight ISO, RelIso04 < 0.15, tight IP
      if not t.Muon_tightId[i]: continue # Tight ID
      if not t.Muon_pfRelIso04_all[i] < 0.15: continue
      dxy = abs(t.Muon_dxy[i])
      dz  = abs(t.Muon_dz[i] )
      if dxy > 0.05 or dz > 0.1: continue

      # pT > 12 GeV, |eta| < 2.4
      if p.Pt() < 12 or abs(p.Eta()) > 2.4: continue
      self.selLeptons.append(lepton(p, charge, 13))
       
    ##### Electrons
    for i in range(t.nElectron):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Electron_pt[i], t.Electron_eta[i], t.Electron_phi[i], t.Electron_mass[i])
      charge = t.Electron_charge[i]
      etaSC    = abs(p.Eta());
      convVeto = t.Electron_convVeto[i]

      # Tight cut-based Id, convVeto, RelIso03 tight, tight IP
      if not t.Electron_cutBased[i] >= 4: continue
      if not convVeto: continue
      relIso03 = t.Electron_pfRelIso03_all[i]
      if   etaSC <= 1.479 and relIso03 > 0.0361: continue
      elif etaSC >  1.479 and relIso03 > 0.094:  continue
      dxy = abs(t.Electron_dxy[i])
      dz  = abs(t.Electron_dz[i] )
      if dxy > 0.05 or dz > 0.1: continue

      # pT > 12 GeV, |eta| < 2.4
      if p.Pt() < 12 or abs(p.Eta()) > 2.4: continue
      self.selLeptons.append(lepton(p, charge, 11))

    leps = self.selLeptons
    pts  = [lep.Pt() for lep in leps]
    self.selLeptons = [lep for _,lep in sorted(zip(pts,leps))]

    ### Calculate the weights
    self.SFelec = 1; self.SFmuon = 1; self.SFelecErr = 0; self. SFmuonErr = 0
    if not self.isData:
      for lep in self.selLeptons:
        if lep.IsMuon():
          sf, err = self.GetSFandErr('MuonIsoSF, MuonIdSF', lep.Pt(), TMath.Abs(lep.Eta()))
          self.SFmuon*=sf
          self.SFmuonErr+=err*err
        else:
          sf, err = self.GetSFandErr('ElecSF', lep.Eta(), lep.Pt())
          self.SFelec*=sf
          self.SFelecErr+=err*err
      self.SFelecErr = sqrt(self.SFelecErr)
      self.SFmuonErr = sqrt(self.SFmuonErr)

    # PU SF --> PLEASE CHECK THAT THE WEIGHTS ARE IN THE TREES THAT YOU'RE USING!
    if not self.isData:
      self.PUSF   = t.puWeight
      self.PUUpSF = t.puWeightUp
      self.PUDoSF = t.puWeightDown
    else:
      self.PUSF   = 1; self.PUUpSF = 1; self.PUDoSF = 1

    ### Event selection
    ###########################################
    
    ### Dilepton pair: 2 leptons, opposite sign, mll > 20 GeV, leading lep pT > 20 GeV
    if not len(leps) >= 2:      return 
    l0 = leps[0]; l1 = leps[1]
    if l0.charge*l1.charge > 0: return 
    if l0.Pt() < 20:            return 
    if InvMass(l0,l1) < 20:     return  

    ### Fill the histograms
    self.FillHistograms(self.selLeptons, self.selJets, self.pmet)
