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
    self.selBJets = []
    self.pmet = TLorentzVector()

    # Create output histograms
    self.CreateTH1F("LepPt",   "", 24, 0, 120)
    self.CreateTH1F("LepEta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("Bjet0_pt",  "", 30, 0, 150)
    self.CreateTH1F("Bjet1_pt",  "", 30, 0, 150)
    self.CreateTH1F("jet0_pt",  "", 30, 0, 150)
    self.CreateTH1F("jet1_pt",  "", 30, 0, 150)
    self.CreateTH1F("Bjet0_eta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("Bjet1_eta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("jet0_eta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("jet1_eta",  "", 50, -2.5, 2.5)
    self.CreateTH1F("Bjets_InvMass",  "", 100, 0,500)
    self.CreateTH1F("jets_InvMass",  "", 100, 0,500)
    self.CreateTH1F("Bjets_DeltaPhi",  "", 100, -3.14/2,3.14/2)
    self.CreateTH1F("Bjets_DeltaPt",  "", 50,0,200)
    self.CreateTH1F("met_pt",  "", 50,0,200)
 
  def resetObjects(self):
    ''' Reset the list where the objects are stored '''
    self.selLeptons = []
    self.selJets = []
    self.selBJets = []
    self.pmet = TLorentzVector()

  def FillHistograms(self, lepton, bjets, jets, pmet):
    ''' Fill all the histograms. Take the inputs from lepton list, jet list, pmet '''
    if not len(lepton) >= 1: return # Just in case
    if not len(bjets) >= 2: return # Just in case
    if not len(jets) >= 2: return # Just in case
    self.weight = self.EventWeight * self.SFmuon * self.SFelec * self.PUSF

    # Re-calculate the observables
    lep_pt  = lepton.Pt()
    lep_eta = lepton.Eta()
    bjet0 = bjets[0] 
    bjet1 = bjets[1] 
    jet0 = jets[0] 
    jet1 = jets[1] 
    bjet0_pt = bjet0.Pt()
    bjet1_pt = bjet1.Pt()
    bjet0_eta = bjet0.Eta()
    bjet1_eta = bjet1.Eta()
    jet0_pt = jet0.Pt()
    jet1_pt = jet1.Pt()
    jet0_eta = jet0.Eta()
    jet1_eta = jet1.Eta()
    bjet_dphi  = DeltaPhi(bjet0, bjet1)
    mbb   = InvMass(bjet0, bjet1)
    mjj   = InvMass(jet0, jet1)
    bjet_dipt = DiPt(bjet0, bjet1)
    met_pt = pmet.Pt()
    
    ### Fill the histograms
    self.obj['Lep_pt'].Fill(lep_pt, self.weight)
    self.obj['Lep_eta'].Fill(lep_eta, self.weight)
    self.obj['Bjet0_pt'].Fill(bjet0_pt, self.weight)
    self.obj['Bjet0_eta'].Fill(bjet0_eta, self.weight)
    self.obj['Bjet1_pt'].Fill(bjet1_pt, self.weight)
    self.obj['Bjet1_eta'].Fill(bjet1_eta, self.weight)
    self.obj['jet0_pt'].Fill(jet0_pt, self.weight)
    self.obj['jet0_eta'].Fill(jet0_eta, self.weight)
    self.obj['jet1_pt'].Fill(jet1_pt, self.weight)
    self.obj['jet1_eta'].Fill(jet1_eta, self.weight)
    self.obj["Bjets_InvMass"].Fill(mbb, self.weight)
    self.obj["jets_InvMass"].Fill(mjj, self.weight)
    self.obj["Bjets_DeltaPhi"].Fill(bjet_dphi, self.weight)
    self.obj["Bjets_DeltaPt"].Fill(bjet_dipt, self.weight)
    self.obj["met_pt"].Fill(met_pt, self.weight)

  def insideLoop(self, t):
    self.resetObjects()

    ### Lepton selection
    ###########################################
    if not self.isData: nGenLep = t.nGenDressedLepton 
    
    ##### Jets
    for i in range (t.nJet):
      p = TLorentzVector()
      p.SetPtEtaPhiM(t.Jet_pt[i], t.Jet_eta[i], t.Jet_phi[i], t.Jet_mass[i])
      
      if p.Pt() < 30 or abs(p.Eta()) > 2.4:continue

      ## medium working point https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
      if t.Jet_btagDeepC > 0.4941 : 
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

    ##### MET 
    self.pmet = TLorentzVector()
    self.pmet.SetPtEtaPhiM(t.MET_pt, 0 , t.MET_phi, 0)

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
    ### We need one lepton, two bjets and two light jets 
    ### Each of them must be at least 10 GeV
    #def returnHighestPt(list_part,n):
    #    list_pt = []
    #    for part in list_part:
    #        list_pt.append(part[0].Pt())
    #    
    #    kept = sorted(range(len(a)), key=lambda i: a[i], reverse=True)[:n]         
    #    return list_part[kept]


    if not len(leps) >= 1:      return 
    #self.selLeptons = returnHighestPt(self.selLeptons,1)
    #self.selBJets= returnHighestPt(self.selBJets,2)
    #self.selJets= returnHighestPt(self.selJets,2)

        
    #if l0.charge*l1.charge > 0: return 
    #if l0.Pt() < 20:            return 
    #if InvMass(l0,l1) < 20:     return  

    ### Fill the histograms
    self.FillHistograms(self.selLeptons[:1], self.selBJets[:2], self.selJets[:2], self.pmet)
