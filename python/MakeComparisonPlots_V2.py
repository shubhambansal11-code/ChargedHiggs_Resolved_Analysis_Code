# -*- coding: utf-8 -*-
#python
import sys
import glob
import math
import re
from ROOT import *
from array import *
#from ROOT import TCanvas, TFile, TLine, TProfile, TNtuple, TH1F, TH2F
#import ROOT
#import array

gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gROOT.LoadMacro("../style/AtlasStyle.C")
gROOT.LoadMacro("../style/AtlasUtils.C")
SetAtlasStyle()

y_axis_label="Event fraction"

c_blue   = TColor.GetColor("#3366ff")
c_red    = TColor.GetColor("#ee3311")
c_orange = TColor.GetColor("#ff9900")

def NormalizeHisto(histo):
     n_events=histo.Integral(-1,histo.GetNbinsX()+1)
     #n_events=histo.Integral()
     if n_events == 0:
         return
     print n_events, histo.Integral(histo.GetNbinsX(),histo.GetNbinsX()+1)
     histo.Scale(1./n_events)
     histo.SetLineWidth(2)
     histo.SetStats(0)
     histo.SetFillStyle(3001)
     histo.SetMarkerColor(histo.GetLineColor())
     histo.SetMarkerSize(0.0)
     histo.GetXaxis().SetTitleOffset(1.2)
     histo.GetYaxis().SetTitleOffset(1.52)
     histo.GetXaxis().SetLabelSize(0.05)
     histo.GetYaxis().SetLabelSize(0.05)
     histo.GetXaxis().SetTitleSize(0.05)
     histo.GetYaxis().SetTitleSize(0.05)
     histo.GetYaxis().SetNdivisions(504)
     histo.GetXaxis().SetNdivisions(504)

c1 = TCanvas("ShapePlots","",500,500)

#for HistoName in ["nBTags"]:
for HistoName in ["nBTags","nJets","HT","HT_bjets","DeltaPhi_HW","mVH","mH","pTWplus","pTH","maxMVAResponse", "mass_resolution"]: # for resolved, histoname set
#for HistoName in ["nBTags","nJets","HT","HT_bjets","DeltaPhi_HW","mVH","mH","pTWplus","pTH","mass_resolution"]:#for boosted, histoname set
#for HistoName in ["mVH"]:
#for HistoName in ["mVH"]:     
#for HistoName in ["maxMVAResponse"]:
## for bTagStrategy in ["Incl","FourPlusTags","ThreeTags","TwoTags"]:
 for Region in ["Resolved_SR"]:
  for btagStrategy in ["Inclusive"]:
   #for btagStrategy in ["Inclusive","FourPlusTags","ThreeTags","TwoTags"]:
     
      ReBin = False
      YAxisScale = 1.4

      if "nBTags" in HistoName:
          Xaxis_label="b-tag multiplicity"
      if "nJets" in HistoName:
          Xaxis_label="Jet Multiplicity"
      if "DeltaPhi_HW" in HistoName:
          Xaxis_label="DeltaPhi_HW"
      if "pTH" in HistoName:
          Xaxis_label="Transverse Momentum of Higgs [GeV]"
      if "pTWplus" in HistoName:
          Xaxis_label="Transverse Momentum of W Boson [GeV]"
      if "mVH" in HistoName:
          Xaxis_label="Mass of Charged Higgs [GeV]"
      if "mH" in HistoName:
          Xaxis_label="Mass of Higgs [GeV]"
      if "pTWplus" in HistoName:
          Xaxis_label="Transverse Momentum Of Wplus [GeV]"
      if "pTH" in HistoName:
          Xaxis_label="Transverse Momentum Of Higgs [GeV]"
      if "mass_resolution" in HistoName:
          Xaxis_label="Mass Resolution"
      if "HT" in HistoName:
          Xaxis_label="H_{T} (Scalar Transverse Momentum Sum of jets) [GeV]"
      if "HT_bjets" in HistoName:
          Xaxis_label="H_{T_{b-jet}} (Scalar Transverse Momentum Sum of b-jets) [GeV]"
      if "maxMVAResponse" in HistoName:
          Xaxis_label="BDT Score (Signal Reconstruction)"    
      #Xaxis_label=""

      #file1       = TFile.Open("../PlotFiles/ForOptimisation/sig_Hplus_Wh_m400-0_70p.root","READ")
      file1       = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/sig_Hplus_Wh_m400-0_70p.root","READ")
      dir1        = file1.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir1        = file1.GetDirectory("Nominal").GetDirectory(HistoName)
      h_sig_Hplus_m400 = dir1.Get("sig_Hplus_Wh_m400-0_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_sig_Hplus_m400.SetLineColor(kRed)
      h_sig_Hplus_m400.SetLineStyle(7)
      if ReBin == True:
          h_sig_Hplus_m400.Rebin(2)

      #file6	  = TFile.Open("../PlotFiles/ForOptimisation/sig_Hplus_Wh_m800-0_70p.root","READ")
      file6	  = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/sig_Hplus_Wh_m800-0_70p.root","READ")
      dir6        = file6.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir6        = file6.GetDirectory("Nominal").GetDirectory(HistoName)
      h_sig_Hplus_m800 = dir6.Get("sig_Hplus_Wh_m800-0_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_sig_Hplus_m800.SetLineColor(kBlack)
      h_sig_Hplus_m800.SetLineStyle(7)
      if ReBin == True:
          h_sig_Hplus_m800.Rebin(2)

      #file7	  = TFile.Open("../PlotFiles/ForOptimisation/sig_Hplus_Wh_m1600-0_70p.root","READ")
      file7	  = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/sig_Hplus_Wh_m1600-0_70p.root","READ")
      dir7        = file7.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir7        = file7.GetDirectory("Nominal").GetDirectory(HistoName)
      h_sig_Hplus_m1600 = dir7.Get("sig_Hplus_Wh_m1600-0_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_sig_Hplus_m1600.SetLineColor(kBlue)
      h_sig_Hplus_m1600.SetLineStyle(7)
      if ReBin == True:
          h_sig_Hplus_m1600.Rebin(2)
        
     # #file2   = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/ttbar_70p.root","READ")
      file2   = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/ttbar_70p.root","READ")
      dir2    = file2.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir2    = file2.GetDirectory("Nominal").GetDirectory(HistoName)
      h_ttbar_background = dir2.Get("ttbar_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_ttbar_background.SetLineColor(kGreen)
      h_ttbar_background.SetLineStyle(3)
      if ReBin == True:
          h_ttbar_background.Rebin(2)

      #file3   = TFile.Open("../PlotFiles/ForOptimisation/Wjets_70p.root","READ")
      file3   = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/Wjets_70p.root","READ")
      dir3    = file3.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir3    = file3.GetDirectory("Nominal").GetDirectory(HistoName)
      h_W_background = dir3.Get("Wjets_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_W_background.SetLineColor(kMagenta)
      h_W_background.SetLineStyle(3)
      if ReBin == True:
          h_W_background.Rebin(2)
      
      #file4   = TFile.Open("../PlotFiles/ForOptimisation/diboson_70p.root","READ")
      file4   = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/diboson_70p.root","READ")
      dir4    = file4.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir4    = file4.GetDirectory("Nominal").GetDirectory(HistoName)
      h_diboson_background = dir4.Get("diboson_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_diboson_background.SetLineColor(kMagenta)
      h_diboson_background.SetLineStyle(3)
      if ReBin == True:
          h_diboson_background.Rebin(2)

      #file5   = TFile.Open("../PlotFiles/ForOptimisation/singleTop_70p.root","READ")
      file5   = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/singleTop_70p.root","READ")
      dir5    = file5.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir5    = file5.GetDirectory("Nominal").GetDirectory(HistoName)
      h_singleTop_background = dir5.Get("singleTop_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_singleTop_background.SetLineColor(kMagenta)
      h_singleTop_background.SetLineStyle(3)
      if ReBin == True:
          h_singleTop_background.Rebin(2)

      #file8   = TFile.Open("../PlotFiles/ForOptimisation/Zjets_70p.root","READ")
      file8   = TFile.Open("/cephfs/user/s6subans/ChargedHiggs_Mock/output_SR_NewBDT_imp2sig/Zjets_70p.root","READ")
      dir8    = file8.GetDirectory("Nominal").GetDirectory(HistoName)
      #dir8    = file8.GetDirectory("Nominal").GetDirectory(HistoName)
      h_Z_background = dir8.Get("Zjets_"+HistoName+"_"+Region+"_"+btagStrategy+"_Nominal")
      h_Z_background.SetLineColor(kMagenta)
      h_Z_background.SetLineStyle(3)
      if ReBin == True:
          h_Z_background.Rebin(2)    

      h_other_background = h_singleTop_background + h_diboson_background + h_W_background + h_Z_background

      nbins=20
      ymax=0
      NormalizeHisto(h_other_background)
      if ymax<h_other_background.GetMaximum():
          ymax=h_other_background.GetMaximum()
      NormalizeHisto(h_ttbar_background)
      if ymax<h_ttbar_background.GetMaximum():
          ymax=h_ttbar_background.GetMaximum()
      NormalizeHisto(h_sig_Hplus_m400)
      if ymax<h_sig_Hplus_m400.GetMaximum():
          ymax=h_sig_Hplus_m400.GetMaximum()
      NormalizeHisto(h_sig_Hplus_m800)
      if ymax<h_sig_Hplus_m800.GetMaximum():
          ymax=h_sig_Hplus_m800.GetMaximum()
      NormalizeHisto(h_sig_Hplus_m1600)
      if ymax<h_sig_Hplus_m1600.GetMaximum():
          ymax=h_sig_Hplus_m1600.GetMaximum()


      h_other_background.Draw("HIST")
      h_ttbar_background.Draw("HISTSAME")
      h_sig_Hplus_m400.Draw("HISTSAME")
      h_sig_Hplus_m800.Draw("HISTSAME")
      h_sig_Hplus_m1600.Draw("HISTSAME")

      if HistoName in "maxMVAResponse":
         leg = TLegend(0.2,0.65,0.725,0.855)
      else:
         leg = TLegend(0.45,0.65,0.925,0.855)
      ATLAS_LABEL(0.20,0.875," Simulation Internal",1,0.19);
      leg.SetShadowColor(kWhite)
      leg.SetFillColor(kWhite)
      leg.SetLineColor(kWhite)
      
      leg.AddEntry(h_sig_Hplus_m400,    "H^{+}#rightarrow W^{+}h(m_{H^{+}} = 400GeV)","L")
      leg.AddEntry(h_sig_Hplus_m800,    "H^{+}#rightarrow W^{+}h(m_{H^{+}} = 800GeV)","L")
      leg.AddEntry(h_sig_Hplus_m1600,   "H^{+}#rightarrow W^{+}h(m_{H^{+}} = 1600GeV)","L")
      leg.AddEntry(h_ttbar_background,  "t#bar{t}","L")
      leg.AddEntry(h_other_background,  "other backgrounds","L")
      leg.SetTextSize(0.0250)
      leg.Draw()
      h_other_background.GetXaxis().SetTitle(Xaxis_label)
      h_other_background.GetYaxis().SetRangeUser(0.001,round(ymax*1.25,3)+0.001)
      
      c1.RedrawAxis()
      c1.Update()
      c1.RedrawAxis()
      c1.SaveAs("../Plots/ShapePlot_%s_qqbb_SR_NewBDT.pdf" % (HistoName+"_"+btagStrategy))
