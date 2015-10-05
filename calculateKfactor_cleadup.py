#!/bin/env python
import ROOT
import array
import glob
from math import ceil,sqrt,fabs
import sys
def GetKeyNames( self, dir = "" ):
        self.cd(dir)
        return [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
ROOT.TFile.GetKeyNames = GetKeyNames

def uniqueName(h):
    h.SetName(h.GetName()+str(id(h)))

def getRatio(h1,h2,name=None):
    if name==None:
        name=h1.GetName()
    ratio=h1.Clone(name)
    uniqueName(ratio)
    for ibin in range(h1.GetNbinsX()):
        sa=h1.GetBinError(ibin)
        a=h1.GetBinContent(ibin)
        sb=h2.GetBinError(ibin)
        b=h2.GetBinContent(ibin)
        if b!=0 and a!=0:
            ratio.SetBinContent(ibin,a/b)
            ratio.SetBinError(ibin, sqrt(fabs((sa/a)**2 + (sb/b)**2 - 2.*sa**2/(a*b)) )   )
        else:
            ratio.SetBinContent(ibin,0)
            ratio.SetBinError(ibin,0)
    return ratio


def main():
    inputFile=ROOT.TFile("StoredN_NLO_hists.root","read")
    binning10=range(0,8001,50)

    hists={}
    for h in inputFile.GetKeyNames():
        hists[h]=inputFile.Get(h).Rebin(len(binning10)-1,inputFile.Get(h).GetName()+str(id(inputFile.Get(h))),array.array("d",binning10))
        print h

    mgFile=ROOT.TFile("Wmass.root","READ")
    LO_MGpPythia=mgFile.Get(mgFile.GetListOfKeys()[0].GetName())

    #this is done for 1 fb and ele + mu correct for one only
    #mgHist.Scale(1./(3.*1000.))
    LO_MGpPythia.Scale(1./3.)
    LO_MGpPythiaReb_h=LO_MGpPythia.Rebin(len(binning10)-1,LO_MGpPythia.GetName()+str(id(LO_MGpPythia)),array.array("d",binning10))
    LO_MGpPythiaReb_h.SetLineColor(ROOT.kOrange)
    for lepton in ["ele","tau","mu"]:
        c1=ROOT.TCanvas()
        ew_ratio=getRatio(hists["QCD_LO_EW_NLO_PhotonPDF_%s_mcsanc"%(lepton)],hists["QCD_LO_EW_LO_PhotonPDF_mcsanc"],name="EW-kfac")
        #qcd_ratio=getRatio(hists["QCD_NNLO_FEWZ"],hists["LO_MGpPythia"],name="QCD-kfac")

        ew_m_lo=hists["QCD_LO_EW_LO_PhotonPDF_mcsanc"].Clone("ew_m_lo")
        uniqueName(ew_m_lo)
        ew_m_lo.Add(hists["QCD_LO_EW_NLO_PhotonPDF_%s_mcsanc"%(lepton)],hists["QCD_LO_EW_LO_PhotonPDF_mcsanc"],1,-1)


        k_fakp=hists["QCD_LO_EW_LO_PhotonPDF_mcsanc"].Clone("k_fac_p")
        uniqueName(k_fakp)
        k_fakp.Add(ew_m_lo,hists["QCD_NNLO_FEWZ"],1,1)


        k_fakm=hists["QCD_LO_EW_LO_PhotonPDF_mcsanc"].Clone("k_fac_m")
        uniqueName(k_fakm)
        k_fakm.Multiply(hists["QCD_NNLO_FEWZ"],ew_ratio,1,1)

        k_fakm.SetLineColor(ROOT.kRed)

        k_fakm.Divide(k_fakm,LO_MGpPythiaReb_h,1,1)
        k_fakp.Divide(k_fakp,LO_MGpPythiaReb_h,1,1)

        k_fakm.Draw()
        #cosmetics
        c1.SetBottomMargin(0.12)
        k_fakm.SetTitle("")
        k_fakm.GetXaxis().SetTitle("M_{W} [GeV]")
        k_fakm.GetYaxis().SetTitleOffset(0.8)
        k_fakm.GetYaxis().SetTitle("(N)NLO/LO")
        k_fakp.Draw("same")

        k_fak_mean=k_fakp.Clone("mean_k_fak")
        uniqueName(k_fak_mean)
        k_fak_mean.SetLineColor(ROOT.kBlue)
        k_fak_mean.Reset()
        k_fak_mean.Add(k_fakp,k_fakm,1,1)
        k_fak_mean.Scale(1./2.)
        k_fak_mean.Draw("same")
        raw_input()
        c1.SaveAs("kfactor_%s.pdf"%(lepton))

        outFile=ROOT.TFile("k_faktors_%s.root"%(lepton),"recreate")
        outFile.cd()
        k_fakp.Write("k_fac_p")
        k_fakm.Write("k_fac_m")
        k_fak_mean.Write("k_fac_mean")
        outFile.Close()





main()
