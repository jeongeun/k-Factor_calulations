#!/bin/env python
import ROOT
import array
#from numpy import log
import glob
from math import ceil,sqrt,fabs
import sys
import helper
import cPickle

helper.setTDRStyle(1)

ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def correctBinning(h,origianlW):
    for ibin in range(h.GetNbinsX()):
        h.SetBinContent(ibin,h.GetBinContent(ibin)/h.GetBinWidth(ibin)*origianlW)
        h.SetBinError(ibin,h.GetBinError(ibin)/h.GetBinWidth(ibin)*origianlW)


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



def frange(start, end=None, inc=None):
    """A range function, that does accept float increments..."""
    import math

    if end == None:
        end = start + 0.0
        start = 0.0
    else: start += 0.0 # force it to be a float

    if inc == None:
        inc = 1.0
    count = int(math.ceil((end - start) / inc))

    L = [None,] * count

    L[0] = start
    for i in xrange(1,count):
        L[i] = L[i-1] + inc
    return L



def uniqueName(h):
    h.SetName(h.GetName()+str(id(h)))


def getmcsanc(file,hist):
    mcsanc_raw=[]
    binning=[]
    mcsancFile=open(file,"r")
    header=""
    for line in mcsancFile:
        if "\t" not in line:
            line=line.strip()
            line=line.strip("----")
            line=line.strip(" ")
            header=line
        if header==hist:
            tmp=line.split()
            if len(tmp)>1:
                #25.00    18.145565       0.69884075E-01
                if "inf" in tmp[2] or "nan" in tmp[2]:
                    mcsanc_raw.append([0.,0.,float(tmp[1])])
                    continue
                mcsanc_raw.append([float(tmp[2]),float(tmp[3]),float(tmp[1])])
                binning.append(float(tmp[0]))
                #last=tmp[1]
    if len(mcsanc_raw)==0:
        return [None,None]
    binning.append(mcsanc_raw[-1][2])

    mcsanc_m_Hist=ROOT.TH1D(file,file,8000,0,8000)
    for i,h in enumerate(mcsanc_raw):
        ibin=mcsanc_m_Hist.FindBin(h[2])
        mcsanc_m_Hist.SetBinContent(ibin,h[0])
        mcsanc_m_Hist.SetBinError(ibin,h[1])
    return mcsanc_m_Hist,binning

def getfewz(file,hist):
    #NNLOQCD
    fewzhist_raw=[]
    mcsancFile=open(file,"r")
    if file.split("_")[1]=="LO":
        file=file.replace("_LO_","_")
    pm=int(file.split("_")[3])


    header=""
    for line in mcsancFile:
        if "  ----" in line:
            line=line.strip()
            line=line.strip("----")
            line=line.strip(" ")
            header=line
        if header==hist:
            tmp=line.split()
            if len(tmp)>2:
                fewzhist_raw.append([float(tmp[0]),float(tmp[1]),float(tmp[2])])

    #ok lets make a 10GeV normalised binning
    width=10
    for ibin in range(len(fewzhist_raw)):
        fewzhist_raw[ibin][1]=fewzhist_raw[ibin][1]/width
        fewzhist_raw[ibin][2]=fewzhist_raw[ibin][2]/width
    return fewzhist_raw,pm



def makeBinningV2(plus,minus):
    combined_raw=[]



    combined_binning=[i for i in range(0,8001,1)]

    fewzHist=ROOT.TH1D("fewz_combined","fewz_combined",8000,0,8000)
    fewzHist=fewzHist.Rebin(len(combined_binning)-1,fewzHist.GetName()+str(id(fewzHist)),array.array("d",combined_binning))

    plusComb=[]
    for i in plus:
        plusComb+=i
    minusComb=[]
    for i in minus:
        minusComb+=i

    minusComb.sort()
    plusComb.sort()
    for p,m in zip(plusComb,minusComb):
        ibin=fewzHist.FindBin(p[0])
        ibin2=fewzHist.FindBin(m[0])
        if m[0]!=p[0]:
            print "That is wrong here ",m,p
        if (fewzHist.GetBinContent(ibin)!=0):
            fewzHist.SetBinContent(ibin,p[1]+m[1]+fewzHist.GetBinContent(ibin))
            fewzHist.SetBinError(ibin,sqrt(p[2]**2 + m[2]**2+fewzHist.GetBinError(ibin)**2.))
        else:
            fewzHist.SetBinContent(ibin,p[1]+m[1])
            fewzHist.SetBinError(ibin,sqrt(p[2]**2 + m[2]**2))


    binning=range(0,8001,1)
    binning10=range(0,8001,10)
    fewzHist=fewzHist.Rebin(len(binning10)-1,fewzHist.GetName()+str(id(fewzHist)),array.array("d",binning10))
    fewzHist.Scale(10)
    return fewzHist








#create different binnings!
binning=range(0,8001,1)
binning10=range(0,8001,10)

logBinning=[1]
for i in range(8000):
    logBinning.append(logBinning[-1]*1.5)
    if logBinning[-1]> 8000:
        break
logBinning=range(0,8001,10)

NNLO_output=ROOT.TFile("StoredN_NLO_hists.root","recreate")
NNLO_output.cd()
#QCD NNLO:

#fewz

fewzHistsp=[]
fewzHistsm=[]
binningp=[]
binningm=[]
for f in glob.glob("fewz/NNLO.fewz*"):
    h,pm=getfewz(f,"Q_ll Invaria")
    if pm==1:
        fewzHistsm.append(h)
    if pm==3:
        fewzHistsp.append(h)
fewzHist=makeBinningV2(fewzHistsp,fewzHistsm)
fewzHist.Write("QCD_NNLO_FEWZ")


#QCD NLO:
fewzNLOHistsp=[]
fewzNLOHistsm=[]
binningNLOp=[]
binningNLOm=[]
for f in glob.glob("fewzNLO/NLO.fewz*"):
    h,pm=getfewz(f,"Q_ll Invaria")
    if pm==1:
        fewzNLOHistsm.append(h)
    if pm==3:
        fewzNLOHistsp.append(h)
fewzHistNLO=makeBinningV2(fewzNLOHistsp,fewzNLOHistsm)
fewzHistNLO.Write("QCD_NLO_FEWZ")


#LO mcsanc with LO PDF
mcsancLOLOHists=[]
binning=[]
for f in glob.glob("mcsancLOLO/*103*.txt"):
    h,binning=getmcsanc(f,"m34-calo-fix")
    mcsancLOLOHists.append(h)
mcsancLOLOEWhist=mcsancLOLOHists[0]
for h in mcsancLOLOHists[1:]:
    mcsancLOLOEWhist.Add(mcsancLOLOEWhist,h,1,1)
mcsancLOLOEWhist=mcsancLOLOEWhist.Rebin(len(binning10)-1,mcsancLOLOEWhist.GetName()+str(id(mcsancLOLOEWhist)),array.array("d",binning10))
correctBinning(mcsancLOLOEWhist,mcsancLOLOEWhist.GetBinWidth(20))
mcsancLOLOEWhist.Write("QCD_LO_EW_LO_mcsanc")
mcsancLOLOEWhist.SetLineColor(ROOT.kGreen)
mcsancLOLOEWhist.Draw("same")





#LOMG  this is the LO change here
mgFile=ROOT.TFile("Wmass.root","READ")
mgHist=None
#for key in mgFile.GetListOfKeys():
    #if "06295cef49794cd580f0563052034260" in key.GetName():
mgHist=mgFile.Get(mgFile.GetListOfKeys()[0].GetName())

#this is done for 1 fb and ele + mu correct for one only
#mgHist.Scale(1./(3.*1000.))
mgHist.Scale(1./(3.))
mgHistReb_h=mgHist.Rebin(len(binning10)-1,mgHist.GetName()+str(id(mgHist)),array.array("d",binning10))
mgHistReb_h.SetLineColor(ROOT.kOrange)
NNLO_output.cd()
mgHistReb_h.Write("LO_MGpPythia")





#mcsanc with NNLO pdf with phonton in the initial state this we need
mcsancPhotonLOHists=[]
binning=[]
for f in glob.glob("mcsancEWphotonLO/*103*.txt"):
    h,binning=getmcsanc(f,"m34-calo-fix")
    mcsancPhotonLOHists.append(h)
mcsancPhotonLOEWhist=mcsancPhotonLOHists[0]
for h in mcsancPhotonLOHists[1:]:
    mcsancPhotonLOEWhist.Add(mcsancPhotonLOEWhist,h,1,1)
mcsancPhotonLOEWhist=mcsancPhotonLOEWhist.Rebin(len(binning10)-1,mcsancPhotonLOEWhist.GetName()+str(id(mcsancPhotonLOEWhist)),array.array("d",binning10))
correctBinning(mcsancPhotonLOEWhist,mcsancPhotonLOEWhist.GetBinWidth(20))
mcsancPhotonLOEWhist.Write("QCD_LO_EW_LO_PhotonPDF_mcsanc")



#EW NLO

#tau
mcsancPhotonHists=[]
binning=[]
for f in glob.glob("mcsancEWphoton/*103*-output.txt"):
    h,binning=getmcsanc(f,"m34-calo-fix")
    mcsancPhotonHists.append(h)
mcsancPhotonEWhist=mcsancPhotonHists[0]
for h in mcsancPhotonHists[1:]:
    mcsancPhotonEWhist.Add(mcsancPhotonEWhist,h,1,1)
mcsancPhotonEWhist=mcsancPhotonEWhist.Rebin(len(binning10)-1,mcsancPhotonEWhist.GetName()+str(id(mcsancPhotonEWhist)),array.array("d",binning10))
correctBinning(mcsancPhotonEWhist,mcsancPhotonEWhist.GetBinWidth(20))
mcsancPhotonEWhist.Write("QCD_LO_EW_NLO_PhotonPDF_tau_mcsanc")

#ele
mcsancPhotonEleHists=[]
binning=[]
for f in glob.glob("mcsancEWphoton_ele/*102*-output.txt"):
    h,binning=getmcsanc(f,"m34-calo-fix")
    mcsancPhotonEleHists.append(h)
mcsancPhotonEleEWhist=mcsancPhotonEleHists[0]
for h in mcsancPhotonEleHists[1:]:
    mcsancPhotonEleEWhist.Add(mcsancPhotonEleEWhist,h,1,1)
mcsancPhotonEleEWhist=mcsancPhotonEleEWhist.Rebin(len(binning10)-1,mcsancPhotonEleEWhist.GetName()+str(id(mcsancPhotonEleEWhist)),array.array("d",binning10))
correctBinning(mcsancPhotonEleEWhist,mcsancPhotonEleEWhist.GetBinWidth(20))
mcsancPhotonEleEWhist.Write("QCD_LO_EW_NLO_PhotonPDF_ele_mcsanc")

#mu
mcsancPhotonMuHists=[]
binning=[]
for f in glob.glob("mcsancEWphoton_mu/*101*-output.txt"):
    h,binning=getmcsanc(f,"m34-bare-fix")
    mcsancPhotonMuHists.append(h)
mcsancPhotonMuEWhist=mcsancPhotonMuHists[0]
for h in mcsancPhotonMuHists[1:]:
    mcsancPhotonMuEWhist.Add(mcsancPhotonMuEWhist,h,1,1)
mcsancPhotonMuEWhist=mcsancPhotonMuEWhist.Rebin(len(binning10)-1,mcsancPhotonMuEWhist.GetName()+str(id(mcsancPhotonMuEWhist)),array.array("d",binning10))
correctBinning(mcsancPhotonMuEWhist,mcsancPhotonMuEWhist.GetBinWidth(20))
mcsancPhotonMuEWhist.Write("QCD_LO_EW_NLO_PhotonPDF_mu_mcsanc")

###Draw controll plot
leg=ROOT.TLegend(0.67,0.67,0.92,0.92)

fewzHist.UseCurrentStyle()
mcsancPhotonEWhist.UseCurrentStyle()
mcsancPhotonLOEWhist.UseCurrentStyle()
fewzHist.SetLineColor(ROOT.kRed)
mcsancPhotonEWhist.SetLineColor(ROOT.kAzure)


mcsancPhotonEWhist.Draw("h")
mcsancLOLOEWhist.Draw("same h")
fewzHist.Draw("same h")



leg.AddEntry(mcsancLOLOEWhist,"LO","l")
leg.AddEntry(mcsancPhotonEWhist,"NLO ew","l")
leg.AddEntry(fewzHist,"NNLO qcd","l")

mcsancPhotonEWhist.GetXaxis().SetRangeUser(30,8000)
mcsancPhotonEWhist.GetXaxis().SetTitle("M_{inv} [GeV]")
mcsancPhotonEWhist.GetYaxis().SetTitle("d#sigma/dM_{#tau#nu} (pb/10 GeV)")


leg.Draw()



raw_input("Does it look ok?")


####### Which lepton do we want?
mcsancPhotonEWhist=mcsancPhotonMuEWhist


fewzHist=fewzHist.Rebin(len(logBinning)-1,fewzHist.GetName()+str(id(fewzHist)),array.array("d",logBinning))
mgHistReb_h=mgHistReb_h.Rebin(len(logBinning)-1,mgHistReb_h.GetName()+str(id(mgHistReb_h)),array.array("d",logBinning))
mcsancPhotonEWhist=mcsancPhotonEWhist.Rebin(len(logBinning)-1,mcsancPhotonEWhist.GetName()+str(id(mcsancPhotonEWhist)),array.array("d",logBinning))
mcsancPhotonLOEWhist=mcsancPhotonLOEWhist.Rebin(len(logBinning)-1,mcsancPhotonLOEWhist.GetName()+str(id(mcsancPhotonLOEWhist)),array.array("d",logBinning))
loHist=mcsancLOLOEWhist

c1=ROOT.TCanvas("c","c",800,800)


ew_ratio=getRatio(mcsancPhotonEWhist,mcsancPhotonLOEWhist,name="EW-kfac")

qcd_ratio=getRatio(fewzHist,loHist,name="QCD-kfac")



qcd_ratio.SetLineColor(ROOT.kRed)
ew_ratio.Draw()
qcd_ratio.Draw("same")



ew_m_lo=mcsancPhotonLOEWhist.Clone("ew_m_lo")
uniqueName(ew_m_lo)
ew_m_lo.Add(mcsancPhotonEWhist,mcsancPhotonLOEWhist,1,-1)

k_fakp=mcsancPhotonLOEWhist.Clone("k_fac_p")
uniqueName(k_fakp)
k_fakp.Add(ew_m_lo,fewzHist,1,1)


k_fakm=mcsancPhotonLOEWhist.Clone("k_fac_m")
uniqueName(k_fakm)
k_fakm.Multiply(fewzHist,ew_ratio,1,1)

k_fakm.SetLineColor(ROOT.kRed)

k_fakm.Divide(k_fakm,loHist,1,1)
k_fakp.Divide(k_fakp,loHist,1,1)

k_fakm.Draw()
k_fakp.Draw("same")

k_fak_mean=k_fakp.Clone("mean_k_fak")
uniqueName(k_fak_mean)
k_fak_mean.Reset()
k_fak_mean.Add(k_fakp,k_fakm,1,1)
k_fak_mean.Scale(1./2.)
k_fak_mean.Draw("same")
#raw_input("finished")


leg=ROOT.TLegend(0.67,0.67,0.92,0.92)
k_fakm.GetXaxis().SetRangeUser(200,7000)
ROOT.gPad.SetLogy(0)
ROOT.gPad.SetLogx(1)
k_fakm.GetXaxis().SetTitle("M_{inv} [GeV]")
k_fakm.GetYaxis().SetTitle("#sigma ((N)NLO)/#sigma (LO)")
leg.Clear()
leg.AddEntry(k_fakm,"QCD #otimes EW")
leg.AddEntry(k_fakp,"QCD #oplus EW")
leg.Draw()


outFile=ROOT.TFile("k_fakNNLO_use.root","RECREATE")
k_fakm.Write("k_fakm")
k_fakp.Write("k_fakp")
k_fak_mean.Write("k_fak_mean")
outFile.Close()

NNLO_output.Close()
raw_input("finished")



















