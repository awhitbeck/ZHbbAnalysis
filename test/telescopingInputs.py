import ROOT
import copy
import math
#from math import *
from array import array

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.16);
ROOT.gStyle.SetPalette(1);


############################################################
############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--reco', action='store_true', dest='reco', default=False, help='no X11 windows')
parser.add_option('--gen', action='store_true', dest='gen', default=False, help='no X11 windows')
(options, args) = parser.parse_args()


def makeROCFromHisto(hsig,hbkg,LtoR):
    
    nbins = hsig.GetNbinsX();
    binsize = hsig.GetBinWidth(1);
    lowedge = hsig.GetBinLowEdge(1);
        
    hsigIntegral = hsig.Integral();
    hbkgIntegral = hbkg.Integral();
    
    xval = array('d', [])
    yval = array('d', [])
    ctr = 0;
    for i in range(1,nbins+1):
        
        effBkg = 0;
        effSig = 0;
        
        if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
        else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;
        
        if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
        else: effSig = hsig.Integral( 1, i )/hsigIntegral;
        
        #print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", effBkg, ", effSig: ", effSig;
        
        xval.append( effSig );
        yval.append( 1-effBkg );
        ctr = ctr + 1;
    
    tg = ROOT.TGraph( nbins, xval, yval );
    return tg;  

if __name__ == '__main__':
    
    
    fns = "";
    fnb = ""; 
    if options.reco: 
        fns = "trees/tree_ZHbb.root";
        fnb = "trees/tree_Zjets.root"; 
    if options.gen: 
        fns = "trees/tree_ZHbb_Gen.root";
        fnb = "trees/tree_Zjets_Gen.root"; 
        

    fileS = ROOT.TFile(fns);
    fileB = ROOT.TFile(fnb);
    tS = fileS.Get("TreeFiller/ZHbbTree");
    tB = fileB.Get("TreeFiller/ZHbbTree");


    xs_bkg = 34.1*1.5*1000.;     #fb, k-factor of 1.5
    xs_sig = 0.415*0.577*1000.;  #fb    
    
    nSig = 99529;
    nBkg = 2626868;
    
    LUMI = 20;
    
    # jet of interest (e.g. ak5);
    joi = 19;
    
    cc_joi_sig = 0;
    cc_joi_bkg = 0;    
    
    mass_lo = 105;
    mass_hi = 135;
    ptcut = 25; 
    
    mpeak = 115;
    
    h_sig_1D = ROOT.TH1F("h_sig_1D",";mass (GeV); N", 50,50,200);
    h_bkg_1D = ROOT.TH1F("h_bkg_1D",";mass (GeV); N", 50,50,200);
    h_sig_2D = ROOT.TH2F("h_sig_2D",";mass (GeV); weight", 50,50,200, 50, 0., 1.);
    h_bkg_2D = ROOT.TH2F("h_bkg_2D",";mass (GeV); weight", 50,50,200, 50, 0., 1.);

    h_wgt_sig = ROOT.TH1F("h_wgt_sig","; weight; N", 50, 0., 1.);
    h_wgt_bkg = ROOT.TH1F("h_wgt_bkg","; weight; N", 50, 0., 1.);    
    
    h_sig_1D_a = ROOT.TH1F("h_sig_1D_a",";mass (GeV); N", 50,50,200);
    h_sig_1D_b = ROOT.TH1F("h_sig_1D_b",";mass (GeV); N", 50,50,200);
    h_sig_1D_c = ROOT.TH1F("h_sig_1D_c",";mass (GeV); N", 50,50,200);
    h_sig_1D_d = ROOT.TH1F("h_sig_1D_d",";mass (GeV); N", 50,50,200);

    h_bkg_1D_a = ROOT.TH1F("h_bkg_1D_a",";mass (GeV); N", 50,50,200);
    h_bkg_1D_b = ROOT.TH1F("h_bkg_1D_b",";mass (GeV); N", 50,50,200);
    h_bkg_1D_c = ROOT.TH1F("h_bkg_1D_c",";mass (GeV); N", 50,50,200);
    h_bkg_1D_d = ROOT.TH1F("h_bkg_1D_d",";mass (GeV); N", 50,50,200);
    
    h_sig_absm = ROOT.TH1F("h_sig_absm",";mass (GeV); N", 100,0,150);
    h_bkg_absm = ROOT.TH1F("h_bkg_absm",";mass (GeV); N", 100,0,150);
                        
    ###################################
    print "tS.GetEntries() = ", tS.GetEntries()                                 
    for i in range(tS.GetEntries()):
        tS.GetEntry(i); 
        if i%10000 == 0: print "i = ", i;
        
        nInterp = tS.Hmass.size();
        nPassed = 0;
        for j in range(nInterp):
            if   j <  25 and j >=  0: h_sig_1D_a.Fill(tS.Hmass[j]);
            elif j <  50 and j >= 25: h_sig_1D_b.Fill(tS.Hmass[j]);
            elif j <  75 and j >= 50: h_sig_1D_c.Fill(tS.Hmass[j]);
            elif j >= 75 and j < 100: h_sig_1D_d.Fill(tS.Hmass[j]);
            else: print "out of range!";
            if tS.Hmass[j] < mass_hi and tS.Hmass[j] > mass_lo and tS.jet1_pt > ptcut and tS.jet2_pt > ptcut: nPassed += 1;
    
        z = nPassed/float(nInterp);
        #print nPassed, nInterp, z
        h_sig_1D.Fill(tS.Hmass[joi]);
        h_sig_2D.Fill(tS.Hmass[joi],z); 
        h_wgt_sig.Fill( z );
        h_sig_absm.Fill( math.fabs( mpeak - tS.Hmass[joi] ) );

    print "tB.GetEntries() = ", tB.GetEntries()                             
    for i in range(tB.GetEntries()):
        tB.GetEntry(i); 
        if i%10000 == 0: print "i = ", i;
        if i == 50000: break;
        
        nInterp = tB.Hmass.size();
        nPassed = 0;
        for j in range(nInterp):
            if   j <  25 and j >=  0: h_bkg_1D_a.Fill(tB.Hmass[j]);
            elif j <  50 and j >= 25: h_bkg_1D_b.Fill(tB.Hmass[j]);
            elif j <  75 and j >= 50: h_bkg_1D_c.Fill(tB.Hmass[j]);
            elif j >= 75 and j < 100: h_bkg_1D_d.Fill(tB.Hmass[j]);
            else: print "out of range!";
            if tB.Hmass[j] < mass_hi and tB.Hmass[j] > mass_lo and tB.jet1_pt > ptcut and tB.jet2_pt > ptcut: nPassed += 1;
    
        z = nPassed/float(nInterp);
        h_bkg_1D.Fill(tB.Hmass[joi]);
        h_bkg_2D.Fill(tB.Hmass[joi],z);         
        h_wgt_bkg.Fill( z );
        h_bkg_absm.Fill( math.fabs( mpeak - tB.Hmass[joi] ) );        
        
    #########################################################
    # save histos
    odir = "";
    if options.reco: odir = "RECO";
    if options.gen: odir = "GEN";
            
    h_sig_1D.SaveAs("histsForAndrew/"+odir+"/h_sig_1D.root");
    h_bkg_1D.SaveAs("histsForAndrew/"+odir+"/h_bkg_1D.root");              
    h_sig_2D.SaveAs("histsForAndrew/"+odir+"/h_sig_2D.root");              
    h_bkg_2D.SaveAs("histsForAndrew/"+odir+"/h_bkg_2D.root"); 

    h_wgt_sig.SaveAs("histsForAndrew/"+odir+"/h_wgt_sig.root"); 
    h_wgt_bkg.SaveAs("histsForAndrew/"+odir+"/h_wgt_bkg.root"); 

    h_sig_1D_a.SaveAs("histsForAndrew/"+odir+"/h_sig_1D_a.root"); 
    h_sig_1D_b.SaveAs("histsForAndrew/"+odir+"/h_sig_1D_b.root"); 
    h_sig_1D_c.SaveAs("histsForAndrew/"+odir+"/h_sig_1D_c.root"); 
    h_sig_1D_d.SaveAs("histsForAndrew/"+odir+"/h_sig_1D_d.root"); 
    h_bkg_1D_a.SaveAs("histsForAndrew/"+odir+"/h_bkg_1D_a.root"); 
    h_bkg_1D_b.SaveAs("histsForAndrew/"+odir+"/h_bkg_1D_b.root"); 
    h_bkg_1D_c.SaveAs("histsForAndrew/"+odir+"/h_bkg_1D_c.root"); 
    h_bkg_1D_d.SaveAs("histsForAndrew/"+odir+"/h_bkg_1D_d.root"); 

    h_sig_absm.SaveAs("histsForAndrew/"+odir+"/h_sig_absm.root"); 
    h_bkg_absm.SaveAs("histsForAndrew/"+odir+"/h_bkg_absm.root"); 

    # draw histos
    sigHistos = [h_sig_1D,h_sig_2D,h_wgt_sig,h_sig_1D_a,h_sig_1D_b,h_sig_1D_c,h_sig_1D_d,h_sig_absm];
    bkgHistos = [h_bkg_1D,h_bkg_2D,h_wgt_bkg,h_bkg_1D_a,h_bkg_1D_b,h_bkg_1D_c,h_bkg_1D_d,h_bkg_absm];       
        
    for i in range(len(bkgHistos)):
        bkgHistos[i].SetLineColor(2);
        sigHistos[i].Scale(1./sigHistos[i].Integral());
        if bkgHistos[i].Integral() > 0: bkgHistos[i].Scale(1./bkgHistos[i].Integral());
        
    c_1D = ROOT.TCanvas("c_1D","c_1D",1000,800);
    h_sig_1D.Draw();
    h_bkg_1D.Draw("sames");
    c_1D.SaveAs("histsForAndrew/"+odir+"/c_1D.eps");
    
    c_2D_sig = ROOT.TCanvas("c_2D_sig","c_2D_sig",1000,800);
    h_sig_2D.SetMaximum(0.25);
    h_sig_2D.Draw("colz");
    ROOT.gPad.SetLogz();
    c_2D_sig.SaveAs("histsForAndrew/"+odir+"/c_2D_sig.eps");

    c_2D_bkg = ROOT.TCanvas("c_2D_bkg","c_2D_bkg",1000,800);
    h_bkg_2D.SetMaximum(0.25);
    h_bkg_2D.Draw("colz");
    ROOT.gPad.SetLogz();
    c_2D_bkg.SaveAs("histsForAndrew/"+odir+"/c_2D_bkg.eps");

    c_wgt = ROOT.TCanvas("c_wgt","c_wgt",1000,800);
    h_wgt_bkg.Draw();
    h_wgt_sig.Draw("sames");
    ROOT.gPad.SetLogy();
    c_wgt.SaveAs("histsForAndrew/"+odir+"/c_wgt.eps");
    
    h_sig_1D_a.SetLineColor(1);
    h_sig_1D_b.SetLineColor(2);
    h_sig_1D_c.SetLineColor(4);
    h_sig_1D_d.SetLineColor(6);

    h_bkg_1D_a.SetLineColor(1);
    h_bkg_1D_b.SetLineColor(2);
    h_bkg_1D_c.SetLineColor(4);
    h_bkg_1D_d.SetLineColor(6);
                
    c_1D_abcd_sig = ROOT.TCanvas("c_1D_abcd_sig","c_1D_abcd_sig",1000,800);
    h_sig_1D_a.Draw();
    h_sig_1D_b.Draw("sames");
    h_sig_1D_c.Draw("sames");
    h_sig_1D_d.Draw("sames");        
    c_1D_abcd_sig.SaveAs("histsForAndrew/"+odir+"/c_1D_abcd_sig.eps");

    c_1D_abcd_bkg = ROOT.TCanvas("c_1D_abcd_bkg","c_1D_abcd_bkg",1000,800);
    h_bkg_1D_a.Draw();
    h_bkg_1D_b.Draw("sames");
    h_bkg_1D_c.Draw("sames");
    h_bkg_1D_d.Draw("sames");        
    c_1D_abcd_bkg.SaveAs("histsForAndrew/"+odir+"/c_1D_abcd_bkg.eps");
        
    c_absm = ROOT.TCanvas("c_absm","c_absm",1000,800);
    h_sig_absm.Draw();
    h_bkg_absm.Draw("sames");    
    c_absm.SaveAs("histsForAndrew/"+odir+"/c_absm.eps");
        
        
    ## make some ROCs
    roc_B = makeROCFromHisto(h_sig_absm,h_bkg_absm,False);
    roc_A = makeROCFromHisto(h_wgt_sig,h_wgt_bkg,True);    

    roc_B.SetLineColor(2);
        
    canRoc = ROOT.TCanvas("canRoc","canRoc",800,800);    
    canRoc.cd();
    hrl = canRoc.DrawFrame(0,0,1.0,1.0);
    hrl.GetXaxis().SetTitle("#epsilon_{sig}");
    hrl.GetYaxis().SetTitle("1 - #epsilon_{bkg}");    
    roc_B.Draw();
    roc_A.Draw();
    canRoc.SaveAs("histsForAndrew/"+odir+"/canRoc.eps");
    
    
    
    
    
    
    
    
        