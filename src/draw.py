import ROOT
from ROOT import TFile, TH1, TCanvas, TGraph
import numpy

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)

canv = TCanvas("a", "a", 900, 1200)

ROOT.gPad.SetLogy()
ROOT.gPad.DrawFrame(0., 1e-7, 50., 1e2)

file_name = ["jets0", "jets1"]
color = [ROOT.kAzure-2, ROOT.kRed-2]

for i in range (len(file_name)) :
    input_file = TFile("../output/" + file_name[i] + ".root")

    hist_njets = input_file.Get("jets_multiplicity")
    hist_njets.SetMarkerStyle(0)

    err_njets_low = TGraph()
    err_njets_up = TGraph()

    err_njets_low.SetLineColor(color[i])
    err_njets_up.SetLineColor(color[i])

    err_njets_low.SetLineWidth(3)
    err_njets_up.SetLineWidth(3)

    err_njets_low.SetLineStyle(2)
    err_njets_up.SetLineStyle(2)

    for k in range (hist_njets.GetNbinsX()) :
        hist_njets.SetBinContent(k, hist_njets.GetBinContent(k)/
            (2.*hist_njets.GetBinCenter(k)*hist_njets.GetBinWidth(k)))
        err_njets_low.AddPoint(hist_njets.GetBinCenter(k), 
            hist_njets.GetBinContent(k) - hist_njets.GetBinError(k));
        err_njets_up.AddPoint(hist_njets.GetBinCenter(k), 
            hist_njets.GetBinContent(k) + hist_njets.GetBinError(k));

    hist_njets.Clone().Draw("SAME E3")
    err_njets_low.Clone().Draw("SAME L")
    err_njets_up.Clone().Draw("SAME L")

canv.SaveAs("../output/spectra.png")
