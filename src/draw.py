import math
import ROOT
from ROOT import TFile, TH1D, TCanvas, TGraph, TLine, TLegend
import numpy

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)
ROOT.TH1.AddDirectory(False)

gen_file_name = ["gen_default", "gen_nohr"]
hist_name_pt = ["jet_mult_pt", "part_mult_pt"]
legend_name = ["pythia8 LO + fastjet3", 
        "#splitline{#splitline{pythia8 LO partons}{Hadronize off}}{#splitline{ISR off}{FSR off}}"]
marker_style = [77, 73]
color = [ROOT.kRed-2, ROOT.kAzure-3]

min_pt = 25.
max_pt = 200.
min_dsigma_dpt = 1e-1
max_dsigma_dpt = 1e6

min_ddy = 0.
max_ddy = 10.
min_dsigma_ddy = 1e-5
max_dsigma_ddy = 1e9

legend = TLegend(0.53, 0.5, 0.83, 0.9)
legend.SetFillColorAlpha(0, 0.)
legend.SetLineColorAlpha(0, 0.)

canv = TCanvas("a", "a", 900, 1200)
ROOT.gPad.SetLogy()

frame = ROOT.gPad.DrawFrame(min_pt, min_dsigma_dpt, max_pt, max_dsigma_dpt)
frame.GetYaxis().SetTitle("d\sigma / dp_{T}, pb \cdot c/GeV")
frame.GetXaxis().SetTitle("p_{T}, GeV/c")
frame.GetYaxis().SetTitleOffset(1.2)
frame.Draw("AXIS")
frame.Draw("SAME AXIS X+ Y+")

analytic_input_file = TFile("../output/analytic.root")
analytic_dsigma_dpt = analytic_input_file.Get("dsigma_dpt").Clone("analytic_dsigma_dpt")
analytic_dsigma_dpt.SetLineColor(ROOT.kBlack)
analytic_dsigma_dpt.SetLineStyle(2)
analytic_dsigma_dpt.SetLineWidth(3)

legend.AddEntry(analytic_dsigma_dpt, "Analytic", "L")
analytic_dsigma_dpt.Draw("SAME L")

for i in range (len(gen_file_name)) :
    input_file = TFile("../output/" + gen_file_name[i] + ".root")

    print("Reading histogram named", hist_name_pt[i])
    gen_dsigma_dpt = input_file.Get(hist_name_pt[i]).Clone(gen_file_name[i] + "_" + "dsigma_dpt")
    gen_dsigma_dpt.SetMarkerStyle(marker_style[i])
    gen_dsigma_dpt.SetMarkerColorAlpha(color[i], 0.7)
    gen_dsigma_dpt.SetLineColorAlpha(color[i], 0.7)
    gen_dsigma_dpt.SetMarkerSize(2.)
    gen_dsigma_dpt.SetLineWidth(2)
    gen_dsigma_dpt.Rebin(5)
    
    for k in range (gen_dsigma_dpt.GetNbinsX()+1) :
        norm = 2.*math.pi*gen_dsigma_dpt.GetBinCenter(k)*gen_dsigma_dpt.GetBinWidth(k)
        gen_dsigma_dpt.SetBinContent(k, gen_dsigma_dpt.GetBinContent(k)/norm)
        gen_dsigma_dpt.SetBinError(k, gen_dsigma_dpt.GetBinError(k)/norm)

    legend.AddEntry(gen_dsigma_dpt, legend_name[i], "P E")
    gen_dsigma_dpt.Draw("SAME E X0")

legend.Draw()

canv.SaveAs("../output/dsigma_dpt.png")
