import math
import ROOT
from ROOT import TFile, TH1D, TCanvas, TGraph, TLine, TLegend
import numpy

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)
ROOT.TH1.AddDirectory(False)

gen_file_name = ["gen_default", "gen_simple"]
legend_name = ["pythia8 LO + fastjet3", 
    "#splitline{#splitline{pythia8 LO partons}{Hadronize off}}{#splitline{ISR off}{FSR off}}"]
marker_style = [77, 73]
color = [ROOT.kRed-2, ROOT.kAzure-3]

legend = TLegend(0.53, 0.5, 0.83, 0.9)
legend.SetFillColorAlpha(0, 0.)
legend.SetLineColorAlpha(0, 0.)

def draw_cross_section(hist_name, titleX, titleY, xmin, ymin, xmax, ymax):
    canv = TCanvas("a", "a", 900, 1200)
    ROOT.gPad.SetLogy()

    frame = ROOT.gPad.DrawFrame(xmin, ymin, xmax, ymax)
    frame.GetXaxis().SetTitle(titleX)
    frame.GetYaxis().SetTitle(titleY)
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

        print("Reading file named", gen_file_name)
        gen_dsigma_dpt = input_file.Get(hist_name).Clone(gen_file_name[i] + "_" + hist_name)
        gen_dsigma_dpt.SetMarkerStyle(marker_style[i])
        gen_dsigma_dpt.SetMarkerColorAlpha(color[i], 0.7)
        gen_dsigma_dpt.SetLineColorAlpha(color[i], 0.7)
        gen_dsigma_dpt.SetMarkerSize(2.)
        gen_dsigma_dpt.SetLineWidth(2)

        legend.AddEntry(gen_dsigma_dpt, legend_name[i], "P E")
        gen_dsigma_dpt.Draw("SAME E X0")

    legend.Draw()

    canv.SaveAs("../output/dsigma_dpt.png")

draw_cross_section("dsigma_dpt", "p_{T}, GeV/c", 
    "d\sigma / dp_{T}, pb \cdot c/GeV", 25., 1e2, 200., 1e9)
