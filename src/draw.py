import math
import ROOT
from ROOT import TFile, TH1D, TCanvas, TGraph, TLine, TLegend
import numpy

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)
ROOT.TH1.AddDirectory(False)

gen_file_name = ["gen_default", "gen_simple"]
legend_name = ["pythia8 LO jets", 
    "pythia8 LO partons"]
marker_style = [77, 73]
color = [ROOT.kRed-2, ROOT.kAzure-3]

def draw_cross_section(hist_name, titleX, titleY, xmin, ymin, xmax, ymax):
    canv = TCanvas("a", "a", 600, 900)
    ROOT.gPad.SetLogy()
    frame = ROOT.gPad.DrawFrame(xmin, ymin, xmax, ymax)
    frame.GetXaxis().SetTitle(titleX)
    frame.GetYaxis().SetTitle(titleY)
    frame.GetYaxis().SetTitleOffset(1.3)

    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)

    ROOT.gPad.SetLeftMargin(0.135)
    ROOT.gStyle.SetLegendTextSize(0.05)
    
    frame.Draw("AXIS")
    frame.Draw("SAME AXIS X+ Y+")

    legend = TLegend(0.35, 0.6, 0.9, 0.9)
    legend.SetFillColorAlpha(0, 0.)
    legend.SetLineColorAlpha(0, 0.)

    analytic_input_file = TFile("../output/analytic.root")
    hist_analytic = analytic_input_file.Get(hist_name).Clone("analytic_" + hist_name)
    hist_analytic.SetLineColor(ROOT.kBlack)
    hist_analytic.SetLineWidth(3)
    hist_analytic.Scale(1.1, "nosw2")
    
    legend.AddEntry(hist_analytic, "Analytic", "L")
    hist_analytic.DrawClone("SAME L")

    for i in range (len(gen_file_name)) :
        input_file = TFile("../output/" + gen_file_name[i] + ".root")

        print("Reading file named", gen_file_name[i])
        hist_gen = input_file.Get(hist_name).Clone(gen_file_name[i] + "_" + hist_name)
        hist_gen.SetMarkerStyle(marker_style[i])
        hist_gen.SetMarkerColorAlpha(color[i], 0.7)
        hist_gen.SetLineColorAlpha(color[i], 0.7)
        hist_gen.SetMarkerSize(2.)
        hist_gen.SetLineWidth(2)
        
        legend.AddEntry(hist_gen.Clone(), legend_name[i], "P E")
        hist_gen.DrawClone("SAME E X0")

    legend.DrawClone()

    canv.SaveAs("../output/" + hist_name + ".png")

def draw_ratio(hist_name, titleX, titleY, xmin, ymin, xmax, ymax):
    canv = TCanvas("a", "a", 600, 900)
    ROOT.gPad.SetLogy(0)
    frame = ROOT.gPad.DrawFrame(xmin, ymin, xmax, ymax)
    frame.GetXaxis().SetTitle(titleX)
    frame.GetYaxis().SetTitle(titleY)
    frame.GetYaxis().SetTitleOffset(1.3)

    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.05)

    ROOT.gPad.SetLeftMargin(0.135)
    ROOT.gStyle.SetLegendTextSize(0.05)
    
    frame.Draw("AXIS")
    frame.Draw("SAME AXIS X+ Y+")

    legend = TLegend(0.15, 0.6, 0.5, 0.9)
    legend.SetFillColorAlpha(0, 0.)
    legend.SetLineColorAlpha(0, 0.)

    analytic_input_file = TFile("../output/analytic.root")
    hist_analytic = analytic_input_file.Get(hist_name).Clone("analytic_" + hist_name)
    hist_analytic.Scale(1.1, "nosw2")

    for i in range (len(gen_file_name)) :
        input_file = TFile("../output/" + gen_file_name[i] + ".root")

        print("Reading file named", gen_file_name[i])
        hist_gen = input_file.Get(hist_name).Clone(gen_file_name[i] + "_" + hist_name)
        hist_gen.SetMarkerStyle(marker_style[i])
        hist_gen.SetMarkerColorAlpha(color[i], 0.7)
        hist_gen.SetLineColorAlpha(color[i], 0.7)
        hist_gen.SetMarkerSize(2.)
        hist_gen.SetLineWidth(2)

        for j in range(hist_gen.GetXaxis().GetNbins()):
            div = hist_analytic.GetBinContent(hist_analytic.GetXaxis().FindBin(hist_gen.GetXaxis().GetBinCenter(j)))
            if div > 0:
                hist_gen.SetBinContent(j, hist_gen.GetBinContent(j)/div)
                hist_gen.SetBinError(j, hist_gen.GetBinError(j)/div)
        
        legend.AddEntry(hist_gen.Clone(), legend_name[i], "P E")
        hist_gen.Rebin(2)
        hist_gen.Scale(1./2.)
        hist_gen.DrawClone("SAME E X0")

    legend.DrawClone()

    canv.SaveAs("../output/ratio_" + hist_name + ".png")
    

draw_cross_section("dsigma_dpt", "p_{T}, GeV/c", 
    "d\sigma / dp_{T}, pb \cdot c/GeV^{3}", 25., 1e2, 200., 1e9)
draw_cross_section("dsigma_ddy", "\Delta y", 
    "d \sigma / d \Delta y, pb \cdot c/GeV^{2}", 0., 1e5, 10., 1e10)

draw_ratio("dsigma_dpt", "p_{T}", 
    "d\sigma / dp_{T} \ Pythia8/Analytic", 25., 0., 200., 2.)

draw_ratio("dsigma_ddy", "\Delta y", 
    "d\sigma / d /Delta y \ Pythia8/Analytic", 0., 0., 10., 2.)
