
import ROOT as rt
from math import sqrt
from helpers.plotting_functions import getOverflowLabel, getUnderflowLabel

#====== flux, xsec, and reint arrays from print_xsec_flux_reint_systematics_individual.py =======#

#Public note results:

#CCnue_flux = [0.1124189613852805, 0.07399749381230854, 0.0530396519260321, 0.05292510563248634, 0.05729280255755563, 0.060194593158459014, 0.06607423469968957, 0.06936226344842056, 0.07858957147804232, 0.08404475729172949, 0.08602330571329511, 0.09536612951237457, 0.10914361624114102, 0.1021708648605673]
#CCnue_xsec = [0.1700153189423411, 0.1437454547043476, 0.12336007113244706, 0.1272410227471294, 0.12954265784821126, 0.13410204573622445, 0.1394031881314876, 0.1327611328744459, 0.13691461647578157, 0.13996966217204238, 0.13417285742316373, 0.14668012670167851, 0.16184288097727648, 0.1474409297702654]
#CCnue_reint = [0.012318643130717073, 0.018260537876033885, 0.01467647988525375, 0.01593452547260682, 0.015375918649560859, 0.013421706923467133, 0.023523628046220678, 0.023014133310368615, 0.018362281410625778, 0.011090737902148241, 0.03769933721116542, 0.009752066830373, 0.01697418711664474, 0.010439961712803071]

#CCnumu_flux = [0.07093712723510005, 0.06547055442149248, 0.06070454176956921, 0.060888754288359756, 0.06343960805555783, 0.06725119888328036, 0.07617871616624353, 0.08493111228768717, 0.09336853048555971, 0.1060081169681605, 0.11358397326452797, 0.11741855042353896, 0.12605248460909774, 0.12959167569465166, 0.1267479205277652, 0.1334334184118158, 0.12069397612602212, 0.11377468767560821, 0.10393028546287716, 0.11445347711927684, 0.07480653818286226]
#CCnumu_xsec = [0.18753090614405046, 0.15498954518547098, 0.12677848981879022, 0.12570561023838936, 0.11776337646230818, 0.13151826378190523, 0.12528865312744172, 0.1409540080036013, 0.13443289697536073, 0.14251637463227038, 0.14206226494189247, 0.14361715698060884, 0.15539852153502753, 0.1432905106421294, 0.14471624683101247, 0.1707787312362085, 0.15800899354011375, 0.1511235698621501, 0.16346222905009822, 0.18297949955709408, 0.15508667263573556]
#CCnumu_reint = [0.0074689601421568845, 0.014016836092753268, 0.01946445836213243, 0.019496155821495817, 0.014398907462933017, 0.013315147844345395, 0.010893695014284917, 0.01103302906677802, 0.009367169671034464, 0.014973552330555955, 0.01558965477268725, 0.012575089971731549, 0.027523813470056793, 0.01629275274499139, 0.014108078159999331, 0.020480043152087013, 0.02668709065685045, 0.031228483366185326, 0.032366689541862316, 0.02075493107371821, 0.03468966176878337]

#Final version results:

CCnue_flux = [0.07816731950484176, 0.0530396519260321, 0.05292510563248634, 0.05729280255755563, 0.060194593158459014, 0.06607423469968957, 0.06936226344842056, 0.07858957147804232, 0.09027455167909595]
CCnue_xsec = [0.1442938677895219, 0.12336007113244706, 0.1272410227471294, 0.12954265784821126, 0.13410204573622445, 0.1394031881314876, 0.1327611328744459, 0.13691461647578157, 0.13695075468546777]
CCnue_reint = [0.017221390152121454, 0.01467647988525375, 0.01593452547260682, 0.015375918649560859, 0.013421706923467133, 0.023523628046220678, 0.023014133310368615, 0.018362281410625778, 0.017410170287470772]

CCnumu_flux = [0.07093712723510005, 0.06547055442149248, 0.06070454176956921, 0.060888754288359756, 0.06343960805555783, 0.06725119888328036, 0.07617871616624353, 0.08493111228768717, 0.09336853048555971, 0.1060081169681605, 0.11358397326452797, 0.11741855042353896, 0.12605248460909774, 0.12959167569465166, 0.1267479205277652, 0.1334334184118158, 0.12069397612602212, 0.11377468767560821, 0.08238989256766768]
CCnumu_xsec = [0.18753090614405046, 0.15498954518547098, 0.12677848981879022, 0.12570561023838936, 0.11776337646230818, 0.13151826378190523, 0.12528865312744172, 0.1409540080036013, 0.13443289697536073, 0.14251637463227038, 0.14206226494189247, 0.14361715698060884, 0.15539852153502753, 0.1432905106421294, 0.14471624683101247, 0.1707787312362085, 0.15800899354011375, 0.1511235698621501, 0.14362094259545946]
CCnumu_reint = [0.0074689601421568845, 0.014016836092753268, 0.01946445836213243, 0.019496155821495817, 0.014398907462933017, 0.013315147844345395, 0.010893695014284917, 0.01103302906677802, 0.009367169671034464, 0.014973552330555955, 0.01558965477268725, 0.012575089971731549, 0.027523813470056793, 0.01629275274499139, 0.014108078159999331, 0.020480043152087013, 0.02668709065685045, 0.031228483366185326, 0.02679460131969319]


#====== detvar arrays from combine_detvar_matrices.py ======#

#Public note results:
#CCnue_detvar = [0.233808, 0.114968, 0.093613, 0.054833, 0.071288, 0.039865, 0.076471, 0.091423, 0.190309, 0.121863, 0.210883, 0.181639, 0.221191, 0.281699]
#CCnumu_detvar = [0.171868, 0.058173, 0.043370, 0.024423, 0.059679, 0.040625, 0.042156, 0.052031, 0.075653, 0.157159, 0.225639, 0.127843, 0.150389, 0.181392, 0.075954, 0.075954, 0.075954, 0.075954, 0.075954, 0.075954, 0.075954]

#Original low stats, no hacks results for public note comparison with/without bkg hacks:
#CCnue_detvar_noHacks = [0.804726, 0.365975, 0.215028, 0.097052, 0.202427, 0.187942, 0.192515, 0.383597, 0.233973, 0.496203, 0.220667, 0.192027, 1.389503, 0.614191]
#CCnumu_detvar_noHacks = [0.171868, 0.058173, 0.043370, 0.024423, 0.059679, 0.040625, 0.042156, 0.052031, 0.075653, 0.157159, 0.225639, 0.127843, 0.150389, 0.181392, 0.215316, 0.449730, 0.234092, 0.248301, 0.421059, 0.492921, 0.384670]

#Final version results:
CCnue_detvar = [0.140786, 0.110968, 0.057446, 0.061724, 0.103982, 0.101497, 0.085781, 0.232140, 0.178371]
CCnumu_detvar = [0.051019, 0.027552, 0.018892, 0.018844, 0.020452, 0.015043, 0.021307, 0.019210, 0.027819, 0.038303, 0.057756, 0.047054, 0.051812, 0.154476, 0.087227, 0.143548, 0.108314, 0.108939, 0.091634]


#====== stats errors from print_stats_errors.py ======#

#Public note results:
#CCnue_stats = [0.17379280466201152, 0.05190113892807601, 0.03324943635028995, 0.025044671978854197, 0.03507966276574342, 0.02339101721420732, 0.049990200544990177, 0.038756480513926134, 0.050216214227916105, 0.030719765928555583, 0.1302613064542679, 0.055776056823431085, 0.0714180081171576, 0.05556963982561719]
#CCnumu_stats = [0.029019912897016286, 0.015193939475021349, 0.012492320217419065, 0.011774485267871119, 0.012015665218863936, 0.013192894222836149, 0.01502540265049065, 0.01769387007307403, 0.02073487011499233, 0.024610470400575738, 0.0308879792812562, 0.03649281742147544, 0.04648162671197284, 0.05203735787782406, 0.0651184805998263, 0.0764727903340569, 0.07523620409101196, 0.09375390463215977, 0.0999331777554916, 0.10923255922716828, 0.07029098118379773]

#Final version results:
CCnue_stats = [0.050647620390681374, 0.03324943635028995, 0.025044671978854197, 0.03507966276574342, 0.02339101721420732, 0.049990200544990177, 0.038756480513926134, 0.050216214227916105, 0.039732096533470394]
CCnumu_stats = [0.029019912897016286, 0.015193939475021349, 0.012492320217419065, 0.011774485267871119, 0.012015665218863936, 0.013192894222836149, 0.01502540265049065, 0.01769387007307403, 0.02073487011499233, 0.024610470400575738, 0.0308879792812562, 0.03649281742147544, 0.04648162671197284, 0.05203735787782406, 0.0651184805998263, 0.0764727903340569, 0.07523620409101196, 0.09375390463215977, 0.053340137632566974]



#======= make plots =======#

rt.gStyle.SetOptStat(0)

#Public note binning:
#nueB_n = 14
#nueB_l = 0.
#nueB_h = 2.8
#numuB_n = 21
#numuB_l = 0.
#numuB_h = 2.1

#Final version binning:
nueB_n = 9
nueB_l = 0.2
nueB_h = 2.0
numuB_n = 19
numuB_l = 0.
numuB_h = 1.9

h_CCnue_flux = rt.TH1F("h_CCnue_flux", "Uncertainties in CCnue Selection", nueB_n, nueB_l, nueB_h)
h_CCnue_xsec = rt.TH1F("h_CCnue_xsec", "Uncertainties in CCnue Selection", nueB_n, nueB_l, nueB_h)
h_CCnue_reint = rt.TH1F("h_CCnue_reint", "Uncertainties in CCnue Selection", nueB_n, nueB_l, nueB_h)
h_CCnumu_flux = rt.TH1F("h_CCnumu_flux", "Uncertainties in CCnumu Selection", numuB_n, numuB_l, numuB_h)
h_CCnumu_xsec = rt.TH1F("h_CCnumu_xsec", "Uncertainties in CCnumu Selection", numuB_n, numuB_l, numuB_h)
h_CCnumu_reint = rt.TH1F("h_CCnumu_reint", "Uncertainties in CCnumu Selection", numuB_n, numuB_l, numuB_h)
h_CCnue_detvar = rt.TH1F("h_CCnue_detvar", "Uncertainties in CCnue Selection", nueB_n, nueB_l, nueB_h)
#h_CCnue_detvar_noHacks = rt.TH1F("h_CCnue_detvar_noHacks", "Detector Systematic Uncertainties in CCnue Selection", nueB_n, nueB_l, nueB_h)
h_CCnumu_detvar = rt.TH1F("h_CCnumu_detvar", "Uncertainties in CCnumu Selection", numuB_n, numuB_l, numuB_h)
#h_CCnumu_detvar_noHacks = rt.TH1F("h_CCnumu_detvar_noHacks", "Detector Systematic Uncertainties in CCnumu Selection", numuB_n, numuB_l, numuB_h)
h_CCnue_stats = rt.TH1F("h_CCnue_stats", "Uncertainties in CCnue Selection", nueB_n, nueB_l, nueB_h)
h_CCnumu_stats = rt.TH1F("h_CCnumu_stats", "Uncertainties in CCnumu Selection", numuB_n, numuB_l, numuB_h)
h_CCnue_tot = rt.TH1F("h_CCnue_tot", "Uncertainties in CCnue Selection", nueB_n, nueB_l, nueB_h)
h_CCnumu_tot = rt.TH1F("h_CCnumu_tot", "Uncertainties in CCnumu Selection", numuB_n, numuB_l, numuB_h)

h_CCnue_tot.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
h_CCnue_tot.GetYaxis().SetTitle("fractional uncertainty")
h_CCnue_tot.GetXaxis().CenterTitle(True);
h_CCnue_tot.SetLineWidth(2)
h_CCnue_tot.SetLineColor(rt.kBlack)
h_CCnue_stats.SetLineWidth(2)
h_CCnue_stats.SetLineColor(40)
h_CCnue_detvar.SetLineWidth(2)
h_CCnue_detvar.SetLineColor(rt.kBlue)
#h_CCnue_detvar_noHacks.SetLineWidth(2)
#h_CCnue_detvar_noHacks.SetLineColor(rt.kRed)
#h_CCnue_detvar_noHacks.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
#h_CCnue_detvar_noHacks.GetYaxis().SetTitle("fractional uncertainty")
#h_CCnue_detvar_noHacks.GetXaxis().CenterTitle(True);
h_CCnue_flux.SetLineWidth(2)
h_CCnue_flux.SetLineColor(8)
h_CCnue_reint.SetLineWidth(2)
h_CCnue_reint.SetLineColor(12)
h_CCnue_xsec.SetLineWidth(2)
h_CCnue_xsec.SetLineColor(rt.kRed)
h_CCnumu_tot.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
h_CCnumu_tot.GetYaxis().SetTitle("fractional uncertainty")
h_CCnumu_tot.GetXaxis().CenterTitle(True);
h_CCnumu_tot.SetLineWidth(2)
h_CCnumu_tot.SetLineColor(rt.kBlack)
h_CCnumu_stats.SetLineWidth(2)
h_CCnumu_stats.SetLineColor(40)
h_CCnumu_detvar.SetLineWidth(2)
h_CCnumu_detvar.SetLineColor(rt.kBlue)
#h_CCnumu_detvar_noHacks.SetLineWidth(2)
#h_CCnumu_detvar_noHacks.SetLineColor(rt.kRed)
#h_CCnumu_detvar_noHacks.GetXaxis().SetTitle("reconstructed neutrino energy (GeV)")
#h_CCnumu_detvar_noHacks.GetYaxis().SetTitle("fractional uncertainty")
#h_CCnumu_detvar_noHacks.GetXaxis().CenterTitle(True);
h_CCnumu_flux.SetLineWidth(2)
h_CCnumu_flux.SetLineColor(8)
h_CCnumu_reint.SetLineWidth(2)
h_CCnumu_reint.SetLineColor(12)
h_CCnumu_xsec.SetLineWidth(2)
h_CCnumu_xsec.SetLineColor(rt.kRed)


for i in range(1,numuB_n+1,1):
  if i <= nueB_n:
    h_CCnue_flux.SetBinContent(i, CCnue_flux[i-1])
    h_CCnue_xsec.SetBinContent(i, CCnue_xsec[i-1])
    h_CCnue_reint.SetBinContent(i, CCnue_reint[i-1])
    h_CCnue_detvar.SetBinContent(i, CCnue_detvar[i-1])
    #h_CCnue_detvar_noHacks.SetBinContent(i, CCnue_detvar_noHacks[i-1])
    h_CCnue_stats.SetBinContent(i, CCnue_stats[i-1])
    h_CCnue_tot.SetBinContent(i, sqrt(CCnue_flux[i-1]**2 + CCnue_xsec[i-1]**2 + CCnue_reint[i-1]**2 + CCnue_detvar[i-1]**2 + CCnue_stats[i-1]**2))
  h_CCnumu_flux.SetBinContent(i, CCnumu_flux[i-1])
  h_CCnumu_xsec.SetBinContent(i, CCnumu_xsec[i-1])
  h_CCnumu_reint.SetBinContent(i, CCnumu_reint[i-1])
  h_CCnumu_detvar.SetBinContent(i, CCnumu_detvar[i-1])
  #h_CCnumu_detvar_noHacks.SetBinContent(i, CCnumu_detvar_noHacks[i-1])
  h_CCnumu_stats.SetBinContent(i, CCnumu_stats[i-1])
  h_CCnumu_tot.SetBinContent(i, sqrt(CCnumu_flux[i-1]**2 + CCnumu_xsec[i-1]**2 + CCnumu_reint[i-1]**2 + CCnumu_detvar[i-1]**2 + CCnumu_stats[i-1]**2))


outfile = rt.TFile("make_gen2_paper_systematics_plots_output.root","RECREATE")

h_CCnue_tot.GetYaxis().SetRangeUser(0,0.4)
h_CCnue_tot.GetXaxis().SetTitleOffset(1.2)
cnv_CCnue = rt.TCanvas("cnv_CCnue","cnv_CCnue")
h_CCnue_tot.Draw("HIST")
h_CCnue_reint.Draw("HISTSAME")
h_CCnue_stats.Draw("HISTSAME")
h_CCnue_detvar.Draw("HISTSAME")
h_CCnue_flux.Draw("HISTSAME")
h_CCnue_xsec.Draw("HISTSAME")
ovfLabel_CCnue = getOverflowLabel(h_CCnue_tot)
ovfLabel_CCnue.Draw()
unfLabel_CCnue = getUnderflowLabel(h_CCnue_tot)
unfLabel_CCnue.Draw()
leg_CCnue = rt.TLegend(0.194842,0.64916,0.648997,0.846639)
leg_CCnue.SetNColumns(2)
leg_CCnue.AddEntry(h_CCnue_stats, "Statistical", "l")
leg_CCnue.AddEntry(h_CCnue_detvar, "Detector", "l")
leg_CCnue.AddEntry(h_CCnue_reint, "Hadron re-int.", "l")
leg_CCnue.AddEntry(h_CCnue_flux, "Flux", "l")
leg_CCnue.AddEntry(h_CCnue_xsec, "Cross Section", "l")
leg_CCnue.AddEntry(h_CCnue_tot, "Total", "l")
leg_CCnue.Draw()
cnv_CCnue.Write()

h_CCnumu_tot.GetYaxis().SetRangeUser(0,0.4)
h_CCnumu_tot.GetXaxis().SetTitleOffset(1.2)
cnv_CCnumu = rt.TCanvas("cnv_CCnumu","cnv_CCnumu")
h_CCnumu_tot.Draw("HIST")
h_CCnumu_reint.Draw("HISTSAME")
h_CCnumu_stats.Draw("HISTSAME")
h_CCnumu_detvar.Draw("HISTSAME")
h_CCnumu_flux.Draw("HISTSAME")
h_CCnumu_xsec.Draw("HISTSAME")
ovfLabel_CCnumu = getOverflowLabel(h_CCnumu_tot)
ovfLabel_CCnumu.Draw()
leg_CCnumu = rt.TLegend(0.170487,0.659664,0.643266,0.861345)
leg_CCnumu.SetNColumns(2)
leg_CCnumu.AddEntry(h_CCnumu_stats, "Statistical", "l")
leg_CCnumu.AddEntry(h_CCnumu_detvar, "Detector", "l")
leg_CCnumu.AddEntry(h_CCnumu_reint, "Hadron re-int.", "l")
leg_CCnumu.AddEntry(h_CCnumu_flux, "Flux", "l")
leg_CCnumu.AddEntry(h_CCnumu_xsec, "Cross Section", "l")
leg_CCnumu.AddEntry(h_CCnumu_tot, "Total", "l")
leg_CCnumu.Draw()
cnv_CCnumu.Write()


#====== Public note comparison with and without background hacks: ======#

#h_CCnue_detvar_noHacks.GetYaxis().SetRangeUser(0,1.46)
#h_CCnue_detvar_noHacks.GetXaxis().SetTitleOffset(1.2)
#cnv_detvar_CCnue = rt.TCanvas("cnv_detvar_CCnue","cnv_detvar_CCnue")
#h_CCnue_detvar_noHacks.Draw("HIST")
#h_CCnue_detvar.Draw("HISTSAME")
#ovfLabel_CCnue_detvar = getOverflowLabel(h_CCnue_detvar_noHacks)
#ovfLabel_CCnue_detvar.Draw()
#leg_CCnue_detvar = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_CCnue_detvar.AddEntry(h_CCnue_detvar, "With Corrections (lower bound)", "l")
#leg_CCnue_detvar.AddEntry(h_CCnue_detvar_noHacks, "No Corrections (upper bound)", "l")
#leg_CCnue_detvar.Draw()
#cnv_detvar_CCnue.Write()

#h_CCnumu_detvar_noHacks.GetXaxis().SetTitleOffset(1.2)
#cnv_detvar_CCnumu = rt.TCanvas("cnv_detvar_CCnumu","cnv_detvar_CCnumu")
#h_CCnumu_detvar_noHacks.Draw("HIST")
#h_CCnumu_detvar.Draw("HISTSAME")
#ovfLabel_CCnumu_detvar = getOverflowLabel(h_CCnumu_detvar_noHacks)
#ovfLabel_CCnumu_detvar.Draw()
#leg_CCnumu_detvar = rt.TLegend(0.7,0.7,0.9,0.9)
#leg_CCnumu_detvar.AddEntry(h_CCnumu_detvar, "With Corrections (lower bound)", "l")
#leg_CCnumu_detvar.AddEntry(h_CCnumu_detvar_noHacks, "No Corrections (upper bound)", "l")
#leg_CCnumu_detvar.Draw()
#cnv_detvar_CCnumu.Write()


