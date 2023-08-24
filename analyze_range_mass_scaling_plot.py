
import ROOT as rt
import numpy as np

import matplotlib.pyplot as plt

from helpers.pionEnergyEstimator import pionRange2T

#from scipy.stats import linregress
#from scipy.interpolate import interp1d
#
#class pionRange2T:
#
#  def __init__(self, r, T):
#    res = linregress(r[60:80], T[60:80])
#    self.slope = res.slope
#    self.intercept = res.intercept
#    self.f_int = interp1d(r[:60], T[:60])
#    self.range = r
#    self.KE = T
#
#  def Eval(self, length):
#    if length < self.range[60]:
#      return self.f_int(length)
#    return self.intercept + self.slope*length


f = rt.TFile("analyze_range_mass_scaling_process_output.root")
t = f.Get("ParticleTree")

h_pi = rt.TH2F("h_pi","h_pi",800,0,800,2000,0,2000)
h_mu = rt.TH2F("h_mu","h_mu",800,0,800,2000,0,2000)

t.Draw("KE:length >> h_pi","pdg == 211")
t.Draw("KE:length >> h_mu","pdg == 13")

range_mu = []
KE_mu = []
range_pi = []
KE_pi = []

for ix in range(1,801):
  range_pi.append(0.5 + (ix-1))
  range_mu.append(0.5 + (ix-1))
  found_pi_val = False
  found_mu_val = False
  for iy in range(1,2001):
    if h_pi.GetBinContent(ix,iy) > 1e-3 and not found_pi_val:
      KE_pi.append(0.5 + (iy-1))
      found_pi_val = True
    if h_mu.GetBinContent(ix,iy) > 1e-3 and not found_mu_val:
      KE_mu.append(0.5 + (iy-1))
      found_mu_val = True
    if found_pi_val and found_mu_val:
      break
  if not found_pi_val:
    KE_pi.append(0)
  if not found_mu_val:
    KE_mu.append(0)

fspline = rt.TFile("/home/matthew/microboone/tufts/ubdl/larflow/larflow/Reco/data/Proton_Muon_Range_dEdx_LAr_TSplines.root")
muSpline = fspline.Get("sMuonRange2T")
prSpline = fspline.Get("sProtonRange2T")
KE_mu_spline = [muSpline.Eval(range_mu[i]) for i in range(800)]
KE_pi_spline1 = [muSpline.Eval(range_pi[i])*1.32 for i in range(800)]
KE_pi_spline2 = [muSpline.Eval(range_pi[i]*1.32) for i in range(800)]

ratio = [KE_pi[i]/KE_mu[i] if (KE_mu[i] > 0) else 0. for i in range(800)]

#slope, intercept, r, p, se = linregress(range_pi[60:80], KE_pi[60:80])
#piKEestimator = pionRange2T(range_pi, KE_pi)
piKEestimator = pionRange2T()
KE_pi_spline3 = [piKEestimator.Eval(range_pi[i]) for i in range(800)]

#print("pion KE estimator results:")
#print("cutoff: ", piKEestimator.range[60])
#print("regression slope: ", piKEestimator.slope)
#print("regression intercept: ", piKEestimator.intercept)
#print("interpolation range points:")
#print(piKEestimator.range[:60])
#print("interpolation KE points:")
#print(piKEestimator.KE[:60])


print("range, Emu, Epi, Epr:")
for t in range(1,11,1):
  try:
    r = 0.3*t
    print(r, muSpline.Eval(r), piKEestimator.Eval(r), prSpline.Eval(r))
  except: 
    continue


g_pi = rt.TGraph(800,np.array(range_pi),np.array(KE_pi))
g_mu = rt.TGraph(800,np.array(range_mu),np.array(KE_mu))
g_pi_spline1 = rt.TGraph(800,np.array(range_pi),np.array(KE_pi_spline1))
g_pi_spline2 = rt.TGraph(800,np.array(range_pi),np.array(KE_pi_spline2))
g_pi_spline3 = rt.TGraph(800,np.array(range_pi),np.array(KE_pi_spline3))
g_mu_spline = rt.TGraph(800,np.array(range_mu),np.array(KE_mu_spline))
g_pi.SetLineColor(rt.kBlue)
g_mu.SetLineColor(rt.kBlue)
g_pi_spline1.SetLineColor(rt.kBlue)
g_pi_spline2.SetLineColor(rt.kRed)
g_pi_spline3.SetLineColor(8)
g_mu_spline.SetLineColor(rt.kBlue)

fout = rt.TFile("foo.root","RECREATE")
h_pi.Write()
h_mu.Write()

cnv_pi = rt.TCanvas("cnv_pi","cnv_pi")
h_pi.Draw()
g_pi.Draw("LSAME")
cnv_pi.Write()

cnv_mu = rt.TCanvas("cnv_mu","cnv_mu")
h_mu.Draw()
g_mu.Draw("LSAME")
cnv_mu.Write()

cnv_pi_spline = rt.TCanvas("cnv_pi_spline","cnv_pi_spline")
h_pi.Draw()
g_pi_spline1.Draw("LSAME")
g_pi_spline2.Draw("LSAME")
g_pi_spline3.Draw("LSAME")
#g_pi.Draw("LSAME")
cnv_pi_spline.Write()

cnv_mu_spline = rt.TCanvas("cnv_mu_spline","cnv_mu_spline")
h_mu.Draw()
g_mu_spline.Draw("LSAME")
cnv_mu_spline.Write()


plt.figure(1)
plt.plot(range_mu[5:85], KE_mu_spline[5:85], 'k-')
plt.plot(range_mu[5:85], KE_mu[5:85], 'b-')
plt.plot(range_pi[5:85], KE_pi[5:85], 'r-')

plt.figure(2)
plt.plot(range_pi[5:85], ratio[5:85], 'g-')

plt.show()


