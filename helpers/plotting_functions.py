
import ROOT as rt

def getMaxBinContent(h):
  hMax = -99.
  for b in range(h.GetNbinsX()):
    if h.GetBinContent(b+1) > hMax:
      hMax = h.GetBinContent(b+1)
  return hMax

def sortHists(hlist):
  iMax = -1
  maxBin = -99.
  for i, h in enumerate(hlist):
    hBinMax = getMaxBinContent(h)
    if hBinMax > maxBin:
      maxBin = hBinMax
      iMax = i
  sortedList = [hlist[iMax]]
  for i, h in enumerate(hlist):
    if i != iMax:
      sortedList.append(h)
  return sortedList

#Must draw histogram before calling this function
def getOverflowLabel(h):
  x = h.GetBinCenter(h.GetNbinsX()) - 0.1*h.GetBinWidth(h.GetNbinsX())
  rt.gPad.Update()
  y = -rt.gPad.GetFrame().GetY2()/27.
  label = rt.TText(x,y,"overflow")
  label.SetTextSize(0.028)
  label.SetTextAngle(-40)
  return label

def getUnderflowLabel(h):
  x = h.GetBinCenter(1) - 0.1*h.GetBinWidth(h.GetNbinsX())
  rt.gPad.Update()
  y = -rt.gPad.GetFrame().GetY2()/30.
  label = rt.TText(x,y,"underflow")
  label.SetTextSize(0.028)
  label.SetTextAngle(-40)
  return label
