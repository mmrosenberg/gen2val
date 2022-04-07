
import ROOT as rt

def sortHists(hlist):
  iMax = -1
  maxBin = -99.
  i = 0 
  for h in hlist:
    for b in range(h.GetNbinsX()):
      if h.GetBinContent(b+1) > maxBin:
        maxBin = h.GetBinContent(b+1)
        iMax = i 
    i = i+1 
  sortedList = [hlist[iMax]]
  for i in range(len(hlist)):
    if i != iMax:
      sortedList.append(hlist[i])
  return sortedList

