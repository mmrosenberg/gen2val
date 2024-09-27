

def isRecoCCnue(evtTree):

  #require that a reco vertex is found in Fiducial volume
  #(this applies the Wire-Cell 3cm-from-SCE-corrected-detector-boundary cut)
  if evtTree.foundVertex == 0 or evtTree.vtxIsFiducial != 1:
    return False

  #require that spacepoints associated with vertex are not all tagged
  #as cosmics by Wire-Cell charge-light matching algorithms
  if evtTree.vtxFracHitsOnCosmic >= (1. - 1e-6):
    return False

  #require no muon tracks or tracks with a high muon score
  nMuons = 0
  maxMuScore = -20.
  for iT in range(evtTree.nTracks):
    if evtTree.trackIsSecondary[iT] == 1 or evtTree.trackClassified[iT] != 1:
      continue
    if evtTree.trackMuScore[iT] > maxMuScore:
      maxMuScore = evtTree.trackMuScore[iT]
    if evtTree.trackPID[iT] == 13:
      nMuons += 1
  if nMuons > 0 or maxMuScore >= -3.7:
    return False

  #requre at least one electron and that the largest is
  #classified as primary and as e- with high confidence
  nElectrons = 0
  elMaxQ = -1.
  elMaxQConf = -1.
  elMaxQProc = -1
  for iS in range(evtTree.nShowers):
    if evtTree.showerIsSecondary[iS] == 1 or evtTree.showerClassified[iS] != 1:
      continue
    if evtTree.showerPID[iS] == 11:
      nElectrons += 1
      if evtTree.showerCharge[iS] > elMaxQ:
        elMaxQ = evtTree.showerCharge[iS]
        elMaxQConf = evtTree.showerElScore[iS] - (evtTree.showerPhScore[iS] + evtTree.showerPiScore[iS])/2.
        elMaxQProc = evtTree.showerProcess[iS]
  if nElectrons < 1 or elMaxQProc != 0 or elMaxQConf <= 7.1:
    return False
   
  #event passed all cuts if we've gotten here
  return True
 
