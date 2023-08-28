import pyLCIO
import ROOT
import glob

# Set up some options
max_events = -1

# Gather input files
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl_v1/muonGun/reco/*.slcio")

# Set up histograms
hists = {}
hists["pfo_pt"] = ROOT.TH1F("pfo_pt", "pfo_pt", 100, 0, 2000)
hists["muon_pt"] = ROOT.TH1F("muon_pt", "muon_pt", 100, 0, 2000)

# Loop over events
i = 0
for f in fnames:
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(f)

    for event in reader:
        if max_events > 0 and i >= max_events: break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        mcpCollection = event.getCollection("MCParticle")
        pfoCollection = event.getCollection("PandoraPFOs")

        # Loop over the reconstructed objects and fill histograms
        for pfo in pfoCollection:
            pfo_p = pfo.getMomentum()
            pfo_tlv = ROOT.TLorentzVector()
            pfo_tlv.SetPxPyPzE(pfo_p[0], pfo_p[1], pfo_p[2], pfo.getEnergy())
            hists["pfo_pt"].Fill(pfo_tlv.Perp())

            if abs(pfo.getType())==13:
                hists["muon_pt"].Fill(pfo_tlv.Perp())
        i+=1

# Make your plots
for i, h in enumerate(hists):
    c = ROOT.TCanvas("c%i"%i, "c%i"%i)
    hists[h].Draw()
    c.SaveAs("%s.png"%h)
