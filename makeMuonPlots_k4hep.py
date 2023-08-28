# Note: This script was to test out opening k4hep files with ROOT. An actual example can be found here:
# https://gitlab.desy.de/ftx-sft-key4hep/tutorials/-/blob/main/edm4hep_analysis/edm4hep_python.ipynb

# Before running this script be sure to run:
# source /cvmfs/ilc.desy.de/key4hep/setup.sh

# Even with this reader I still couldn't get this to work. RDataFrame seems to be the only thing that lets you do more than simple Draw() commands.

import glob
import ROOT
from ROOT import edm4hep
from podio.root_io import Reader
#from ROOT.edm4hep import utils
#USE_ENERGY = utils.detail.UseEnergyTag()
#from edm4hep_path import get_edm4hep_path


# Load some key libraries
#ROOT.gSystem.Load("libedm4hep")
#ROOT.gSystem.Load("libpodio")

# Set up some options
max_events = 1

# Open the edm4hep files with ROOT
fnames = glob.glob("/collab/project/snowmass21/data/muonc/fmeloni/DataMuC_MuColl_v1/muonGun/edm4hep/*.root")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl_v1/muonGun/edm4hep/*.root")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/photonGun_1000/reco_k4/*.root")
#t = ROOT.TChain("events")
#for f in fnames: t.Add(f)

#reader = Reader(fnames)
#for i, event in enumerate(reader.get('events')):
#    print(i)



# Trying opening with a TFile instead of a TChain in case this is the issue (it's not)
#print(fnames[0])
#f = ROOT.TFile.Open(fnames[0], "READ")
#t = f.Get("events")

# Trying opening with RDataFrames
rdf = ROOT.RDataFrame("events", fnames[0])

#h = rdf.Histo1D("MCParticle.PDG")
#h = rdf.Histo1D("AllTracks.trackerHits_begin")
#h = rdf.Histo1D("AllTracks.getTrackerHits()")
h.Draw()
input("...")

# Check how many events are in the file
#n_events = t.GetEntries()
#print("Found %i events."%n_events)

# This works fine
#t.Draw("MCParticle.PDG")
#input("...")

# Can successfully Draw() but running into issues when doing any kind of loop. Not sure of the cause.
# Loop over events
#for e in t: print("!")

'''
for i, event in enumerate(t):
    print(i)
    #if max_events > 0 and i >= max_events: break
    #if i%100 == 0: print("Processing event %i."%i)
'''



