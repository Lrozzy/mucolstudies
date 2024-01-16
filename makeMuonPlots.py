import pyLCIO
import ROOT
import glob
import json
from math import *
import numpy as np


# ############## SETUP #############################
# Prevent ROOT from drawing while you're running -- good for slow remote servers
# Instead, save files and view them with an sftp client like Fetch (feel free to ask me for my UTK license)
ROOT.gROOT.SetBatch()

# Set up some options
max_events = -1

# Gather input files
# Note: these are using the path convention from the singularity command in the MuCol tutorial (see README)
fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_*/*.slcio") #found a way to run through all directories at once
#fnames = glob.glob("/data/fmeloni/LegacyProductions/before29Jul23/DataMuC_MuColl_v1/muonGun/reco/*.slcio")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/gen_muonGun/recoBIB/*.slcio")
#fnames = glob.glob("/data/fmeloni/DataMuC_MuColl10_v0A/muonGun_1000/recoBIB/*.slcio")
exclude_pattern = "/data/fmeloni/DataMuC_MuColl10_v0A/reco/muonGun_pT_1000_5000/"
fnames = [fname for fname in fnames if not fname.startswith(exclude_pattern)] #excluding pt_1000_5000
print("Found %i files."%len(fnames))

# ############## CREATE EMPTY HISTOGRAM OBJECTS  #############################
# Set up histograms
# This is an algorithmic way of making a bunch of histograms and storing them in a dictionary
variables = {}
variables["pt"] =  {"nbins": 20, "xmin": 0, "xmax": 2000}
variables["eta"] = {"nbins": 50, "xmin": -5, "xmax": 5}
variables["phi"] = {"nbins": 20, "xmin": -3.5, "xmax": 3.5}
variables["n"] = {"nbins": 20, "xmin": 0, "xmax": 20}
hists = {}
for obj in ["pfo", "pfo_mu", "mcp", "mcp_mu", "mcp_mu_match"]:
    for var in variables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, obj+"_"+var, variables[var]["nbins"], variables[var]["xmin"], variables[var]["xmax"])

# Adding histograms for D0 and Z0 resolutions
variables1 = {}
variables1["d0_res"] = {"nbins": 100, "xmin": -0.1, "xmax": 0.1}
variables1["z0_res"] = {"nbins": 100, "xmin": -0.1, "xmax": 0.1}
for var in ["d0_res", "z0_res"]:
    hists[var] = ROOT.TH1F(var, var, variables1[var]["nbins"], variables1[var]["xmin"], variables1[var]["xmax"])

# Adding histograms for pT resolution, reconstructed hit counts
variables1["nhits"] = {"nbins": 50, "xmin": 0, "xmax": 50}
variables1["pt_res_hits"] = {"nbins": 100, "xmin": -0.5, "xmax": 0.5}
for var in ["nhits", "pt_res_hits"]:
    hists[var] = ROOT.TH1F(var, var, variables1[var]["nbins"], variables1[var]["xmin"], variables1[var]["xmax"])

# Adding 2D histograms
hists_2d = {}
variables_2d = {
    "d0_res_vs_pt": {"xbins": 100, "xmin": 0, "xmax": 2000, "ybins": 100, "ymin": -0.1, "ymax": 0.1,
                    "ylabel": "D0 Resolution", "xlabel": "pT"},
    "d0_res_vs_eta": {"xbins": 100, "xmin": -3, "xmax": 3, "ybins": 100, "ymin": -0.1, "ymax": 0.1,
                    "ylabel": "D0 Resolution", "xlabel": "Eta"},
    "z0_res_vs_pt": {"xbins": 100, "xmin": 0, "xmax": 2000, "ybins": 100, "ymin": -0.1, "ymax": 0.1,
                    "ylabel": "Z0 Resolution", "xlabel": "pT"},
    "z0_res_vs_eta": {"xbins": 100, "xmin": -3, "xmax": 3, "ybins": 100, "ymin": -0.1, "ymax": 0.1,
                    "ylabel": "Z0 Resolution", "xlabel": "Eta"},
    "pt_res_vs_eta": {"xbins": 100, "xmin": -3, "xmax": 3, "ybins": 100, "ymin": -0.5, "ymax": 0.5,
                    "ylabel": "pT Resolution", "xlabel": "Eta"},
    "pt_res_vs_pt": {"xbins": 100, "xmin": 0, "xmax": 2000, "ybins": 100, "ymin": -0.5, "ymax": 0.5,
                    "ylabel": "pT Resolution", "xlabel": "pT"},
}

for var in variables_2d:
    hists_2d[var] = ROOT.TH2F(var, var, variables_2d[var]["xbins"], variables_2d[var]["xmin"],
                              variables_2d[var]["xmax"], variables_2d[var]["ybins"], variables_2d[var]["ymin"],
                              variables_2d[var]["ymax"])


# Making a separate set of binning conventions for plots showing resolutions
# these plots will all be filled with the difference between a pfo and a mcp object value
dvariables = {}
dvariables["dpt"] =     {"nbins": 100, "xmin": -500, "xmax": 500}
dvariables["drelpt"] =  {"nbins": 100, "xmin": -0.5, "xmax": 0.5}
dvariables["dphi"] =    {"nbins": 100, "xmin": -0.001, "xmax": 0.001}
dvariables["deta"] =    {"nbins": 100, "xmin": -0.001, "xmax": 0.001}
for obj in ["d_mu"]:
    for var in dvariables:
        hists[obj+"_"+var] = ROOT.TH1F(obj+"_"+var, obj+"_"+var, dvariables[var]["nbins"], dvariables[var]["xmin"], dvariables[var]["xmax"])

# Finally making one 2D histogram non-algorithmically; this is what I'll use for a
# pT resolution vs. pT plot.
h_2d_relpt = ROOT.TH2F("h_2d_relpt", "h_2d_relpt", 20, 0, 1000, 500, -0.5, 0.5)

def check_hard_radiation(mcp, fractional_threshold):
    had_hard_rad = False
    daughters = mcp.getDaughters() 
    for d in daughters:
        if d.getPDG() == 22 or d.getPDG() == 23 or d.getPDG() == 24:        
            if d.getEnergy() > fractional_threshold*mcp.getEnergy():
                had_hard_rad = True   
    return had_hard_rad

# Create empty lists for each variable
mcp_pt = [] #mcp = MCParticle (truth)
mcp_phi = []
mcp_eta = []

pfo_pt = [] #pfo = Particle FLow Object (reconstructed)
pfo_phi = []
pfo_eta = []
pfo_mu_pt = []
pfo_mu_phi = []
pfo_mu_eta = []

mcp_mu_pt = [] #mcp_mu = MCParticle muon
mcp_mu_phi = []
mcp_mu_eta = []

mcp_mu_match_pt = [] #mcp_mu_match = MCParticle muon that was matched to a PFO muon
mcp_mu_match_phi = []
mcp_mu_match_eta = []

d_mu_dpt = [] #d_mu = difference between PFO muon and MCParticle muon
d_mu_drelpt = []
d_mu_dphi = []
d_mu_deta = []

# TRACK
d0_res = [] #d0_res = d0 resolution
d0_res_match = [] #d0_res_match = d0 resolution for matched muons
z0_res = [] #z0_res = z0 resolution
z0_res_match = [] #z0_res_match = z0 resolution for matched muons

nhits = []
pixel_nhits = []
inner_nhits = []
outer_nhits = []
pt_res_hits = []

d0_res_vs_pt = [] 
d0_res_vs_eta = []
z0_res_vs_pt = []
z0_res_vs_eta = []
pt_res_vs_eta = [] #track muon pt resolution vs pt
pt_res_vs_pt = []
pt_res = [] 

# Truth matched 
pt_match = [] #This is truth pt
track_pt = [] #This is track pt
track_eta = [] #This is track eta
eta_match = [] #This is truth eta
theta_match = []
phi_match = []
ndf = []
chi2 = []

# LC Relation track 
LC_pt_match = []
LC_track_pt = []
LC_track_eta = []
LC_eta_match = []
LC_track_theta = []
LC_phi_match = []
LC_ndf = []
LC_chi2 = []
LC_d0 = []
LC_z0 = []
LC_nhits = []
LC_pixel_nhits = []
LC_inner_nhits = []
LC_outer_nhits = []
LC_pt_res = [] 
LC_dr = []

# Fake Tracks
fake_pt = []
fake_theta = []
fake_eta = []
fake_phi = []
fake_d0 = []
fake_z0 = []
fake_ndf = []
fake_chi2 = []
fake_nhits = []
fake_pixel_nhits = []
fake_inner_nhits = []
fake_outer_nhits = []

h2d_relpt = [] #pfo muon pt resolution vs pt

with open("bad_res.txt", "w") as textfile:
    pass
with open("no_inner_hits.txt", "w") as textfile:
    pass
no_inner_hits = 0
i = 0
num_matched_tracks = 0
num_dupes = 0
num_fake_tracks = 0
hard_rad_discard = 0
total_n_pfo_mu = 0
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.setReadCollectionNames(["MCParticle", "PandoraPFOs", "SiTracks", "SiTracks_Refitted", "MCParticle_SiTracks", "MCParticle_SiTracks_Refitted", "IBTrackerHits", "IETrackerHits", "OBTrackerHits", "OETrackerHits", "VBTrackerHits", "VETrackerHits"])
# ############## LOOP OVER EVENTS AND FILL HISTOGRAMS  #############################
# Loop over events
for f in fnames:
    if f == "/data/fmeloni/DataMuC_MuColl10_v0A/reco/merged/muonGun_pT_1000_5000.slcio":
        continue
    if max_events > 0 and i >= max_events: break
    #if no_inner_hits > np.abs(max_events): break
    # if f != "/data/fmeloni/DataMuC_MuColl10_v0A/reco/merged/muonGun_pT_250_1000.slcio":
    #     print("Skipped")
    #     continue
    reader.open(f)
    for event in reader: 
        if max_events > 0 and i >= max_events: break
        #if no_inner_hits > np.abs(max_events): break
        if i%100 == 0: print("Processing event %i."%i)

        # Get the collections we care about
        relationCollection = event.getCollection('MCParticle_SiTracks_Refitted')
        relation = pyLCIO.UTIL.LCRelationNavigator(relationCollection)

        mcpCollection = event.getCollection("MCParticle")
        pfoCollection = event.getCollection("PandoraPFOs")
        trackCollection = event.getCollection("SiTracks_Refitted")

        hit_collections = []
        IBTrackerHits = event.getCollection('IBTrackerHits')
        hit_collections.append(IBTrackerHits)
        IETrackerHits = event.getCollection('IETrackerHits')
        hit_collections.append(IETrackerHits)
        OBTrackerHits = event.getCollection('OBTrackerHits')
        hit_collections.append(OBTrackerHits)
        OETrackerHits = event.getCollection('OETrackerHits')
        hit_collections.append(OETrackerHits)
        VBTrackerHits = event.getCollection('VBTrackerHits')
        hit_collections.append(VBTrackerHits)
        VETrackerHits = event.getCollection('VETrackerHits')
        hit_collections.append(VETrackerHits)

        # Make counter variables
        n_mcp_mu = 0
        n_pfo_mu = 0
        has_pfo_mu = False
        my_pfo_mu = 0

        # Pfos 
        ipfo_pt = []
        ipfo_eta = []
        ipfo_phi = []
        ipfo_mu_pt = []
        ipfo_mu_eta = []
        ipfo_mu_phi = []

        # Tracks
        id0_res_vs_pt = []
        id0_res_vs_eta = []
        iz0_res_vs_pt = []
        iz0_res_vs_eta = []
        ipt_res_vs_eta = []
        ipt_res_vs_pt = []
        id0_res = []
        iz0_res = []
        ipt_res = []
        ipt_match = []
        itrack_pt = []
        itrack_eta = []
        ieta_match = []
        itheta_match = []
        iphi_match = []
        indf = []
        ichi2 = []
        id0_res_match = []
        iz0_res_match = []
        inhits = []

        # Fake Tracks
        ifake_pt = []
        ifake_theta = []
        ifake_eta = []
        ifake_phi = []
        ifake_d0 = []
        ifake_z0 = []
        ifake_ndf = []
        ifake_chi2 = []
        ifake_nhits = []

        # Loop over the reconstructed objects and fill histograms
        for pfo in pfoCollection:
            pfo_p = pfo.getMomentum()
            pfo_tlv = ROOT.TLorentzVector()
            pfo_tlv.SetPxPyPzE(pfo_p[0], pfo_p[1], pfo_p[2], pfo.getEnergy())
            hists["pfo_pt"].Fill(pfo_tlv.Perp())
            hists["pfo_eta"].Fill(pfo_tlv.Eta())
            hists["pfo_phi"].Fill(pfo_tlv.Phi())

            ipfo_pt.append(pfo_tlv.Perp())
            ipfo_eta.append(pfo_tlv.Eta())
            ipfo_phi.append(pfo_tlv.Phi())

            if abs(pfo.getType())==13:
                hists["pfo_mu_pt"].Fill(pfo_tlv.Perp())
                hists["pfo_mu_eta"].Fill(pfo_tlv.Eta())
                hists["pfo_mu_phi"].Fill(pfo_tlv.Phi())

                ipfo_mu_pt.append(pfo_tlv.Perp())
                ipfo_mu_eta.append(pfo_tlv.Eta())
                ipfo_mu_phi.append(pfo_tlv.Phi())
                n_pfo_mu += 1
                has_pfo_mu = True
                my_pfo_mu = pfo_tlv     # Storing this to use for matching in the next loop

        # Loop over the truth objects and fill histograms
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
            hists["mcp_pt"].Fill(mcp_tlv.Perp())
            hists["mcp_eta"].Fill(mcp_tlv.Eta())
            hists["mcp_phi"].Fill(mcp_tlv.Phi())

            imcp_pt = []
            imcp_eta = []
            imcp_phi = []

            imcp_pt.append(mcp_tlv.Perp())
            imcp_eta.append(mcp_tlv.Eta())
            imcp_phi.append(mcp_tlv.Phi())

            if abs(mcp.getPDG())==13 and mcp.getGeneratorStatus()==1:
                # Check if the muon radiated significant energy
                hard_rad = check_hard_radiation(mcp, 0.00)
                # if hard_rad:
                #     #print("radiated significant energy, discarding")
                #     hard_rad_discard += 1

                # else:
                hists["mcp_mu_pt"].Fill(mcp_tlv.Perp())
                hists["mcp_mu_eta"].Fill(mcp_tlv.Eta())
                hists["mcp_mu_phi"].Fill(mcp_tlv.Phi())

                imcp_mu_pt = []
                imcp_mu_eta = []
                imcp_mu_phi = []

                imcp_mu_pt.append(mcp_tlv.Perp())
                imcp_mu_eta.append(mcp_tlv.Eta())
                imcp_mu_phi.append(mcp_tlv.Phi())
                n_mcp_mu += 1
                # print("Truth pt, eta, phi:", mcp_tlv.Perp(), mcp_tlv.Eta(), mcp_tlv.Phi())

                if (mcp_tlv.Perp() > 0.5): # Remove ultra-low pt tracks    
                    tracks = relation.getRelatedToObjects(mcp)
                    for track in tracks:
                        Bfield = 5 #T, 3.57 for legacy
                        theta = np.pi/2- np.arctan(track.getTanLambda())
                        phi = track.getPhi()
                        eta = -np.log(np.tan(theta/2))
                        pt  = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
                        track_tlv = ROOT.TLorentzVector()
                        track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
                        dr = mcp_tlv.DeltaR(track_tlv)
                        nhitz = track.getTrackerHits().size()
                        ptres = (mcp_tlv.Perp() - pt) / mcp_tlv.Perp()
                        d0 = track.getD0()
                        z0 = track.getZ0()
                        LC_pt_match.append(imcp_mu_pt)
                        LC_track_pt.append([pt])
                        LC_track_eta.append([eta])
                        LC_eta_match.append(imcp_mu_eta)
                        LC_track_theta.append([theta])
                        LC_phi_match.append([phi])
                        LC_ndf.append([track.getNdf()])
                        LC_chi2.append([track.getChi2()])
                        LC_d0.append([d0])
                        LC_z0.append([z0])
                        LC_nhits.append([nhitz])
                        LC_pt_res.append([ptres])
                        LC_dr.append([dr])

                        LC_pixel_nhit = 0
                        LC_inner_nhit = 0
                        LC_outer_nhit = 0
                        for hit in track.getTrackerHits():
                        # now decode hits
                            encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                            decoder = pyLCIO.UTIL.BitField64(encoding)
                            cellID = int(hit.getCellID0())
                            decoder.setValue(cellID)
                            detector = decoder["system"].value()
                            if detector == 1 or detector == 2:
                                LC_pixel_nhit += 1
                            if detector == 3 or detector == 4:
                                LC_inner_nhit += 1
                            if detector == 5 or detector == 6:
                                LC_outer_nhit += 1
                        LC_pixel_nhits.append([LC_pixel_nhit])
                        LC_inner_nhits.append([LC_inner_nhit])
                        LC_outer_nhits.append([LC_outer_nhit])
                        num_matched_tracks += 1
                        # if LC_outer_nhit == 0:
                        #     no_inner_hits += 1
                        #     with open("no_inner_hits.txt", "a") as textfile:
                        #         textfile.write("Filename: " + str(f) + "\n")
                        #         textfile.write("Event: " + str(event.getEventNumber()) + "\n")
                        #         textfile.write("Truth theta: " + str(mcp_tlv.Theta()) + "\n")
                        
                        # print("Reco pt, eta, phi, nhits, dr:", pt, eta, phi, nhitz, dr)
                        if hard_rad:
                            #print("radiated significant energy, discarding")
                            hard_rad_discard += 1
                            # for d in mcp.getDaughters():
                            #     print(d.getPDG())
                            #print("Nhits (total, pixel, inner, outer):", nhitz, LC_pixel_nhit, LC_inner_nhit, LC_outer_nhit)
                                                
                    # For events in which a PFO mu was reconstructed, fill histograms that will
                    # be used for efficiency. Both numerator and denominator must be filled with truth values!
                    # Also fill resolution histograms
                    if has_pfo_mu:
                        hists["mcp_mu_match_pt"].Fill(mcp_tlv.Perp())
                        hists["mcp_mu_match_eta"].Fill(mcp_tlv.Eta())
                        hists["mcp_mu_match_phi"].Fill(mcp_tlv.Phi())

                        imcp_mu_match_pt = []
                        imcp_mu_match_eta = []
                        imcp_mu_match_phi = []

                        id_mu_dpt = []
                        id_mu_drelpt = []
                        id_mu_deta = []
                        id_mu_dphi = []
                        ih2d_relpt = []

                        imcp_mu_match_pt.append(mcp_tlv.Perp())
                        imcp_mu_match_eta.append(mcp_tlv.Eta())
                        imcp_mu_match_phi.append(mcp_tlv.Phi())

                        hists["d_mu_dpt"].Fill(my_pfo_mu.Perp() - mcp_tlv.Perp())
                        hists["d_mu_drelpt"].Fill((my_pfo_mu.Perp() - mcp_tlv.Perp())/mcp_tlv.Perp())
                        hists["d_mu_deta"].Fill(my_pfo_mu.Eta() - mcp_tlv.Eta())
                        hists["d_mu_dphi"].Fill(my_pfo_mu.Phi() - mcp_tlv.Phi())
                        h_2d_relpt.Fill(mcp_tlv.Perp(), (my_pfo_mu.Perp() - mcp_tlv.Perp())/mcp_tlv.Perp())

                        id_mu_dpt.append(my_pfo_mu.Perp() - mcp_tlv.Perp())
                        id_mu_drelpt.append((my_pfo_mu.Perp() - mcp_tlv.Perp())/mcp_tlv.Perp())
                        id_mu_deta.append(my_pfo_mu.Eta() - mcp_tlv.Eta())
                        id_mu_dphi.append(my_pfo_mu.Phi() - mcp_tlv.Phi())
                        ih2d_relpt.append([mcp_tlv.Perp(), (my_pfo_mu.Perp() - mcp_tlv.Perp())/mcp_tlv.Perp()])

        # Loop over the track objects and fill histograms for D0, Z0, and hit counts
        counter = 0
        max_hits = 0
        best_track = None
        for track in trackCollection:
            Bfield = 5 #T, 3.57 for legacy
            theta = np.pi/2- np.arctan(track.getTanLambda())
            phi = track.getPhi()
            eta = -np.log(np.tan(theta/2))
            pt  = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
            track_tlv = ROOT.TLorentzVector()
            track_tlv.SetPtEtaPhiE(pt, eta, phi, 0)
            dr = mcp_tlv.DeltaR(track_tlv)
            nhitz = track.getTrackerHits().size()
            # print("Reco pt, eta, phi, nhits, dr:", pt, eta, phi, nhitz, dr)

            # Fake tracks
            if len(relation.getRelatedFromObjects(track)) < 1: # If there's no associated truth muon
                ifake_pt.append(pt)
                ifake_theta.append(theta)
                ifake_eta.append(eta)
                ifake_phi.append(phi)
                ifake_d0.append(track.getD0())
                ifake_z0.append(track.getZ0())
                ifake_ndf.append(track.getNdf())
                ifake_chi2.append(track.getChi2())
                ifake_nhits.append(nhitz)
                fake_pixel_nhit = 0
                fake_inner_nhit = 0
                fake_outer_nhit = 0
                for hit in track.getTrackerHits():
                    # now decode hits
                        encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                        decoder = pyLCIO.UTIL.BitField64(encoding)
                        cellID = int(hit.getCellID0())
                        decoder.setValue(cellID)
                        detector = decoder["system"].value()
                        if detector == 1 or detector == 2:
                            fake_pixel_nhit += 1
                        if detector == 3 or detector == 4:
                            fake_inner_nhit += 1
                        if detector == 5 or detector == 6:
                            fake_outer_nhit += 1
                fake_pixel_nhits.append(fake_pixel_nhit)
                fake_inner_nhits.append(fake_inner_nhit)
                fake_outer_nhits.append(fake_outer_nhit)
                num_fake_tracks += 1

            min_dr = 0.005
            if dr < min_dr: # Do dR check first, then do nhits check
                if nhitz > max_hits:
                    max_hits = nhitz
                    best_track = track
                counter += 1
                if counter > 1:
                    num_dupes += 1
                    # print("More than one track in event! # of dupes:", num_dupes)
        # print("# of reco for event (after dr matching):", counter)
        if best_track is not None:
            track = best_track
            d0 = track.getD0()
            z0 = track.getZ0()

            hists["d0_res"].Fill(d0)
            hists["z0_res"].Fill(z0)
            hists["nhits"].Fill(max_hits)

            id0_res.append(d0)
            iz0_res.append(z0)
            inhits.append(max_hits)

            for j, particle_pt in enumerate(imcp_mu_pt):
                particle_eta = imcp_mu_eta[j]
                particle_phi = imcp_mu_phi[j]
                mcp_tlv = ROOT.TLorentzVector()
                mcp_tlv.SetPtEtaPhiE(particle_pt, particle_eta, particle_phi, 0)
                #print(particle_pt, pt)
                ptres = (particle_pt - pt) / particle_pt
                #print(j, dr)
                # Fill 2D histograms
                hists_2d["d0_res_vs_pt"].Fill(particle_pt, d0)
                hists_2d["d0_res_vs_eta"].Fill(particle_eta, d0)
                hists_2d["z0_res_vs_pt"].Fill(particle_pt, z0)
                hists_2d["z0_res_vs_eta"].Fill(particle_eta, z0)
                hists_2d["pt_res_vs_eta"].Fill(particle_eta, ptres)
                hists_2d["pt_res_vs_pt"].Fill(particle_pt, ptres)
                #num_matched_tracks += 1

                id0_res_vs_pt.append([particle_pt, d0])
                id0_res_vs_eta.append([particle_eta, d0])
                iz0_res_vs_pt.append([particle_pt, z0])
                iz0_res_vs_eta.append([particle_eta, z0])
                ipt_res_vs_eta.append([particle_eta, ptres])
                ipt_res_vs_pt.append([particle_pt, ptres])
                ipt_match.append(particle_pt)
                itrack_pt.append(pt)
                itrack_eta.append(eta)
                ieta_match.append(particle_eta)
                iphi_match.append(track.getPhi())
                id0_res_match.append(d0)
                iz0_res_match.append(z0)
                ipt_res.append(ptres)
                #itheta_match.append(mcp_mu_theta[j])
                indf.append(track.getNdf())
                ichi2.append(track.getChi2())
            
                if np.abs(ptres) > 25:
                    with open("bad_res.txt", "a") as textfile:
                        textfile.write("Filename: " + str(f) + "\n")
                        textfile.write("Event: " + str(event.getEventNumber()) + "\n")
                        textfile.write("pt_res: " + str(ptres) + "\n")

                pixel_nhit = 0
                inner_nhit = 0
                outer_nhit = 0
                for hit in track.getTrackerHits():
                    # now decode hits
                        encoding = hit_collections[0].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                        decoder = pyLCIO.UTIL.BitField64(encoding)
                        cellID = int(hit.getCellID0())
                        decoder.setValue(cellID)
                        detector = decoder["system"].value()
                        if detector == 1 or detector == 2:
                            pixel_nhit += 1
                        if detector == 3 or detector == 4:
                            inner_nhit += 1
                        if detector == 5 or detector == 6:
                            outer_nhit += 1
                        
        #print("End of tracks")
        # This is here to check that we never reconstruct multiple muons
        # If we did, we'd have to match the correct muon to the MCP object to do eff/res plots
        # But since we don't, we can skip that step
        if n_pfo_mu > 1: print(n_pfo_mu)
        hists["mcp_n"].Fill(len(mcpCollection))
        hists["pfo_n"].Fill(len(pfoCollection))
        hists["mcp_mu_n"].Fill(n_mcp_mu)
        hists["pfo_mu_n"].Fill(n_pfo_mu)
        hists["mcp_mu_match_n"].Fill(n_pfo_mu)
        i+=1
        if len(id0_res_vs_pt) > 0:
            d0_res.append(id0_res)
            z0_res.append(iz0_res)
            nhits.append(inhits)
            pixel_nhits.append([pixel_nhit])
            inner_nhits.append([inner_nhit])
            outer_nhits.append([outer_nhit])
            #pt_res_hits.append(ipt_res_hits)
            d0_res_vs_pt.append(id0_res_vs_pt)
            d0_res_vs_eta.append(id0_res_vs_eta)
            z0_res_vs_pt.append(iz0_res_vs_pt)
            z0_res_vs_eta.append(iz0_res_vs_eta)
            pt_res_vs_eta.append(ipt_res_vs_eta)
            pt_res_vs_pt.append(ipt_res_vs_pt)
            pt_res.append(ipt_res)
            pt_match.append(ipt_match)
            track_pt.append(itrack_pt)
            track_eta.append(itrack_eta)
            eta_match.append(ieta_match)
            theta_match.append(itheta_match)
            phi_match.append(iphi_match)
            ndf.append(indf)
            chi2.append(ichi2)
            d0_res_match.append(id0_res_match)
            z0_res_match.append(iz0_res_match)
        mcp_pt.append(imcp_pt)
        mcp_eta.append(imcp_eta)
        mcp_phi.append(imcp_phi)
        if hard_rad == False:
            mcp_mu_pt.append(imcp_mu_pt)
            mcp_mu_eta.append(imcp_mu_eta)
            mcp_mu_phi.append(imcp_mu_phi)
            if has_pfo_mu:    
                mcp_mu_match_pt.append(imcp_mu_match_pt)
                mcp_mu_match_eta.append(imcp_mu_match_eta)
                mcp_mu_match_phi.append(imcp_mu_match_phi)
                d_mu_dpt.append(id_mu_dpt)
                d_mu_drelpt.append(id_mu_drelpt)
                d_mu_deta.append(id_mu_deta)
                d_mu_dphi.append(id_mu_dphi)
                h2d_relpt.append(ih2d_relpt)
        pfo_pt.append(ipfo_pt)
        pfo_eta.append(ipfo_eta)
        pfo_phi.append(ipfo_phi)
        pfo_mu_pt.append(ipfo_mu_pt)
        pfo_mu_eta.append(ipfo_mu_eta)
        pfo_mu_phi.append(ipfo_mu_phi)
        if len(ifake_pt) > 0:
            fake_pt.append(ifake_pt)
            fake_theta.append(ifake_theta)
            fake_eta.append(ifake_eta)
            fake_phi.append(ifake_phi)
            fake_d0.append(ifake_d0)
            fake_z0.append(ifake_z0)
            fake_ndf.append(ifake_ndf)
            fake_chi2.append(ifake_chi2)
            fake_nhits.append(ifake_nhits)
        #     # print(fake_pt)
    reader.close()

# ############## MANIPULATE, PRETTIFY, AND SAVE HISTOGRAMS #############################
print("\nSummary statistics:")
print("Ran over %i events."%i)
print("Found:")
print("\t%i MCPs"%hists["mcp_pt"].GetEntries())
print("\t%i mu MCPs"%hists["mcp_mu_pt"].GetEntries())
# print("\tSanity check mcp_mu_pt:", len(mcp_mu_pt))
# print("\t%i PFOs"%hists["pfo_pt"].GetEntries())
# print("\t%i mu PFOs"%hists["pfo_mu_pt"].GetEntries())
print('\t%i matched muon tracks'%(num_matched_tracks))
print('\t%i duplicates eliminated'%num_dupes)
print('\t%i hard radiations discarded'%hard_rad_discard)
print('\t%i fake tracks'%num_fake_tracks)
# print('\t%i GeV'%np.max(mcp_mu_pt))


# Make a list of all the data you want to save
data_list = {}
data_list["mcp_pt"] = mcp_pt
data_list["mcp_eta"] = mcp_eta
data_list["mcp_phi"] = mcp_phi
data_list["mcp_mu_pt"] = mcp_mu_pt
data_list["mcp_mu_eta"] = mcp_mu_eta
data_list["mcp_mu_phi"] = mcp_mu_phi

# data_list["pfo_pt"] = pfo_pt
# data_list["pfo_eta"] = pfo_eta
# data_list["pfo_phi"] = pfo_phi
# data_list["pfo_mu_pt"] = pfo_mu_pt
# data_list["pfo_mu_eta"] = pfo_mu_eta
# data_list["pfo_mu_phi"] = pfo_mu_phi

# data_list["mcp_mu_match_pt"] = mcp_mu_match_pt
# data_list["mcp_mu_match_eta"] = mcp_mu_match_eta
# data_list["mcp_mu_match_phi"] = mcp_mu_match_phi
# data_list["d_mu_dpt"] = d_mu_dpt
# data_list["d_mu_drelpt"] = d_mu_drelpt
# data_list["d_mu_dphi"] = d_mu_dphi
# data_list["d_mu_deta"] = d_mu_deta

data_list["d0_res"] = d0_res
data_list["z0_res"] = z0_res
data_list["nhits"] = nhits
data_list["pixel_nhits"] = pixel_nhits
data_list["inner_nhits"] = inner_nhits
data_list["outer_nhits"] = outer_nhits
data_list["pt_res_hits"] = pt_res_hits
data_list["d0_res_vs_pt"] = d0_res_vs_pt
data_list["d0_res_vs_eta"] = d0_res_vs_eta
data_list["z0_res_vs_pt"] = z0_res_vs_pt
data_list["z0_res_vs_eta"] = z0_res_vs_eta
data_list["pt_res_vs_eta"] = pt_res_vs_eta
data_list["pt_res_vs_pt"] = pt_res_vs_pt
data_list["pt_res"] = pt_res
data_list["pt_match"] = pt_match
data_list["track_pt"] = track_pt
data_list["track_eta"] = track_eta
data_list["eta_match"] = eta_match
data_list["theta_match"] = theta_match
data_list["phi_match"] = phi_match
data_list["ndf"] = ndf
data_list["chi2"] = chi2
data_list["d0_res_match"] = d0_res_match
data_list["z0_res_match"] = z0_res_match

data_list["LC_pt_match"] = LC_pt_match
data_list["LC_track_pt"] = LC_track_pt
data_list["LC_track_eta"] = LC_track_eta
data_list["LC_eta_match"] = LC_eta_match
data_list["LC_track_theta"] = LC_track_theta
data_list["LC_phi_match"] = LC_phi_match
data_list["LC_ndf"] = LC_ndf
data_list["LC_chi2"] = LC_chi2
data_list["LC_d0"] = LC_d0
data_list["LC_z0"] = LC_z0
data_list["LC_nhits"] = LC_nhits
data_list["LC_pixel_nhits"] = LC_pixel_nhits
data_list["LC_inner_nhits"] = LC_inner_nhits
data_list["LC_outer_nhits"] = LC_outer_nhits
data_list["LC_pt_res"] = LC_pt_res
data_list["LC_dr"] = LC_dr

data_list["fake_pt"] = fake_pt
data_list["fake_theta"] = fake_theta
data_list["fake_eta"] = fake_eta
data_list["fake_phi"] = fake_phi
data_list["fake_d0"] = fake_d0
data_list["fake_z0"] = fake_z0
data_list["fake_ndf"] = fake_ndf
data_list["fake_chi2"] = fake_chi2
data_list["fake_nhits"] = fake_nhits
data_list["fake_pixel_nhits"] = fake_pixel_nhits
data_list["fake_inner_nhits"] = fake_inner_nhits
data_list["fake_outer_nhits"] = fake_outer_nhits

data_list["h_2d_relpt"] = h2d_relpt

# After the loop is finished, save the data_list to a .json file
output_json = "v0_noBIB_all.json"
with open(output_json, 'w') as fp:
    json.dump(data_list, fp)

# Draw all the 1D histograms you filled
for i, h in enumerate(hists):
    c = ROOT.TCanvas("c%i"%i, "c%i"%i)
    hists[h].Draw()
    hists[h].GetXaxis().SetTitle(h)
    hists[h].GetYaxis().SetTitle("Entries")

    # For resolution plots, fit them and get the mean and sigma
    if h.startswith("d_mu"):
        f = ROOT.TF1("f%i"%i, "gaus")
        f.SetLineColor(ROOT.kRed)
        hists[h].Fit("f%i"%i)
        c.SetLogy()
        latex = ROOT.TLatex()
        p = f.GetParameters()
        latex.DrawLatexNDC(.64, .85, "Mean: %f"%p[1])
        latex.DrawLatexNDC(.64, .78, "Sigma: %f"%p[2])
    c.SaveAs("plots/%s.png"%h)

# Make efficiency plots
# In these files, there are at most 1 PFO mu, so matching isn't needed
for v in variables:
    if v=="n": continue
    c = ROOT.TCanvas("c%s"%v, "c%s"%v)
    eff = ROOT.TEfficiency(hists["mcp_mu_match_"+v], hists["mcp_mu_"+v])
    eff.Draw("ape")
    ROOT.gPad.Update()
    eff.SetLineWidth(2)
    eff.GetPaintedGraph().SetMinimum(0)
    eff.GetPaintedGraph().SetMaximum(1)
    eff.SetTitle(";%s;Efficiency"%v)
    c.SaveAs("plots/eff_%s.png"%v)

# Make 2D plot and a TProfile to understand pT resolution v pT
c = ROOT.TCanvas("crelpt2d", "crelpt2d")
h_2d_relpt.Draw("colz")
h_2d_relpt.GetXaxis().SetTitle("pt")
h_2d_relpt.GetYaxis().SetTitle("drelpt")
c.SaveAs("plots/d_mu_relpt_2d.png")

c = ROOT.TCanvas("crelpt2dprof", "crelpt2dprof")
h_prof = h_2d_relpt.ProfileX("_pfx", 1, -1, "s")
h_prof.GetXaxis().SetTitle("pt")
h_prof.GetYaxis().SetTitle("drelpt")
h_prof.Draw()
c.SaveAs("plots/d_mu_relpt_prof.png")

for var in hists_2d:
    c = ROOT.TCanvas("c_" + var, "c_" + var)
    hists_2d[var].Draw("colz")
    hists_2d[var].GetXaxis().SetTitle(variables_2d[var]["xlabel"])
    hists_2d[var].GetYaxis().SetTitle(variables_2d[var]["ylabel"])
    c.SaveAs("plots/" + var + ".png")



