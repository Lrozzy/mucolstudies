from pyLCIO import IOIMPL, EVENT, UTIL
import ROOT
from math import*
import numpy as np
#from drawtrack import getCurvature, getTransverseImpact, getX, getY
import drawtrack
ROOT.gStyle.SetOptStat(0)

# Adding 2D histograms
hists_2d = {}
variables_2d = {
    "z_vs_r": {"xbins": 100, "xmin": -1000, "xmax": 1000, "ybins": 100, "ymin": 0, "ymax": 1000,
                    "ylabel": "R (mm)", "xlabel": "Z (mm)"},
    "x_vs_y": {"xbins": 100, "xmin": -500, "xmax": 500, "ybins": 100, "ymin": -500, "ymax": 500,
                    "ylabel": "Y (mm)", "xlabel": "X (mm)"},
}

# Read the .txt file and extract filenames and event numbers
def read_txt_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    fnames = []
    event_numbers = []
    for idx, line in enumerate(lines):
        if line.startswith("Filename:"):
            fnames.append(line.split(":")[1].strip())
        elif line.startswith("Event:"):
            event_numbers.append(int(line.split(":")[1].strip()))
    return fnames, event_numbers

Bfield = 5 # T
filename = 'bad_res.txt'
fnames, event_numbers = read_txt_file(filename)
for f, event_number in zip(fnames, event_numbers):
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open([f])
    for event in reader:
        if event.getEventNumber() != event_number:
            continue
        counter = 0
        theta_l = []
        rho_l = []
        pt_l = []
        phi_l = []
        s_l = []
        w_l = []
        coslambda_l = []
        truth_vt_l = []
        truth_phi_l = []
        truth_pt_l = []
        truth_mass_l = []
        detector_l = []
        layer_l = []
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
        for var in variables_2d:
            hists_2d[var] = ROOT.TH2F(var+"_" +str(event_number), var, variables_2d[var]["xbins"], variables_2d[var]["xmin"],
                              variables_2d[var]["xmax"], variables_2d[var]["ybins"], variables_2d[var]["ymin"],
                              variables_2d[var]["ymax"])
        trackCollection = event.getCollection("SiTracks")
        mcpCollection = event.getCollection("MCParticle")
        for mcp in mcpCollection:
            mcp_p = mcp.getMomentum()
            #mcp_b = mcp.beta()
            mcp_tlv = ROOT.TLorentzVector()
            mcp_tlv.SetPxPyPzE(mcp_p[0], mcp_p[1], mcp_p[2], mcp.getEnergy())
            truth_pt = mcp_tlv.Perp() #pT in GeV
            if abs(mcp.getPDG())==13 and mcp.getGeneratorStatus()==1:
                truth_phi = mcp_tlv.Phi()
                truth_eta = mcp_tlv.Eta()
                truth_mass = mcp_tlv.M() #rest mass in GeV
                truth_vt = truth_pt/truth_mass
                truth_phi_l.append(truth_phi)
                truth_pt_l.append(truth_pt)
                truth_vt_l.append(truth_vt)
                truth_mass_l.append(truth_mass)
        for track in trackCollection:
            d0 = track.getD0()
            z0 = track.getZ0()
            phi = track.getPhi()
            theta = np.pi/2- np.arctan(track.getTanLambda())
            theta_l.append(theta)
            eta = -np.log(np.tan(theta/2))
            w = track.getOmega()
            pt = 0.3 * Bfield / fabs(w* 1000.) #pT in GeV
            # rho = (1/0.33) * pt / (Bfield) # THIS IS R, NOT RHO; alternative formula rho = 1/fabs(track.getOmega())
            rho = 0.33*pt/Bfield 
            s = (0.3*Bfield*1.5**2)/(8*pt) #sagitta?
            coslambda = np.cos(np.arctan(track.getTanLambda()))
            rho_l.append(rho)   
            pt_l.append(pt) 
            phi_l.append(phi)
            s_l.append(s)
            w_l.append(w)
            coslambda_l.append(coslambda)
            for hit in track.getTrackerHits():
                # now decode hits
                encoding = hit_collections[0].getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID)
                detector = decoder["system"].value()
                layer = decoder["layer"].value()
                if detector == 1:
                    detector = "VB"
                elif detector == 2:
                    detector = "VE"
                elif detector == 3:
                    detector = "IB"
                elif detector == 4:
                    detector = "IE"
                elif detector == 5:
                    detector = "OB"
                elif detector == 6:
                    detector = "OE"
                detector_l.append(detector)
                layer_l.append(layer)
                #print(detector, layer)
                x = hit.getPosition()[0]
                y = hit.getPosition()[1]
                z = hit.getPosition()[2]
                # Fill the histograms with the particle positions
                hists_2d["z_vs_r"].Fill(z, np.sqrt(x**2 + y**2))
                hists_2d["x_vs_y"].Fill(x, y)
        for var in hists_2d:    
            c = ROOT.TCanvas("c_" + var + "_" + str(event_number), "c_" + var) 
            hists_2d[var].Draw()
            hists_2d[var].SetMarkerStyle(ROOT.kFullCircle)
            funcs = []
            funcs1 = []
            for i, (rho, theta, pt, phi, s, coslambda, truth_phi, truth_pt, truth_vt, truth_mass, w) in enumerate(zip(rho_l, theta_l,pt_l, phi_l, s_l, coslambda_l, truth_phi_l, truth_pt_l, truth_vt_l, truth_mass_l, w_l)):
                r_func = '2*{rho}*sin(x*tan({theta})/(2*{rho}))'.format(rho=rho, theta=theta)
                if var == "z_vs_r":
                    legend = ROOT.TLegend(0.1,0.5,0.6,0.9)
                    tf1 = ROOT.TF1("truth_track_zr_" + str(i), r_func, -1000, 500)
                    tf1.Draw('same')
                    tf1.SetLineColor(ROOT.kRed)
                    funcs.append(tf1)
                elif var == "x_vs_y":
                    legend = ROOT.TLegend(0.1,0.7,0.6,0.9)
                    r = 1000.*truth_pt/(0.3*Bfield)
                    truth_omega = Bfield/truth_mass
                    w = w
                    numPoints = 10000
                    sMin = -1.125
                    sMax = 1.125
                    tMin = -0.05
                    tMax = 0.05
                    helix = ROOT.TGraph()
                    for j in range(numPoints):
                        # s = (sMin + (sMax - sMin) * j / (numPoints - 1))
                        # x = r * (ROOT.TMath.Cos(truth_phi + s*coslambda/r)-ROOT.TMath.Cos(truth_phi))
                        # y = r * (ROOT.TMath.Sin(truth_phi + s*coslambda/r)-ROOT.TMath.Sin(truth_phi))
                        # print(r,s,coslambda,truth_phi)
                        t = (tMin + (tMax - tMin) * j / (numPoints - 1))
                        x = (truth_vt/w)*np.sin(w*t)
                        y = (truth_vt/w)*np.cos(w*t)
                        helix.SetPoint(j, x, y)
                    print(truth_omega, w)
                    helix.Draw("same")
                    helix.SetLineColor(ROOT.kBlue)
                    funcs1.append(helix)  
            legend.SetBorderSize(0)
            legend.SetFillColor(0)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.04)
            legend.SetHeader
            legend.SetNColumns(2)
            legend.AddEntry(0, "pT: {:.2e} GeV".format(pt), "")
            legend.AddEntry(0, "Truth pT: {:.2e}".format(truth_pt), "")
            legend.AddEntry(hists_2d[var], "Tracker Hits", "P")     
            # legend.AddEntry(0, "", "")
            # legend.AddEntry(0, "eta: {}".format(eta), "")
            # legend.AddEntry(0, "", "")
            # legend.AddEntry(0, "Truth eta: {}".format(truth_eta), "")
            # legend.AddEntry(0, "", "")
            # legend.AddEntry(0, "phi: {}".format(phi), "")
            # legend.AddEntry(0, "", "")
            # legend.AddEntry(0, "Truth phi: {}".format(truth_phi), "")
            if var == "z_vs_r":
                legend.AddEntry(funcs[0], "Truth Tracks", "L")
                for x,y in zip(detector_l, layer_l):
                    legend.AddEntry(0, "Detector: {}".format(x), "")
                    legend.AddEntry(0, "Layer: {}".format(y), "")
            elif var == "x_vs_y":
                legend.AddEntry(funcs1[0], "Truth Tracks", "L") 
            legend.Draw()     
            hists_2d[var].GetXaxis().SetTitle(variables_2d[var]["xlabel"])
            hists_2d[var].GetYaxis().SetTitle(variables_2d[var]["ylabel"])
            c.SaveAs("bad_res_plots/" + var + "_" + str(event_number) + ".pdf")
    reader.close()



