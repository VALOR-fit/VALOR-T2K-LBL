#!/usr/bin/env python

# A bunch of utilities used for testing the T2K run 10 raw mtuples and processed DSTs

import ROOT

#________________________________________________________________________________________
# Compare debug trees made when processing VALOR MC DSTs against the same from p-theta's
def compare_trees(file_v, file_p, attributes = ["enu_true", "enu_reco"], maxn = 5, precision_we_care_about = 1e-6, only_print_mismatches=False, print_list_of_erec_differences=False, just_print_summary=False):
    attribute_name_map = {"enu_true":"enu_true", "enu_reco":"enu_reco", "Theta_reco":"Theta_reco", "wght":"ptheta_weight_sans_mcscale", "neut_mode":"neut_mode", "valor_reac_mode":"ptheta_int_cat"}
    ptheta_mode_to_valor_mode = {0:"CCQE", 1:"CCMEC", 2:"CC1pi", 3:"CCcoherent", 4:"CCDIS", 5:"CCMpi", 6:"CCmisc", 7:"NC1pi0", 8:"NC1piC", 9:"NCcoherent", 10:"NC1Gamma", 11:"NCother"}
    ptheta_channel_to_valor_channel = {"numu":"NuMu to NuMu", "numubar":"NuMuBar to NuMuBar", "nue":"NuE to NuE", "nuebar":"NuEBar to NuEBar", "nue_sig":"NuMu to NuE", "nuebar_sig":"NuMuBar to NuEBar"}
    nev_mismatch = 0
    attr_mismatch = dict()
    for attr in attributes:
        attr_mismatch[attr] = 0
    list_diff_erec = []
    tree_valor = file_v.Get("mode_tree")
    tree_ptheta = file_p.Get("outtree")
    if tree_valor.GetEntries() != tree_ptheta.GetEntries():
        print "ERROR: valor tree has", tree_valor.GetEntries(), "events, while the p-theta tree has", tree_ptheta.GetEntries(), ". Will compare what we can."
        print "  valor filename:", file_v.GetName(), ",  ptheta filename:", file_p.GetName()
        return
        maxn_temp = min(tree_valor.GetEntries(), tree_ptheta.GetEntries())
        if maxn > 0:
            maxn_temp = min(maxn, maxn_temp)
        maxn = maxn_temp
    if maxn < 1:
        maxn = tree_valor.GetEntries()
    for i in range(0, maxn):
        if i > maxn: break
        tree_valor.GetEntry(i)
        tree_ptheta.GetEntry(i)
        found_mismatch = False
        print_str = str(i) + " \t\t\t valor \t\t\t\t p-theta \n"
        for attr in attributes:
            val_v = getattr(tree_valor, attr)
            val_p = getattr(tree_ptheta, attribute_name_map[attr])
            spacing = " \t\t\t"
            name_spacing = " \t"
            attr_print = attr
            if attr == "valor_reac_mode":
                spacing = " \t"
                fname_p = file_p.GetName().rsplit("/")[-1].replace(".root","")  # remove the file path and ".root" from the filename
                channel_p = ptheta_channel_to_valor_channel[fname_p]
                mode_p = ptheta_mode_to_valor_mode[val_p]
                val_p = channel_p +  " "  + mode_p
                channel_v = util.FlavourAsString(util.ReacModeToInitFlavor(val_v)) + " to " + util.FlavourAsString(util.ReacModeToFinalFlavor(val_v))
                mode_v = util.ReacTypeAsString(util.ReacModeToType(val_v))
                if mode_v == "CC1piC" or mode_v == "CC1pi0":
                    mode_v = "CC1pi"
                val_v = channel_v + " " + mode_v
                attr_print = "reac_mode" # just want a shorter print name, like the other attributes
            elif attr == "wght":
                name_spacing = " \t\t"
            elif attr == "enu_reco":
                if ROOT.TMath.Abs(val_v - val_p) > precision_we_care_about:
                    list_diff_erec.append( ROOT.TMath.Abs(val_v - val_p))
            if type(val_p) == type(getattr(tree_valor, "neut_mode")):
                spacing = " \t\t\t\t"
            print_str += " \t" + attr_print + name_spacing + str(val_v) + spacing + str(val_p)
            if type(val_v) != type(getattr(tree_valor, "enu_true")):
                if val_v != val_p:
                    if not found_mismatch:
                        nev_mismatch += 1
                    found_mismatch = True
                    print_str += "\t\t\t --> !!MISMATCH!!"
                    attr_mismatch[attr] += 1
            else:
                if ROOT.TMath.Abs(val_v - val_p) > precision_we_care_about:
                    if not found_mismatch:
                        nev_mismatch += 1
                    found_mismatch = True
                    print_str += "\t\t\t --> !!MISMATCH!!"
                    attr_mismatch[attr] += 1
            print_str += "\n"
        print_str += "\n"
        if only_print_mismatches:
            if found_mismatch and not just_print_summary:
                print print_str
        else:
            print print_str
    if not just_print_summary:
        print "\n"
    print "After comparing", maxn, "events,", nev_mismatch, "of which contained mismatches. Broken down by attribute, this is:", attr_mismatch
    if print_list_of_erec_differences:
        list_diff_erec.sort(reverse=True)
        print list_diff_erec

#________________________________________________________________________________________
# Run compare_trees() on the VALOR and p-theta tree pairs for all beams, samples and osc channels
def compare_all_trees():
    import sys
    valor_beams = ["JPARC_FHC", "JPARC_RHC"]
    ptheta_beams = ["numode", "antinu"]
    valor_samples = ["1rmu", "1re", "multire"]
    ptheta_samples = ["numu", "nue", "cc1pi"]
    valor_channels = ["NuMu_no_osc", "NuMuBar_no_osc", "NuE_no_osc", "NuE_osc", "NuEBar_no_osc", "NuEBar_osc"]
    ptheta_channels = ["numu", "numubar", "nue", "nue_sig", "nuebar", "nuebar_sig"]
    for ib in range(0, len(valor_beams)):
        vb = valor_beams[ib]
        pb = ptheta_beams[ib]
        for isa in range(0, len(valor_samples)):
            vs = valor_samples[isa]
            ps = ptheta_samples[isa]
            if vb == "JPARC_RHC" and vs == "multire":
                # We don't use this beam + sample combination in the analysis
                continue
            for ic in range(0, len(valor_channels)):
                vc = valor_channels[ic]
                pc = ptheta_channels[ic]
                vfn = "valor_trees/"+vs+"_"+vb+"_"+vc+".root"               # TODO: WARNING: hardcoded paths!
                pfn = "p-theta_trees/"+pb+"_"+ps+"_binning/"+pc+".root"
                vf = ROOT.TFile(vfn, "read")
                pf = ROOT.TFile(pfn, "read")
                if not vf or not vf.IsOpen():
                    print "ERROR: file", vfn, "not found!"
                    sys.exit()
                if not pf or not pf.IsOpen():
                    print "ERROR: file", pfn, "not found!"
                    sys.exit()
                print "Comparing:", vb, vs, vc
                compare_trees(vf, pf, ["enu_true","enu_reco","Theta_reco","wght","neut_mode","valor_reac_mode"], -1, 6e-5, True)
                print "\n\n\n"
                vf.Close()
                pf.Close()

#________________________________________________________________________________________
# Compare the number of events in the VALOR and p-theta debug trees to the number of events in the raw mtuples
def compare_nevents():
    valor_beams = ["JPARC_FHC", "JPARC_RHC"]
    ptheta_beams = ["numode", "antinu"]
    raw_beams = ["fhc", "rhc"]
    valor_samples = ["1rmu", "1re", "multire"]
    ptheta_samples = ["numu", "nue", "cc1pi"]
    raw_samples = ["numu", "nue", "cc1pi"]
    valor_channels = ["NuMu_no_osc", "NuMuBar_no_osc", "NuE_no_osc", "NuE_osc", "NuEBar_no_osc", "NuEBar_osc"]
    ptheta_channels = ["numu", "numubar", "nue", "nue_sig", "nuebar", "nuebar_sig"]
    raw_channels = ["numu_x_numu", "numubar_x_numubar", "nue_x_nue", "numu_x_nue", "nuebar_x_nuebar", "numubar_x_nuebar"]
    print "\t\t\t\t mtuple \t valor \t p-theta"
    for ib in range(0, len(valor_beams)):
        vb = valor_beams[ib]
        pb = ptheta_beams[ib]
        rb = raw_beams[ib]
        for isa in range(0, len(valor_samples)):
            vs = valor_samples[isa]
            ps = ptheta_samples[isa]
            rs = raw_samples[isa]
            if vb == "JPARC_RHC" and vs == "multire":
                # We don't use this beam + sample combination in the analysis
                continue
            for ic in range(0, len(valor_channels)):
                vc = valor_channels[ic]
                pc = ptheta_channels[ic]
                rc = raw_channels[ic]
                vfn = "valor_trees/"+vs+"_"+vb+"_"+vc+".root"               # TODO: WARNING: hardcoded paths!
                pfn = "p-theta_trees/"+pb+"_"+ps+"_binning/"+pc+".root"
                rfn = "../t2ksk19b.fqv4r0."+rb+"."+rc+"_"+rs+"selec.root"
                vf = ROOT.TFile(vfn, "read")
                pf = ROOT.TFile(pfn, "read")
                rf = ROOT.TFile(rfn, "read")
                if not vf or not vf.IsOpen():
                    print "ERROR: file", vfn, "not found!"
                    sys.exit()
                if not pf or not pf.IsOpen():
                    print "ERROR: file", pfn, "not found!"
                    sys.exit()
                if not rf or not rf.IsOpen():
                    print "ERROR: file", rfn, "not found!"
                    sys.exit()
                npadding = len("cc1pi numubar_x_numubar") - len(rs+rc) +1
                padding = " "*npadding
                print "For", rb, rs, rc, padding, rf.Get("mtuple").GetEntries(), ",\t", vf.Get("mode_tree").GetEntries(), ",\t", pf.Get("outtree").GetEntries()
                vf.Close()
                pf.Close()
                rf.Close()