#!/usr/bin/env python2

# Run plot_contours_TN.py many times to produce all of the constant dchisq plots needed for the T2K run 10 analysis.
# Also constructs LaTex slides from the produced figures and confidence regions.
#
# BEWARE: There are many hardcoded things in this macro. Most of the hardcoding is done in a block just below
#         this comment, but there is also hardcoding throughout this macro, with e.g. file and directory names
#         being assumed to follow a certain structure etc.
#
# Francis Bench 2020-11-11

import argparse
import os
from datetime import datetime

# TODO: Maybe put chisq calc w.r.t. MO on the plots, or is okay just in the slide title?

#________________________________________________________________________________________________________
# HARDCODED STUFF HERE!
fit_params    = ["cpv", "mixing13", "mixing23", "dm23sq", "mixing13_cpv", "mixing23_cpv", "mixing23_dm23sq"]
fit_datasets  = ["AsimovA", "AsimovB", "data"]
#fit_datasets  = ["AsimovA", "AsimovB"]
fit_reactors  = ["reactor", "t2k"]
fit_orderings = ["normal", "inverted"]

fit_params_ptheta = {
    "cpv"             : "dCP",
    "mixing13"        : None,          # Don't have p-theta's Erec-theta files for this
    "mixing23"        : None,
    "dm23sq"          : None,
    "mixing13_cpv"    : "th13_dCP",
    "mixing23_cpv"    : None,
    "mixing23_dm23sq" : "dm2_32_th23"
}
reac_ptheta = {
    "reactor" : "",
    "t2k"     : "_noRC"
}

fit_params_mach3 = {
    "cpv"             : "1D_dcp",
    "mixing13"        : "1D_th13",
    "mixing23"        : "1D_th23",
    "dm23sq"          : "1D_dm23",
    "mixing13_cpv"    : "th13dcp",
    "mixing23_cpv"    : None,          # Don't usually make comparisons for this plot, as we already do 2D comparisons with mixing23 and cpv separately in other plots
    "mixing23_dm23sq" : "th23dm23"
}
reac_mach3 = {
    "reactor" : "_wRC",
    "t2k"     : ""
}
mo_mach3 = {
    "normal"   : "NH",
    "inverted" : "IH"
}

colours = [1, 856, 629]

script = "python " +str(os.environ['VALOR_EXPT'])+"/macros/contour/plot_contours_TN.py"
standard_args    = " --line-width 2 --text-size 0.05 --title-size 0.055"

reac_out_name = {
    "reactor" : "wRC",
    "t2k"     : "woRC"
}
mo_out_name = {
    "normal"   : "NO",
    "inverted" : "IO"
}

is_2d = {
    "cpv"             : False,
    "mixing13"        : False,
    "mixing23"        : False,
    "dm23sq"          : False,
    "mixing13_cpv"    : True,
    "mixing23_cpv"    : True,
    "mixing23_dm23sq" : True
}

title_args = {
    "cpv"             : " --title-x \"#delta_{CP} (Radians)\" --title-y \"#Delta #chi^{2}\"",
    "mixing13"        : " --title-x \"sin^{2} #theta_{13}\" --title-y \"#Delta #chi^{2}\"",
    "mixing23"        : " --title-x \"sin^{2} #theta_{23}\" --title-y \"#Delta #chi^{2}\"",
    "dm23sq"          : {
        "both"     : " --title-x \"#Delta m^{2}_{23} (NO) / |#Delta m^{2}_{31}| (IO) (eV^{2} / c^{4})\" --title-y \"#Delta #chi^{2}\"",
        "normal"   : " --title-x \"#Delta m^{2}_{23} (eV^{2} / c^{4})\" --title-y \"#Delta #chi^{2}\"",
        "inverted" : " --title-x \"|#Delta m^{2}_{31}| (eV^{2} / c^{4})\" --title-y \"#Delta #chi^{2}\""
    },
    "mixing13_cpv"    : " --title-x \"sin^{2} #theta_{13}\" --title-y \"#delta_{CP} (Radians)\"",
    "mixing23_cpv"    : " --title-x \"sin^{2} #theta_{23}\" --title-y \"#delta_{CP} (Radians)\"",
    "mixing23_dm23sq" : {
        "both"     : " --title-x \"sin^{2} #theta_{23}\" --title-y \"#Delta m^{2}_{23} (NO) / |#Delta m^{2}_{31}| (IO) (eV^{2} / c^{4})\"",
        "normal"   : " --title-x \"sin^{2} #theta_{23}\" --title-y \"#Delta m^{2}_{23} (eV^{2} / c^{4})\"",
        "inverted" : " --title-x \"sin^{2} #theta_{23}\" --title-y \"|#Delta m^{2}_{31}| (eV^{2} / c^{4})\""
    }
}

range_args = {
    "cpv"             : None,    # None = use full range
    "mixing13"        : None,
    "mixing23"        : None,
    "dm23sq"          : None,
    "mixing13_cpv"    : {
        "reactor" : " --x-range 19e-3 25e-3",
        "t2k"     : " --x-range 7e-3 50e-3"
    },
    "mixing23_cpv"    : None,
    "mixing23_dm23sq" : " --x-range 0.4 0.65 --y-range 2.25e-3 2.75e-3"
}

# Look on my Work, ye Mighty, and despair, for it is long and poorly organised...
leg_coords_args = {
    "cpv"             : {
        "AsimovA" : {
            "reactor" : {
                "valor"   : " --leg2-coords 0.15 0.80 0.45 0.90",
                "group"   : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.15 0.60 0.45 0.70"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.15 0.80 0.45 0.90",
                "group" : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.15 0.60 0.45 0.70"
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.65 0.45 0.75"
        },
        "AsimovB" : {
            "reactor" : {
                "valor" : " --leg2-coords 0.15 0.80 0.45 0.90",
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.15 0.80 0.45 0.90",
            },
            "both_RC" : " --leg-coords 0.50 0.83 0.80 0.93 --leg2-coords 0.15 0.83 0.45 0.93"
        },
        "data"    : {
            "reactor" : {
                "valor" : " --leg2-coords 0.15 0.60 0.45 0.70",
                "group" : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.15 0.60 0.45 0.70"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.15 0.60 0.45 0.70",
                "group" : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.15 0.60 0.45 0.70"
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.65 0.45 0.75"
        }
    },
    "mixing13"        : {
        "AsimovA" : {
            "t2k"     : {
                "valor" : " --leg-coords  0.45 0.85 0.75 0.90 --leg2-coords 0.45 0.70 0.75 0.80",
                "group" : " --leg-coords  0.45 0.75 0.75 0.90 --leg2-coords 0.45 0.60 0.75 0.70"
            }
        },
        "AsimovB" : {
            "t2k"     : {
                "valor" : " --leg-coords  0.45 0.85 0.75 0.90 --leg2-coords 0.45 0.70 0.75 0.80",
                "group" : " --leg-coords  0.45 0.75 0.75 0.90 --leg2-coords 0.45 0.60 0.75 0.70"
            }
        },
        "data"    : {
            "t2k"     : {
                "valor" : " --leg-coords  0.45 0.80 0.75 0.90 --leg2-coords 0.45 0.65 0.75 0.75",
                "group" : " --leg-coords  0.45 0.80 0.75 0.90 --leg2-coords 0.45 0.65 0.75 0.75"
            }
        }
    },
    "mixing23"        : {
        "AsimovA" : {
            "reactor" : {
                "valor" : " --leg2-coords 0.45 0.80 0.75 0.90",
                "group" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.45 0.80 0.75 0.90",
                "group" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
            },
            "both_RC" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
        },
        "AsimovB" : {
            "reactor" : {
                "valor" : " --leg2-coords 0.45 0.80 0.75 0.90",
                "group" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.45 0.80 0.75 0.90",
                "group" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
            },
            "both_RC" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
        },
        "data"    : {
            "reactor" : {
                "valor" : " --leg2-coords 0.45 0.80 0.75 0.90",
                "group" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.45 0.80 0.75 0.90",
                "group" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
            },
            "both_RC" : " --leg-coords  0.45 0.80 0.75 0.90  --leg2-coords 0.45 0.65 0.75 0.75"
        }
    },
    "dm23sq"          : {
        "AsimovA" : {
            "reactor" : {
                "valor" : " --leg2-coords 0.40 0.80 0.70 0.90",
                "group" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.40 0.80 0.70 0.90",
                "group" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
            },
            "both_RC" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
        },
        "AsimovB" : {
            "reactor" : {
                "valor" : " --leg2-coords 0.40 0.80 0.70 0.90",
                "group" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.40 0.80 0.70 0.90",
                "group" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
            },
            "both_RC" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
        },
        "data"    : {
            "reactor" : {
                "valor" : " --leg2-coords 0.40 0.80 0.70 0.90",
                "group" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
            },
            "t2k"     : {
                "valor" : " --leg2-coords 0.40 0.80 0.70 0.90",
                "group" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
            },
            "both_RC" : " --leg-coords  0.40 0.80 0.70 0.90 --leg2-coords 0.40 0.65 0.70 0.75"
        }
    },
    "mixing13_cpv"    : {
        "AsimovA" : {
            "reactor" : {
                "valor" : " --leg-coords  0.15 0.70 0.45 0.90 --leg2-coords 0.70 0.70 0.95 0.90",
                "group" : " --leg-coords  0.15 0.70 0.45 0.90 --leg2-coords 0.70 0.70 0.95 0.90"
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.65 0.45 0.85 --leg2-coords 0.70 0.65 0.95 0.85",
                "group" : " --leg-coords  0.15 0.65 0.45 0.85 --leg2-coords 0.70 0.65 0.95 0.85"
            },
            "both_RC" : " --leg-coords  0.15 0.65 0.45 0.85 --leg2-coords 0.70 0.65 0.95 0.85"
        },
        "AsimovB" : {
            "reactor" : {
                "valor" : " --leg-coords  0.70 0.40 0.95 0.50 --leg2-coords 0.70 0.15 0.95 0.35",
                "group" : " --leg-coords  0.70 0.40 0.95 0.50 --leg2-coords 0.70 0.15 0.95 0.35"
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.70 0.40 0.95 0.50 --leg2-coords 0.70 0.15 0.95 0.35",
                "group" : " --leg-coords  0.70 0.40 0.95 0.50 --leg2-coords 0.70 0.15 0.95 0.35"
            },
            "both_RC" : " --leg-coords  0.65 0.40 0.95 0.50 --leg2-coords 0.65 0.15 0.95 0.35"
        },
        "data" : {
            "reactor" : {
                "valor" : " --leg-coords  0.70 0.80 0.95 0.90 --leg2-coords 0.70 0.55 0.95 0.75",
                "group" : " --leg-coords  0.70 0.80 0.95 0.90 --leg2-coords 0.70 0.60 0.95 0.75"
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.55 0.45 0.75",
                "group" : " --leg-coords  0.15 0.65 0.45 0.85 --leg2-coords 0.70 0.65 0.95 0.85"
            },
            "both_RC" : " --leg-coords  0.15 0.70 0.45 0.80 --leg2-coords 0.70 0.65 0.95 0.85"
        }
    },
    "mixing23_cpv"    : {
        "AsimovA" : {
            "reactor" : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75",
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75",
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75"
        },
        "AsimovB" : {
            "reactor" : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75",
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75",
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75"
        },
        "data"    : {
            "reactor" : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75",
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75",
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.15 0.57 0.45 0.75"
        }
    },
    "mixing23_dm23sq" : {
        "AsimovA" : {
            "reactor" : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90",
                "group" : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90",
                "group" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
        },
        "AsimovB" : {
            "reactor" : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90",
                "group" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90",
                "group" : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
        },
        "data"    : {
            "reactor" : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90",
                "group" : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
            },
            "t2k"     : {
                "valor" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90",
                "group" : " --leg-coords  0.15 0.75 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
            },
            "both_RC" : " --leg-coords  0.15 0.80 0.45 0.90 --leg2-coords 0.75 0.70 0.95 0.90"
        }
    }
}

_tex_slide_delimiter = "\n\n%------------------------------------------------------------\n"

#________________________________________________________________________________________________________
def main(base_dir, output_tex_name, conf_levels, prelim, skip_group, tex_image_path):

    tex_output = init_tex_slides(output_tex_name)
    
    prelim_arg = ""
    if prelim:
        prelim_arg = " --prelim"

    print_BFPs_and_CRs = []
    for params in fit_params:
        for dataset in fit_datasets:

            # Group comparisons (only for Asimov A)
            if not skip_group and dataset == "AsimovA":
                if is_2d[params]:
                    plot_group_comparison_2d(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, tex_image_path)
                else:
                    plot_group_comparison_1d(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, tex_image_path)

            # NO vs. IO; wRC and woRC on separate plots
            if is_2d[params]:
                # chisq min in true MO
                plot_valor_2d_fixed_rc(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, "global", tex_image_path)

                # chisq min across both MO
                plot_valor_2d_fixed_rc(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, "MHpair", tex_image_path)
            else:
                plot_valor_1d_fixed_rc(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, print_BFPs_and_CRs, "MHpair", tex_image_path)

            # wRC vs. woRC; NO and IO on separate plots
            if is_2d[params]:
                # chisq min in true MO
                plot_valor_2d_fixed_mo(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, "global", tex_image_path)
            else:
                plot_valor_1d_fixed_mo(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, "global", tex_image_path)

            # wRC vs. woRC; NO vs. IO (all on the same plot, rather than separately, as above)
            if not is_2d[params]:
                plot_valor_1d(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, "MHpair", tex_image_path)

            # Data vs. Asimov
            if dataset == "data":
                if is_2d[params]:
                    # One MO and RC per plot
                    # chisq min in true MO
                    plot_valor_data_vs_asimov_2d(base_dir, params, dataset, "AsimovA", conf_levels, prelim_arg, tex_output, "global", tex_image_path)
                else:
                    # Both MOs on one plot, both separate RCs
                    plot_valor_data_vs_asimov_1d(base_dir, params, dataset, "AsimovA", conf_levels, prelim_arg, tex_output, "MHpair", tex_image_path)

    # Finally print all of the best fit points and confidence regions for all fits
    slide_title = "1D best fit points and confidence regions"
    print "================================================================================"
    print "=\t", slide_title, "\t="
    print "================================================================================"
    for p in print_BFPs_and_CRs:
        print p
        tex_BFPs_and_CRs = os.linesep.join([s for s in p.splitlines() if s])    # Get rid of empty lines
        tex_BFPs_and_CRs = tex_BFPs_and_CRs.replace("\n"  , "\\\\\n")
        tex_BFPs_and_CRs = tex_BFPs_and_CRs.replace("-->" , "$\\rightarrow$")
        tex_BFPs_and_CRs = tex_BFPs_and_CRs.replace("    ", "\\hspace{0.5cm}")
        tex_BFPs_and_CRs = tex_BFPs_and_CRs.replace("\t"  , "\\hspace{0.5cm}")
        tex_output.write(_tex_slide_delimiter)
        tex_output.write("\\begin{frame}{" + slide_title +"}\n")
        tex_output.write(tex_BFPs_and_CRs +"\n")
        tex_output.write("\\end{frame}\n")

    tex_output.write("\\end{document}\n")
    print "Wrote plots to Tex file \'{}\'".format(output_tex_name)
    print "Done!"

#________________________________________________________________________________________________________
# 2D group comparisons - each plot has a single MO and RC
def plot_group_comparison_2d(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    # For the 2D group comparison plots, we want all contours from the same group to be the same colour, and each CL corresponds to a different line style
    styles = ""
    for icl in range(1, ncl+1):
        styles += " " +str(icl)

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    for reactor in fit_reactors:
        if params == "mixing13" and reactor == "reactor":
            continue
        for mo in fit_orderings:
            filename = base_dir +"/" +params +"/" +dataset +"/" +"fit_" +dataset + \
                       "_run1-10_" +mo +"_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"
            colours_valor = (" " +str(colours[0]))*ncl
            args_valor  = " --filename " +filename +" --label \"VALOR\" --colour" +colours_valor +" --style" +styles +" "   # The number of colours and styles here depends on the number of confidence levels, so is variable
            args_ptheta = ""
            args_mach3  = ""
            if fit_params_ptheta[params]:
                # p-theta stores both NO and IO fits in the same file, so no need to include the MO in the filename
                filename_ptheta = base_dir +"/ptheta/cont-20200511-01-02-eETh/contour_" + \
                                  fit_params_ptheta[params] +reac_ptheta[reactor] +"_20200511-01-eETh-" +dataset +".root"
                colours_ptheta = (" " +str(colours[2]))*ncl
                args_ptheta = " --compare-ptheta " +filename_ptheta +" --ptheta-label \"p-theta (E-#theta)\" --ptheta-colour" +colours_ptheta +" --ptheta-style" +styles +" "
            if fit_params_mach3[params]:
                filename_mach3 = base_dir +"/mach3/contours_" +dataset +"_" +fit_params_mach3[params] +"_" +mo_mach3[mo] +reac_mach3[reactor] +".root"
                colours_mach3 = (" " +str(colours[1]))*ncl
                args_mach3 = " --compare-mach3 " +filename_mach3 + " --mach3-label \"MaCh3\" --mach3-colour" +colours_mach3 +" --mach3-style" +styles +" "

            if args_ptheta != "" or args_mach3 != "":
                # Make args
                axis_title_args = None
                if "dm23sq" in params:
                    axis_title_args = title_args[params][mo]
                else:
                    axis_title_args = title_args[params]
                axes_range_args = ""
                if range_args[params]:
                    if type(range_args[params]) == type({}):
                        axes_range_args = range_args[params][reactor]
                    else:
                        axes_range_args = range_args[params]
                leg_coords = ""
                if leg_coords_args[params][dataset][reactor]["group"]:
                    leg_coords = leg_coords_args[params][dataset][reactor]["group"]
                dchisq_calc = "global"
                dchisq_calc_arg = " --chisq-calc " +dchisq_calc
                group_names = "VALOR"
                if args_mach3 != "":
                    group_names += "_vs_MaCh3"
                if args_ptheta != "":
                    group_names += "_vs_p-theta"
                output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name[reactor] +"__" +mo_out_name[mo] +"__" +group_names
                output_name_arg = " --output " +output_name
                standard_args_2d = " --best-fit-label \"Best Fit\"  --best-fit-style 22"
                command = script +args_valor +args_ptheta +args_mach3 +standard_args + \
                          standard_args_2d +axis_title_args +axes_range_args +leg_coords + \
                          conf_levels_arg +dchisq_calc_arg +prelim_arg +output_name_arg +reac_band_arg

                # Run plotting script
                slide_title = "2D " +params.replace("_","-") +" " +dataset +" " +reactor +" " +mo +" " +"Group Comparison"
                print "================================================================================"
                print "=\t", slide_title, "\t="
                print "================================================================================"
                print "Running:\n", command, "\n"
                captured_output = ""
                with suppress_stdout_stderr():
                    captured_output = os.popen(command).read()
                print get_filtered_output(captured_output, dchisq_calc)    # Just print the output, as we don't calculate the CRs and BF points for 2D plots
                print "\n\n\n"

                # Make a slide with the produced plot
                make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
# 1D group comparisons - each plot has a single RC, but both MOs
def plot_group_comparison_1d(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    mo = "normal"   # Although we specify NO here, the --chisq-calc MHpair argument we use later makes the plotting script automatically look for the IO file too.
    for reactor in fit_reactors:
        if params == "mixing13" and reactor == "reactor":
            continue
        filename = base_dir +"/" +params +"/" +dataset +"/" +"fit_" +dataset + \
                   "_run1-10_" +mo +"_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"
        colours_valor = (" " +str(colours[0]))*2   # The number of colours and styles here depends on the number of MOs being plot, so 2
        args_valor  = " --filename " +filename +" --label \"VALOR\" --colour" +colours_valor +" --style 1 2"
        args_ptheta = ""
        args_mach3  = ""
        if fit_params_ptheta[params]:
            # p-theta stores both NO and IO fits in the same file, so no need to include the MO in the filename
            filename_ptheta = base_dir +"/ptheta/cont-20200511-01-02-eETh/contour_" + \
                              fit_params_ptheta[params] +reac_ptheta[reactor] +"_20200511-01-eETh-" +dataset +".root"
            colours_ptheta = (" " +str(colours[2]))*2
            args_ptheta = " --compare-ptheta " +filename_ptheta +" --ptheta-label \"p-theta (E-#theta)\" --ptheta-colour" +colours_ptheta +" --ptheta-style 1 2 "
        if fit_params_mach3[params]:
            filename_mach3 = base_dir +"/mach3/contours_" +dataset +"_" +fit_params_mach3[params] +"_" +mo_mach3[mo] +reac_mach3[reactor] +".root"
            colours_mach3 = (" " +str(colours[1]))*2
            args_mach3 = " --compare-mach3 " +filename_mach3 + " --mach3-label \"MaCh3\" --mach3-colour" +colours_mach3 +" --mach3-style 1 2"

        if args_ptheta != "" or args_mach3 != "":
            # Make args
            axis_title_args = None
            if "dm23sq" in params:
                axis_title_args = title_args[params]["both"]
            else:
                axis_title_args = title_args[params]
            axes_range_args = ""
            if range_args[params]:
                if type(range_args[params]) == type({}):
                    axes_range_args = range_args[params][reactor]
                else:
                    axes_range_args = range_args[params]
            leg_coords = ""
            if leg_coords_args[params][dataset][reactor]["group"]:
                leg_coords = leg_coords_args[params][dataset][reactor]["group"]
            dchisq_calc = "MHpair"
            dchisq_calc_arg = " --chisq-calc " +dchisq_calc
            group_names = "VALOR"
            if args_mach3 != "":
                group_names += "_vs_MaCh3"
            if args_ptheta != "":
                group_names += "_vs_p-theta"
            output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name[reactor] +"__NO_vs_IO__" +group_names
            output_name_arg = " --output " +output_name
            command = script +args_valor +args_ptheta +args_mach3 +standard_args + \
                      axis_title_args +axes_range_args +leg_coords + \
                      conf_levels_arg +dchisq_calc_arg +prelim_arg +output_name_arg +reac_band_arg

            # Run plotting script
            slide_title = "1D " +params.replace("_","-") +" " +dataset +" " +reactor +" " +"Group Comparison"
            print "================================================================================"
            print "=\t", slide_title, "\t="
            print "================================================================================"
            print "Running:\n", command, "\n"
            captured_output = ""
            with suppress_stdout_stderr():
                captured_output = os.popen(command).read()
            print get_filtered_output(captured_output, dchisq_calc)
            print "\n\n\n"

            # Make a slide with the produced plot
            make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def plot_valor_2d_fixed_rc(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, dchisq_calc, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    # For the 2D plots, we want all contours from the same fit to be the same colour, and each CL corresponds to a different line style
    # The order is fit0+CL0, fit0+CL1, ..., fit0+CLN, fit1+CL0, ..., fitN+CLN
    # Here fit0 is NO and fit1 is IO
    styles = ""
    cols   = ""
    for imo in range(0, 2):
        for icl in range(1, ncl+1):
            cols   += " " +str(colours[imo])
            styles += " " +str(icl)

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    for reactor in fit_reactors:
        if params == "mixing13" and reactor == "reactor":
            continue

        filename = base_dir +"/" +params +"/" +dataset +"/" +"fit_" +dataset + "_run1-10_normal_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"
        if dchisq_calc == "global" or dchisq_calc == "local":
            filename += " " +filename.replace("normal", "inverted")
        args_valor  = " --filename " +filename +" --label \"Normal\" \"Inverted\" --colour" +cols +" --style" +styles

        # Make args
        axis_title_args = None
        if "dm23sq" in params:
            axis_title_args = title_args[params]["both"]
        else:
            axis_title_args = title_args[params]

        axes_range_args = ""
        if range_args[params]:
            if type(range_args[params]) == type({}):
                axes_range_args = range_args[params][reactor]
            else:
                axes_range_args = range_args[params]

        leg_coords = ""
        if leg_coords_args[params][dataset][reactor]["valor"]:
            leg_coords = leg_coords_args[params][dataset][reactor]["valor"]

        dchisq_calc_arg = " --chisq-calc " +dchisq_calc
        dchisq_calc_desc = ""
        if dchisq_calc == "MHpair":
            dchisq_calc_desc = "__dchisq_min_across_MH"
        elif dchisq_calc == "local":
            dchisq_calc_desc = "__dchisq_min_wrt_min_of_all_inputs"

        output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name[reactor] +"__" +mo_out_name["normal"] +"_vs_" +mo_out_name["inverted"] +dchisq_calc_desc
        output_name_arg = " --output " +output_name

        standard_args_2d = " --best-fit-label \"Best Fit\"  --best-fit-style 22 22"     # The number of bf style numbers depends on the number of MOs or RCs plot

        command = script +args_valor +standard_args + \
                  standard_args_2d +axis_title_args +axes_range_args +leg_coords + \
                  conf_levels_arg +dchisq_calc_arg +reac_band_arg +prelim_arg +output_name_arg

        # Run plotting script
        slide_title = "2D " +params.replace("_","-") +" " +dataset +" " +reactor +" NO vs. IO" +dchisq_calc_desc.replace("_", " ")
        print "================================================================================"
        print "=\t", slide_title, "\t="
        print "================================================================================"
        print "Running:\n", command, "\n"
        captured_output = ""
        with suppress_stdout_stderr():
            captured_output = os.popen(command).read()
        print get_filtered_output(captured_output, dchisq_calc)    # Just print the output, as we don't calculate the CRs and BF points for 2D plots
        print "\n\n\n"

        # Make a slide with the produced plot
        make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def plot_valor_1d_fixed_rc(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, print_BFPs_and_CRs, dchisq_calc, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    mo = "normal"   # Although we specify NO here, the --chisq-calc MHpair argument we use later makes the plotting script automatically look for the IO file too.
    for reactor in fit_reactors:
        if params == "mixing13" and reactor == "reactor":
            continue
        filename = base_dir +"/" +params +"/" +dataset +"/" +"fit_" +dataset + \
                   "_run1-10_" +mo +"_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"
        colours_valor = (" " +str(colours[0]))*2   # The number of colours and styles here depends on the number of MOs being plot, so 2
        args_valor  = " --filename " +filename +" --colour" +colours_valor +" --style 1 2"

        axis_title_args = None
        if "dm23sq" in params:
            axis_title_args = title_args[params]["both"]
        else:
            axis_title_args = title_args[params]

        axes_range_args = ""
        if range_args[params]:
            if type(range_args[params]) == type({}):
                axes_range_args = range_args[params][reactor]
            else:
                axes_range_args = range_args[params]

        leg_coords = ""
        if leg_coords_args[params][dataset][reactor]["valor"]:
            leg_coords = leg_coords_args[params][dataset][reactor]["valor"]

        dchisq_calc_arg = " --chisq-calc " +dchisq_calc

        output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name[reactor] +"__NO_vs_IO"
        output_name_arg = " --output " +output_name
        command = script +args_valor +standard_args + \
                  axis_title_args +axes_range_args +leg_coords + \
                  conf_levels_arg +dchisq_calc_arg +prelim_arg +output_name_arg +reac_band_arg

        # Run plotting script
        slide_title = "1D " +params.replace("_","-") +" " +dataset +" " +reactor +" NO vs. IO"
        print "================================================================================"
        print "=\t", slide_title, "\t="
        print "================================================================================"
        print "Running:\n", command, "\n"
        captured_output = ""
        with suppress_stdout_stderr():
            captured_output = os.popen(command).read()
        text_to_return = get_filtered_output(captured_output, dchisq_calc)
        print "\n\n\n"

        if text_to_return and text_to_return != "":
            print_BFPs_and_CRs.append(text_to_return)

        # Make a slide with the produced plot
        make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def plot_valor_2d_fixed_mo(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, dchisq_calc, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    # For the 2D plots, we want all contours from the same fit to be the same colour, and each CL corresponds to a different line style
    # The order is fit0+CL0, fit0+CL1, ..., fit0+CLN, fit1+CL0, ..., fitN+CLN
    # Here fit0 is wRC and fit1 is woRC
    styles = ""
    cols   = ""
    for ireac in range(0, 2):
        for icl in range(1, ncl+1):
            cols   += " " +str(colours[ireac])
            styles += " " +str(icl)

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    for mo in fit_orderings:

        filename = base_dir +"/" +params +"/" +dataset +"/" +"fit_" +dataset + "_run1-10_" +mo +"_" +params +"_reactor_0_1_result" +smear_arg +".root"
        filename += " " +filename.replace("reactor", "t2k")
        args_valor  = " --filename " +filename +" --label \"T2K + Reactor\" \"T2K-only\" --colour" +cols +" --style" +styles

        # Make args
        axis_title_args = None
        if "dm23sq" in params:
            axis_title_args = title_args[params][mo]
        else:
            axis_title_args = title_args[params]

        axes_range_args = ""
        if range_args[params]:
            if type(range_args[params]) == type({}):
                axes_range_args = range_args[params]["t2k"]
            else:
                axes_range_args = range_args[params]

        leg_coords = ""
        if leg_coords_args[params][dataset]["both_RC"]:
            leg_coords = leg_coords_args[params][dataset]["both_RC"]

        dchisq_calc_arg = " --chisq-calc " +dchisq_calc

        output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name["reactor"]+"_vs_" +reac_out_name["t2k"] +"__" +mo_out_name[mo]
        output_name_arg = " --output " +output_name

        standard_args_2d = " --best-fit-label \"Best Fit\"  --best-fit-style 22 22"     # The number of bf style numbers depends on the number of MOs or RCs plot

        command = script +args_valor +standard_args + \
                  standard_args_2d +axis_title_args +axes_range_args +leg_coords + \
                  conf_levels_arg +dchisq_calc_arg +reac_band_arg +prelim_arg +output_name_arg

        # Run plotting script
        slide_title = "2D " +params.replace("_","-") +" " +dataset +" wRC vs. woRC " +mo_out_name[mo]
        print "================================================================================"
        print "=\t", slide_title, "\t="
        print "================================================================================"
        print "Running:\n", command, "\n"
        captured_output = ""
        with suppress_stdout_stderr():
            captured_output = os.popen(command).read()
        print get_filtered_output(captured_output, dchisq_calc)    # Just print the output, as we don't calculate the CRs and BF points for 2D plots
        print "\n\n\n"

        # Make a slide with the produced plot
        make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def plot_valor_1d_fixed_mo(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, dchisq_calc, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    for mo in fit_orderings:
        if params == "mixing13":
            continue
        filename = base_dir +"/" +params +"/" +dataset +"/" +"fit_" +dataset + \
                   "_run1-10_" +mo +"_" +params +"_reactor_0_1_result" +smear_arg +".root"
        filename += " " +filename.replace("reactor", "t2k")
        colours_valor = (" " +str(colours[0]))*2   # The number of colours and styles here depends on the number of RCs being plot, so 2
        args_valor  = " --filename " +filename +" --label \"T2K + Reactor\" \"T2K-only\" --colour" +colours_valor +" --style 1 2"

        axis_title_args = None
        if "dm23sq" in params:
            axis_title_args = title_args[params][mo]
        else:
            axis_title_args = title_args[params]

        axes_range_args = ""
        if range_args[params]:
            if type(range_args[params]) == type({}):
                axes_range_args = range_args[params]["t2k"]
            else:
                axes_range_args = range_args[params]

        leg_coords = ""
        if leg_coords_args[params][dataset]["both_RC"]:
            leg_coords = leg_coords_args[params][dataset]["both_RC"]

        dchisq_calc_arg = " --chisq-calc " +dchisq_calc

        output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name["reactor"] +"_vs_" +reac_out_name["t2k"] +"__" +mo_out_name[mo]
        output_name_arg = " --output " +output_name
        command = script +args_valor +standard_args + \
                  axis_title_args +axes_range_args +leg_coords + \
                  conf_levels_arg +dchisq_calc_arg +prelim_arg +output_name_arg +reac_band_arg

        # Run plotting script
        slide_title = "1D " +params.replace("_","-") +" " +dataset +" " +reac_out_name["reactor"] +" vs " +reac_out_name["t2k"] +" " +mo_out_name[mo]
        print "================================================================================"
        print "=\t", slide_title, "\t="
        print "================================================================================"
        print "Running:\n", command, "\n"
        captured_output = ""
        with suppress_stdout_stderr():
            captured_output = os.popen(command).read()
        print get_filtered_output(captured_output, dchisq_calc)
        print "\n\n\n"

        # Make a slide with the produced plot
        make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def plot_valor_1d(base_dir, params, dataset, conf_levels, prelim_arg, tex_output, dchisq_calc, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    if params == "mixing13":
        return
    filename = base_dir +"/" +params +"/" +dataset +"/" +"fit_" +dataset + \
               "_run1-10_normal_" +params +"_reactor_0_1_result" +smear_arg +".root"
    filename += " " +filename.replace("reactor", "t2k")
    colours_valor = (" " +str(colours[0]))*2 + (" " +str(colours[1]))*2  # The number of colours and styles here depends on both the number of RCs and MOs being plot, so 4
    args_valor  = " --filename " +filename +" --label \"T2K + Reactor\" \"T2K-only\" --colour" +colours_valor +" --style 1 2 1 2"

    axis_title_args = None
    if "dm23sq" in params:
        axis_title_args = title_args[params]["both"]
    else:
        axis_title_args = title_args[params]

    axes_range_args = ""
    if range_args[params]:
        if type(range_args[params]) == type({}):
            axes_range_args = range_args[params]["t2k"]
        else:
            axes_range_args = range_args[params]

    leg_coords = ""
    if leg_coords_args[params][dataset]["both_RC"]:
        leg_coords = leg_coords_args[params][dataset]["both_RC"]

    dchisq_calc_arg = " --chisq-calc " +dchisq_calc

    output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name["reactor"] +"_vs_" +reac_out_name["t2k"] +"__" +mo_out_name["normal"] +"_vs_" +mo_out_name["inverted"]
    output_name_arg = " --output " +output_name
    command = script +args_valor +standard_args + \
              axis_title_args +axes_range_args +leg_coords + \
              conf_levels_arg +dchisq_calc_arg +prelim_arg +output_name_arg +reac_band_arg

    # Run plotting script
    slide_title = "1D " +params.replace("_","-") +" " +dataset +" " +reac_out_name["reactor"] +" vs " +reac_out_name["t2k"] +", " +mo_out_name["normal"] +" vs " +mo_out_name["inverted"]
    print "================================================================================"
    print "=\t", slide_title, "\t="
    print "================================================================================"
    print "Running:\n", command, "\n"
    captured_output = ""
    with suppress_stdout_stderr():
        captured_output = os.popen(command).read()
    print get_filtered_output(captured_output, dchisq_calc)
    print "\n\n\n"

    # Make a slide with the produced plot
    make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def plot_valor_data_vs_asimov_2d(base_dir, params, dataset, asimov_dataset, conf_levels, prelim_arg, tex_output, dchisq_calc, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    # For the 2D plots, we want all contours from the same fit to be the same colour, and each CL corresponds to a different line style
    # The order is fit0+CL0, fit0+CL1, ..., fit0+CLN, fit1+CL0, ..., fitN+CLN
    # Here fit0 is NO and fit1 is IO
    styles = ""
    cols   = ""
    for imo in range(0, 2):
        for icl in range(1, ncl+1):
            cols   += " " +str(colours[imo])
            styles += " " +str(icl)

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    for reactor in fit_reactors:
        if params == "mixing13" and reactor == "reactor":
            continue

        for mo in fit_orderings:


            filename_data   = base_dir +"/" +params +"/" +dataset        +"/" +"fit_" +dataset        +"_run1-10_" +mo +"_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"
            filename_asimov = base_dir +"/" +params +"/" +asimov_dataset +"/" +"fit_" +asimov_dataset +"_run1-10_" +mo +"_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"

            args_valor  = " --filename " +filename_data +" " +filename_asimov +" --label \"Data\" \"Asimov\" --colour" +cols +" --style" +styles

            # Make args
            axis_title_args = None
            if "dm23sq" in params:
                axis_title_args = title_args[params][mo]
            else:
                axis_title_args = title_args[params]

            axes_range_args = ""
            if range_args[params]:
                if type(range_args[params]) == type({}):
                    axes_range_args = range_args[params][reactor]
                else:
                    axes_range_args = range_args[params]

            leg_coords = ""
            if leg_coords_args[params][dataset][reactor]["valor"]:
                leg_coords = leg_coords_args[params][dataset][reactor]["valor"]

            dchisq_calc_arg = " --chisq-calc " +dchisq_calc
            dchisq_calc_desc = ""
            if dchisq_calc == "MHpair":
                dchisq_calc_desc = "__dchisq_min_across_MH"
            elif dchisq_calc == "local":
                dchisq_calc_desc = "__dchisq_min_wrt_min_of_all_inputs"

            output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name[reactor] +"__" +mo_out_name[mo] +"__data_vs_Asimov__" +dchisq_calc_desc
            output_name_arg = " --output " +output_name

            standard_args_2d = " --best-fit-label \"Best Fit\"  --best-fit-style 22 22"     # The number of bf style numbers depends on the number of MOs or RCs plot

            command = script +args_valor +standard_args + \
                      standard_args_2d +axis_title_args +axes_range_args +leg_coords + \
                      conf_levels_arg +dchisq_calc_arg +reac_band_arg +prelim_arg +output_name_arg

            # Run plotting script
            slide_title = "2D " +params.replace("_","-") +" " +dataset +" " +reactor + " " + mo +" data vs. Asimov " +dchisq_calc_desc.replace("_", " ")
            print "================================================================================"
            print "=\t", slide_title, "\t="
            print "================================================================================"
            print "Running:\n", command, "\n"
            captured_output = ""
            with suppress_stdout_stderr():
                captured_output = os.popen(command).read()
            print get_filtered_output(captured_output, dchisq_calc)    # Just print the output, as we don't calculate the CRs and BF points for 2D plots
            print "\n\n\n"

            # Make a slide with the produced plot
            make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def plot_valor_data_vs_asimov_1d(base_dir, params, dataset, asimov_dataset, conf_levels, prelim_arg, tex_output, dchisq_calc, tex_image_path):
    conf_levels_arg = " --conf-levels "
    ncl = len(conf_levels)
    for cl in conf_levels:
        conf_levels_arg += cl +" "

    smear_arg = ""
    if "dm23sq" in params:
        smear_arg = "_smeared"

    reac_band_arg = ""
    if "mixing13" in params:
        reac_band_arg = " --draw-reactor"

    mo = "normal"   # Although we specify NO here, the --chisq-calc MHpair argument we use later makes the plotting script automatically look for the IO file too.
    for reactor in fit_reactors:
        if params == "mixing13" and reactor == "reactor":
            continue
        filename_data   = base_dir +"/" +params +"/" +dataset        +"/" +"fit_" +dataset        +"_run1-10_" +mo +"_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"
        filename_asimov = base_dir +"/" +params +"/" +asimov_dataset +"/" +"fit_" +asimov_dataset +"_run1-10_" +mo +"_" +params +"_" +reactor +"_0_1_result" +smear_arg +".root"
        colours_valor = (" " +str(colours[0]))*2 +(" " +str(colours[1]))*2
        args_valor  = " --filename " +filename_data +" " +filename_asimov +" --label \"Data\" \"Asimov\" --colour" +colours_valor +" --style 1 2 1 2"

        axis_title_args = None
        if "dm23sq" in params:
            axis_title_args = title_args[params]["both"]
        else:
            axis_title_args = title_args[params]

        axes_range_args = ""
        if range_args[params]:
            if type(range_args[params]) == type({}):
                axes_range_args = range_args[params][reactor]
            else:
                axes_range_args = range_args[params]

        leg_coords = ""
        if params == "mixing13" and leg_coords_args[params][dataset]["t2k"]["valor"]:
            leg_coords = leg_coords_args[params][dataset]["t2k"]["valor"]
        elif leg_coords_args[params][dataset]["both_RC"]:     # Yeah, we've only got one RC per plot here, but this does the job anyway...
                leg_coords = leg_coords_args[params][dataset]["both_RC"]

        dchisq_calc_arg = " --chisq-calc " +dchisq_calc

        output_name = "contour_run10__" +params +"__" +dataset +"__" +reac_out_name[reactor] +"__NO_vs_IO__data_vs_Asimov"
        output_name_arg = " --output " +output_name
        command = script +args_valor +standard_args + \
                  axis_title_args +axes_range_args +leg_coords + \
                  conf_levels_arg +dchisq_calc_arg +prelim_arg +output_name_arg +reac_band_arg

        # Run plotting script
        slide_title = "1D " +params.replace("_","-") +" " +dataset +" " +reactor +" NO vs. IO, data vs. Asimov"
        print "================================================================================"
        print "=\t", slide_title, "\t="
        print "================================================================================"
        print "Running:\n", command, "\n"
        captured_output = ""
        with suppress_stdout_stderr():
            captured_output = os.popen(command).read()
        print get_filtered_output(captured_output, dchisq_calc)
        print "\n\n\n"

        # Make a slide with the produced plot
        make_tex_image_slide(tex_output, slide_title, output_name, tex_image_path)

#________________________________________________________________________________________________________
def get_filtered_output(captured_tex, dchisq_calc_type):    # Might not always want to print CRs (e.g. for 2D)
    print_str = ""

    chisq_calc_str = None
    if dchisq_calc_type == "global":
        chisq_calc_str = "w.r.t. true MO:"
    elif dchisq_calc_type == "local":
        chisq_calc_str = "w.r.t. min chisq across all input files:"
    elif dchisq_calc_type == "MHpair":
        chisq_calc_str = "w.r.t. min chisq across both MOs:"

    ignore_lines_containing = ["Number of chisq hists to plot", "======", "has been created"]

    captured = captured_tex.split("\n")
    for l in captured:
        if "fit from file" in l:
            fit_name = l.split("/")[-1]
            fit_name = fit_name.split("_0_1_result")[0]
            fit_name = fit_name.split("fit_")[-1]
            fit_name = fit_name.replace("_run1-10","")
            fit_name = fit_name.replace("_"," ")
            print_str += "\n\n"+ fit_name +"; " +chisq_calc_str +"\n"
        elif not any(ignore in l for ignore in ignore_lines_containing):
            if "ERROR" in l or "WARN" in l or "NOTE" in l or "FATAL" in l:
                print l
            else:
                print_str += "\t" +l +"\n"
    return print_str



#________________________________________________________________________________________________________
# For me, running pyvalor produces a large amount of meaningless errors about redefining `TColor`s.
# These don't impact the macros we're running in any meaningful way expect for flooding the terminal output with rubbish.
# So, define a context manager to suppress stdout and stderr.
# Source: https://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

#________________________________________________________________________________________________________
def init_tex_slides(output_tex_name):
    tex_output = open(output_tex_name,"w")
    tex_output.write(
        """\\documentclass{beamer}
\\usepackage{pgfpages}
\\usepackage{graphicx}
\\usepackage{booktabs}
\\usepackage{multirow}
\\usetheme{Madrid}
\\pgfpagesuselayout{resize to}[a4paper,landscape,border shrink=5mm]
\\beamertemplatenavigationsymbolsempty   % get rid of navigation bar

\\begin{document}
           
\\title[VALOR T2K Run 1-10 Fits] { \\begin{center} {VALOR T2K Run 1-10 Fits} \\end{center} }
\\author[F. Bench] {Francis Bench, Maria Antonova, Jaiden Parlone}
\\date{ """ + datetime.today().strftime('%Y-%m-%d') + " \\ \\vspace{2cm}}"
    )
    tex_output.write(_tex_slide_delimiter)
    tex_output.write(
        """\\begin{frame}[plain]
\\titlepage
\\end{frame}"""
    )
    return tex_output

#________________________________________________________________________________________________________
def make_tex_image_slide(tex_output, slide_title, pdf_name, tex_image_path):
    tex_output.write(_tex_slide_delimiter)
    tex_output.write("\\begin{frame}{" + slide_title +"}\n")
    tex_output.write("\\begin{center}\n")
    tex_output.write("\\includegraphics[width=0.9\\textwidth]{{" +tex_image_path +pdf_name +"}.pdf}\n")
    tex_output.write("\\end{center}\n")
    tex_output.write("\\end{frame}\n")

#________________________________________________________________________________________________________
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("--base-dir"       , dest="base_dir"                                            , help="Directory containing the fit directories for each param: e.g. base_dir/cpv/AsimovA/fit.root.")
    parser.add_argument("--output-tex"     , dest="output_tex"            , default="fits.tex"          , help="Name of output Tex file containing a slide for each plot produced.")
    parser.add_argument("--tex-image-path" , dest="tex_image_path"        , default="images/run10_fits/", help="Where you put your images in relation to the tex file")
    parser.add_argument("--conf-levels"    , dest="conf_levels", nargs="+", default=["90pct", "1sigma"] , help="The confidence levels to plot const chisq lines/contours for. See ConfLevel::AsString() for the choices")
    parser.add_argument("--prelim"         , dest="prelim"     , action="store_true"                    , help="Make a box with a caveat about the plot being a preliminary result."    )
    parser.add_argument("--skip-group-comp", dest="skip_group" , action="store_true"                    , help="If true, doesn't plot the group comparisons (useful if e.g. you want to plot at CLs which the other groups haven't made contours for)"    )

    args = parser.parse_args()

    main(args.base_dir, args.output_tex, args.conf_levels, args.prelim, args.skip_group, args.tex_image_path)

