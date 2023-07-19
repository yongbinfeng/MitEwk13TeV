from collections import OrderedDict
import json

def getGenXsecs(dir_xsec, sqrtS = "13TeV", gen = "amcPythia", doFiducial = False):
    # total inclusive xsec, the xsec is printed in the 1st line (starting with 0)
    # total fiducial xsec, the xsec is printed in the last line.
    idx = -1 if doFiducial else 1
    results = OrderedDict()
    for cat in ["born", "dressed"]:
        results[cat] = OrderedDict()
        for ch in ["zmm", "wmp", "wmm", "zee", "wep", "wem"]:
            fname = "%s/%s/GEN_%s_%s_%s_%s/acceptance.txt" % (dir_xsec, sqrtS, ch, sqrtS, gen, cat)

            # loop over file and get the last line
            with open(fname, "r") as f_in:
                lines = filter(None, [line.rstrip() for line in f_in])
                xsec = lines[idx].split()[1]

                print("xsec for %s is %s" % (fname, xsec))
                
                results[cat][ch] = float(xsec)
    
    return results

if __name__ == "__main__":
    
    doFiducial = True
    suffix = "fid" if doFiducial else "inc"
    
    results = OrderedDict()
    dir_xsec = "/eos/uscms/store/user/yofeng/LowPUResults/TestAccept/"

    results["13TeV"] = getGenXsecs(dir_xsec, "13TeV", doFiducial=doFiducial)
    results["5TeV"] = getGenXsecs(dir_xsec, "5TeV", doFiducial=doFiducial)
    print(results)
    
    with open("xsecs_%s.json"%suffix, "w") as f_out:
        json.dump(results, f_out, indent=4, sort_keys=False)
