from collections import OrderedDict
import json

def getGenXsecs(dir_xsec, sqrtS = "13TeV", gen = "amcPythia"):
    results = OrderedDict()
    for cat in ["born", "dressed"]:
        results[cat] = OrderedDict()
        for ch in ["zmm", "wmp", "wmm", "zee", "wep", "wem"]:
            fname = "%s/%s/GEN_%s_%s_%s_%s/acceptance.txt" % (dir_xsec, sqrtS, ch, sqrtS, gen, cat)

            # loop over file and get the last line
            with open(fname, "r") as f_in:
                lines = filter(None, [line.rstrip() for line in f_in])
                xsec = lines[-1].split()[1]

                print("xsec for %s is %s" % (fname, xsec))
                
                results[cat][ch] = float(xsec)
    
    return results

if __name__ == "__main__":
    
    results = OrderedDict()
    dir_xsec = "/eos/uscms/store/user/yofeng/LowPUResults/TestAccept/"

    results["13TeV"] = getGenXsecs(dir_xsec, "13TeV")
    results["5TeV"] = getGenXsecs(dir_xsec, "5TeV")
    print(results)
    
    with open("xsecs.json", "w") as f_out:
        json.dump(results, f_out, indent=4, sort_keys=False)