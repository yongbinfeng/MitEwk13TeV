#include "EnergyScaleCorrection.h"

//#include "FWCore/ParameterSet/interface/FileInPath.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>

EnergyScaleCorrection::EnergyScaleCorrection(const std::string& correctionFileName, unsigned int genSeed)
    : smearingType_(ECALELF)
{

    if (!correctionFileName.empty()) {
        std::string filename = correctionFileName + "_scales.dat";
        readScalesFromFile(filename);
        if (scales_.empty()) {
            std::cerr << "[ERROR] scale correction map empty" << std::endl;
            exit(1);
        }
    }

    if (!correctionFileName.empty()) {
        std::string filename = correctionFileName + "_smearings.dat";
        readSmearingsFromFile(filename);
        if (smearings_.empty()) {
            std::cerr << "[ERROR] scale correction map empty" << std::endl;
            exit(1);
        }
    }
}

EnergyScaleCorrection::EnergyScaleCorrection(const std::string& correctionFileName, FileFormat type)
    : smearingType_(type)
{
    // cout << "smearing type  " << smearingType_ << endl;
    cout << "Read the scales " << endl;
    if (!correctionFileName.empty()) {
        std::string filename = correctionFileName + "_scales.dat";
        readScalesFromFile(filename);
        if (scales_.empty()) {
            std::cerr << "[ERROR] scale correction map empty" << std::endl;
            exit(1);
        }
    }

    cout << "Read the smearings " << endl;
    if (!correctionFileName.empty()) {
        std::string filename = correctionFileName + "_smearings.dat";
        readSmearingsFromFile(filename);
        if (smearings_.empty()) {
            std::cerr << "[ERROR] scale correction map empty" << std::endl;
            exit(1);
        }
    }
}

float EnergyScaleCorrection::scaleCorr(unsigned int runNumber, double et, double eta, double r9,
    unsigned int gainSeed, std::bitset<kErrNrBits> uncBitMask) const
{
    // std::cout << "scale corr" << std::endl;
    const ScaleCorrection* scaleCorr = getScaleCorr(runNumber, et, eta, r9, gainSeed);
    if (scaleCorr != nullptr)
        return scaleCorr->scale();
    else
        return kDefaultScaleVal_;
}

float EnergyScaleCorrection::scaleCorrUncert(unsigned int runNumber, double et, double eta, double r9,
    unsigned int gainSeed, int uncBitMask) const
// unsigned int gainSeed, std::bitset<kErrNrBits> uncBitMask) const
{

    const ScaleCorrection* scaleCorr = getScaleCorr(runNumber, et, eta, r9, gainSeed);
    // std::cout << "bits " << uncBitMask << std::endl;
    if (scaleCorr != nullptr)
        return scaleCorr->scaleErr(uncBitMask);
    else
        return 0.;
}

float EnergyScaleCorrection::smearingSigma(int runnr, double et, double eta, double r9,
    unsigned int gainSeed, ParamSmear par,
    float nSigma) const
{
    // std::cout << "do smear" << std::endl;
    // smearCorr->print();
    if (par == kRho)
        return smearingSigma(runnr, et, eta, r9, gainSeed, nSigma, 0.);
    if (par == kPhi)
        return smearingSigma(runnr, et, eta, r9, gainSeed, 0., nSigma);
    return smearingSigma(runnr, et, eta, r9, gainSeed, 0., 0.);
}

float EnergyScaleCorrection::smearingSigma(int runnr, double et, double eta, double r9,
    unsigned int gainSeed, float nrSigmaRho,
    float nrSigmaPhi) const
{
    const SmearCorrection* smearCorr = getSmearCorr(runnr, et, eta, r9, gainSeed);
    // smearCorr->print();
    // std::cout << "rho " << nrSigmaRho << "  phi " << nrSigmaPhi << std::endl;
    // std::cout << "smear correction " << smearCorr << std::endl;
    if (smearCorr != nullptr)
        return smearCorr->sigma(et, nrSigmaRho, nrSigmaPhi);
    else
        return kDefaultSmearVal_;
}

const EnergyScaleCorrection::ScaleCorrection*
EnergyScaleCorrection::getScaleCorr(unsigned int runnr, double et, double eta, double r9,
    unsigned int gainSeed) const
{

    // buld the category based on the values of the object
    CorrectionCategory category(runnr, et, eta, r9, gainSeed);
    auto result = std::equal_range(scales_.begin(), scales_.end(), category, Sorter<CorrectionCategory, ScaleCorrection>());
    auto nrFound = std::distance(result.first, result.second);
    if (nrFound == 0) {
        // std::cout << "hello" << std::endl;
        std::cout << "Scale category not found: " << category << " Returning uncorrected value." << std::endl;
        return nullptr;
    } else if (nrFound > 1) {
        std::ostringstream foundCats;
        for (auto it = result.first; it != result.second; ++it) {
            foundCats << "    " << it->first << std::endl;
        }
        //throw cms::Exception("ConfigError") <<" scale error category "<<category<<" has "<<nrFound<<" entries "<<std::endl<<foundCats.str();
        std::cerr << "Config error" << std::endl;
    }
    //validate the result, just to be sure
    if (!result.first->first.inCategory(runnr, et, eta, r9, gainSeed)) {
        //throw cms::Exception("LogicError") <<" error found scale category "<<result.first->first<<" that does not contain run "<<runnr<<" et "<<et<<" eta "<<eta<<" r9 "<<r9<<" gain seed "<<gainSeed;
        std::cerr << "Logic error" << std::endl;
    }
    return &result.first->second;
}

const EnergyScaleCorrection::SmearCorrection*
EnergyScaleCorrection::getSmearCorr(unsigned int runnr, double et, double eta, double r9,
    unsigned int gainSeed) const
{

    // buld the category based on the values of the object
    CorrectionCategory category(runnr, et, eta, r9, gainSeed);
    // std::cout << "category we're correcting for " << category << std::endl;
    auto result = std::equal_range(smearings_.begin(), smearings_.end(), category, Sorter<CorrectionCategory, SmearCorrection>());
    auto nrFound = std::distance(result.first, result.second);
    if (nrFound == 0) {
        // std::cout << "hello" << std::endl;
        //edm::LogInfo("EnergyScaleCorrection") << "Smear category not found: " << category << " Returning uncorrected value.";
        return nullptr;
    } else if (nrFound > 1) {
        std::ostringstream foundCats;
        for (auto it = result.first; it != result.second; ++it) {
            foundCats << "    " << it->first << std::endl;
        }
        //throw cms::Exception("ConfigError") <<" error smear category "<<category<<" has "<<nrFound<<" entries "<<std::endl<<foundCats.str();
        std::cerr << "ConfigError" << std::endl;
    }
    //validate the result, just to be sure
    if (!result.first->first.inCategory(runnr, et, eta, r9, gainSeed)) {
        //throw cms::Exception("LogicError") <<" error found smear category "<<result.first->first<<" that does not contain run "<<runnr<<" et "<<et<<" eta "<<eta<<" r9 "<<r9<<" gain seed "<<gainSeed;
        std::cerr << "LogicError" << std::endl;
    }
    return &result.first->second;
}

void EnergyScaleCorrection::addScale(const std::string& category, int runMin, int runMax,
    double energyScale, double energyScaleErrStat,
    double energyScaleErrSyst, double energyScaleErrGain)
{

    CorrectionCategory cat(category, runMin, runMax); // build the category from the string

    // std::cout << energyScale << " " << energyScaleErrStat << " " << energyScaleErrSyst << " " << energyScaleErrGain << " " << std::endl;
    ScaleCorrection corr(energyScale, energyScaleErrStat, energyScaleErrSyst, energyScaleErrGain);
    scales_.push_back({ cat, corr });
    std::sort(scales_.begin(), scales_.end(), Sorter<CorrectionCategory, ScaleCorrection>());
}

// this method adds the correction values read from the txt file to the map
void EnergyScaleCorrection::addScale(int runMin_, int runMax_,
    float etaMin, float etaMax,
    float r9Min, float r9Max,
    float etMin, float etMax,
    unsigned int gain,
    double energyScale, double energyScaleErrStat,
    double energyScaleErrSyst, double energyScaleErrGain)
{
    //if(gain==0) gain=12;
    CorrectionCategory cat(runMin_, runMax_, etaMin, etaMax, r9Min, r9Max, etMin, etMax, gain); // build the category from the string

    // std::cout << "adding category " << cat << std::endl;
    auto result = std::equal_range(scales_.begin(), scales_.end(), cat, Sorter<CorrectionCategory, ScaleCorrection>());
    if (result.first != result.second) {
        //throw cms::Exception("ConfigError") << "Category already defined! "<<cat;
        std::cerr << "Category already defined" << std::endl;
    }
    ScaleCorrection corr(energyScale, energyScaleErrStat, energyScaleErrSyst, energyScaleErrGain);
    // scales[cat] = corr;
    scales_.push_back({ cat, corr });
    std::sort(scales_.begin(), scales_.end(), Sorter<CorrectionCategory, ScaleCorrection>());

    std::cout << "[INFO:scale correction] " << cat << "\t" << corr << std::endl;

    return;
}

void EnergyScaleCorrection::addSmearing(const std::string& category, int runMin, int runMax,
    double rho, double errRho,
    double phi, double errPhi,
    double eMean, double errEMean)
{
    CorrectionCategory cat(category);
    // std::cout << "add smearing category " << cat << std::endl;

    auto res = std::equal_range(smearings_.begin(), smearings_.end(), cat, Sorter<CorrectionCategory, SmearCorrection>());

    if (res.first != res.second) {
        //throw cms::Exception("EnergyScaleCorrection") << "Smearing category already defined "<<cat;
        std::cerr << "Smearing category already defined" << std::endl;
    }
    // std::cout << rho << "  "<<  errRho << "  "<< phi << "  "<< errPhi << "  "<< eMean << "  "<< errEMean << std::endl;
    SmearCorrection corr(rho, errRho, phi, errPhi, eMean, errEMean);
    smearings_.push_back({ cat, corr });
    // std::cout << "category " << cat << "   corr " << corr << std::endl;
    std::sort(smearings_.begin(), smearings_.end(), Sorter<CorrectionCategory, SmearCorrection>());
}

void EnergyScaleCorrection::setSmearingType(FileFormat value)
{
    if (value >= 0 && value <= 1) {
        smearingType_ = value;
    } else {
        smearingType_ = UNKNOWN;
    }
}

void EnergyScaleCorrection::readScalesFromFile(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file.good()) {
        //throw cms::Exception("EnergyScaleCorrection") << "file " << filename << " not readable.";
        std::cerr << "not readable" << std::endl;
    }

    int runMin, runMax;
    std::string category, region2;
    double energyScale, energyScaleErr, energyScaleErrStat, energyScaleErrSyst, energyScaleErrGain;
    float etaMin; ///< Min eta value for the bin
    float etaMax; ///< Max eta value for the bin
    float r9Min; ///< Min R9 vaule for the bin
    float r9Max; ///< Max R9 value for the bin
    float etMin; ///< Min Et value for the bin
    float etMax; ///< Max Et value for the bin
    unsigned int gain; ///< 12, 6, 1, 61 (double gain switch)

    // std::cout << "reading file" << std::endl;
    // for(file >> category; file.good(); file >> category) {

    // file >> region2
    // >> runMin >> runMax
    // >> energyScale >> energyScaleErr >> energyScaleErrStat >> energyScaleErrSyst >> energyScaleErrGain;
    // std::cout << category << " " << region2 << std::endl;
    // std::cout << energyScale << " " << energyScaleErr << " " << energyScaleErrStat << " " << energyScaleErrSyst << " " << energyScaleErrGain << std::endl;
    // addScale(category, runMin, runMax, energyScale, energyScaleErrStat, energyScaleErrSyst, energyScaleErrGain);
    // }

    if (smearingType_ == ECALELF) {
        for (file >> category; file.good(); file >> category) {
            file >> region2
                >> runMin >> runMax
                >> energyScale >> energyScaleErr >> energyScaleErrStat >> energyScaleErrSyst >> energyScaleErrGain;

            addScale(category, runMin, runMax, energyScale, energyScaleErrStat, energyScaleErrSyst, energyScaleErrGain);
        }
    } else {
        if (file.peek() == 'r')
            file.ignore(1000, 10);
        for (file >> runMin; file.good(); file >> runMin) {
            file >> runMax >> etaMin >> etaMax >> r9Min >> r9Max >> etMin >> etMax >> gain
                >> energyScale >> energyScaleErr;
            std::cout << energyScaleErr << " ##" << (char)file.peek() << "$$" << std::endl;
            file.ignore(1000, 10); //			if(file.peek()!=10) file>> err_deltaP_stat >> energyScaleErrSyst >> energyScaleErrGain;
            //else
            energyScaleErrStat = energyScaleErr;
            std::cout << runMin << "\t" << runMax << std::endl;
            addScale(runMin, runMax, etaMin, etaMax, r9Min, r9Max, etMin, etMax, gain, energyScale, energyScaleErrStat, energyScaleErrSyst, energyScaleErrGain);
        }
    }

    file.close();
    return;
}

//also more or less untouched function from the orginal package
void EnergyScaleCorrection::readSmearingsFromFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.good()) {
        //throw cms::Exception("EnergyScaleCorrection") << "file " << filename << " not readable";
        std::cerr << "not readable" << std::endl;
    }

    int runMin = 0;
    int runMax = 900000;
    int unused = 0;
    std::string category, region2;
    double rho, phi, eMean, errRho, errPhi, errEMean;
    double etaMin, etaMax, r9Min, r9Max;
    std::string phiString, errPhiString;

    while (file.peek() != EOF && file.good()) {
        if (file.peek() == 10) { // 10 = \n
            file.get();
            continue;
        }

        if (file.peek() == 35) { // 35 = #
            file.ignore(1000, 10); // ignore the rest of the line until \n
            // std::cout << "hello " << std::endl;
            continue;
        }

        if (smearingType_ == UNKNOWN) { // trying to guess: not recommended
            //throw cms::Exception("ConfigError") <<"unknown smearing type";
            // std::cerr << "unknown smearing type" << std::endl;

        } else if (smearingType_ == GLOBE) {
            file >> category >> unused >> etaMin >> etaMax >> r9Min >> r9Max >> runMin >> runMax >> eMean >> errEMean >> rho >> errRho >> phi >> errPhi;

            addSmearing(category, runMin, runMax, rho, errRho, phi, errPhi, eMean, errEMean);

        } else if (smearingType_ == ECALELF || smearingType_ == TABLE) {
            // std::cout << "hello TABLE" << std::endl;
            file >> category >> eMean >> errEMean >> rho >> errRho >> phiString >> errPhiString;
            // std::cout << "emean " << eMean << "  errEMean " << errEMean << "  ro " << rho << std::endl;
            // std::cout << "category name " << category << std::endl;
            if (phiString == "M_PI_2")
                phi = M_PI_2;
            else
                phi = std::stod(phiString);

            if (errPhiString == "M_PI_2")
                errPhi = M_PI_2;
            else
                errPhi = std::stod(errPhiString);
            // std::cout << "phi " << phi << "   errPhi " << errPhi << std::endl;

            // std::cout << runMin << "  "<< runMax << "  "<< rho << "  "<<  errRho << "  "<< phi << "  "<< errPhi << "  "<< eMean << "  "<< errEMean << std::endl;
            addSmearing(category, runMin, runMax, rho, errRho, phi, errPhi, eMean, errEMean);

        } else {
            file >> category >> rho >> phi;
            errRho = errPhi = eMean = errEMean = 0;
            addSmearing(category, runMin, runMax, rho, errRho, phi, errPhi, eMean, errEMean);
        }
    }

    file.close();
    return;
}

std::ostream& EnergyScaleCorrection::ScaleCorrection::print(std::ostream& os) const
{
    os << "( " << scale_ << " +/- " << scaleErrStat_ << " +/- " << scaleErrSyst_ << " +/- " << scaleErrGain_ << ")";
    return os;
}

// float EnergyScaleCorrection::ScaleCorrection::scaleErr(const std::bitset<kErrNrBits>& uncBitMask)const
float EnergyScaleCorrection::ScaleCorrection::scaleErr(int type) const
{
    double totErr(0);
    auto pow2 = [](const double& x) { return x * x; };
    // std::cout << "bitmask " << type << std::endl;
    // std::cout << scaleErrStat_ << " " << scaleErrSyst_ << " " << scaleErrGain_ << std::endl;
    if (type == kErrStatBitNr)
        totErr += pow2(scaleErrStat_);
    if (type == kErrSystBitNr)
        totErr += pow2(scaleErrSyst_);
    if (type == kErrGainBitNr)
        totErr += pow2(scaleErrGain_);
    // if(uncBitMask.test(kErrStatBitNr)) totErr+=pow2(scaleErrStat_);
    // if(uncBitMask.test(kErrSystBitNr)) totErr+=pow2(scaleErrSyst_);
    // if(uncBitMask.test(kErrGainBitNr)) totErr+=pow2(scaleErrGain_);

    return std::sqrt(totErr);
}

void EnergyScaleCorrection::SmearCorrection::print() const
{
    std::cout << rho_ << " +/- " << rhoErr_
              << "\t"
              << phi_ << " +/- " << phiErr_
              << "\t"
              << eMean_ << " +/- " << eMeanErr_;
    return;
}

//here be dragons
//this function is nasty and needs to be replaced
EnergyScaleCorrection::CorrectionCategory::CorrectionCategory(const std::string& category, int runnrMin, int runnrMax)
    : runMin_(runnrMin)
    , runMax_(runnrMax)
    , etaMin_(0)
    , etaMax_(3)
    , r9Min_(-1)
    , r9Max_(999)
    , etMin_(0)
    , etMax_(9999999)
    , gain_(0)
{
    size_t p1, p2; // boundary

    // std::cout << "before: " << r9Max_ << " " << etaMax_ << " " << etMax_ << " " << gain_ << std::endl;
    // eta region
    p1 = category.find("absEta_");

    if (category.find("absEta_0_1.4442") != std::string::npos) {
        // cout << "cat 0a" << endl;
        etaMin_ = 0;
        etaMax_ = 1.479;
    } else if (category.find("absEta_1.566_2.5") != std::string::npos) {

        // cout << "cat 0b" << endl;
        etaMin_ = 1.479;
        etaMax_ = 3;
    } else if (category.find("absEta_0_1") != std::string::npos) {
        // cout << "cat a" << endl;
        etaMin_ = 0;
        etaMax_ = 1;
    } else if (category.find("absEta_1_1.4442") != std::string::npos) {

        // cout << "cat b" << endl;
        etaMin_ = 1;
        etaMax_ = 1.479;
    } else if (category.find("absEta_1.566_2") != std::string::npos) {
        // cout << "cat c" << endl;
        etaMin_ = 1.479;
        etaMax_ = 2;
    } else if (category.find("absEta_2_2.5") != std::string::npos) {
        // cout << "cat d" << endl;
        etaMin_ = 2;
        etaMax_ = 3;
    } else {
        if (p1 != std::string::npos) {
            // std::cout << "hello" << std::endl;
            p1 = category.find("_", p1);
            p2 = category.find("_", p1 + 1);
            etaMin_ = std::stof(category.substr(p1 + 1, p2 - p1 - 1));
            p1 = p2;
            p2 = category.find("-", p1);
            etaMax_ = std::stof(category.substr(p1 + 1, p2 - p1 - 1));
        }
    }
    // cout << "eta min " << etaMin_ << "  eta Max " << etaMax_ << endl;

    if (category.find("EBlowEta") != std::string::npos) {
        etaMin_ = 0;
        etaMax_ = 1;
    };
    if (category.find("EBhighEta") != std::string::npos) {
        etaMin_ = 1;
        etaMax_ = 1.479;
    };
    if (category.find("EElowEta") != std::string::npos) {
        etaMin_ = 1.479;
        etaMax_ = 2;
    };
    if (category.find("EEhighEta") != std::string::npos) {
        etaMin_ = 2;
        etaMax_ = 7;
    };

    // Et region
    p1 = category.find("-Et_");

    if (p1 != std::string::npos) {
        p1 = category.find("_", p1);
        p2 = category.find("_", p1 + 1);
        etMin_ = std::stof(category.substr(p1 + 1, p2 - p1 - 1));
        p1 = p2;
        p2 = category.find("-", p1);
        etMax_ = std::stof(category.substr(p1 + 1, p2 - p1 - 1));
    }

    if (category.find("gold") != std::string::npos || category.find("Gold") != std::string::npos || category.find("highR9") != std::string::npos) {
        r9Min_ = 0.94;
        r9Max_ = std::numeric_limits<float>::max();
    } else if (category.find("bad") != std::string::npos || category.find("Bad") != std::string::npos || category.find("lowR9") != std::string::npos) {
        r9Min_ = -1;
        r9Max_ = 0.94;
    };
    // R9 region
    p1 = category.find("-R9");
    p2 = p1 + 1;
    if (p1 != std::string::npos) {
        p1 = category.find("_", p1);
        p2 = category.find("_", p1 + 1);
        r9Min_ = std::stof(category.substr(p1 + 1, p2 - p1 - 1));
        // If there is one value, just set lower bound
        if (p2 != std::string::npos) {
            p1 = p2;
            p2 = category.find("-", p1);
            r9Max_ = std::stof(category.substr(p1 + 1, p2 - p1 - 1));
            if (r9Max_ >= 1.0)
                r9Max_ = std::numeric_limits<float>::max();
        }
    }
    // cout << "r9 min " << r9Min_ << "  r9 max " << r9Max_ << endl;
    //------------------------------
    p1 = category.find("gainEle_"); // Position of first character
    if (p1 != std::string::npos) {
        p1 += 8; // Position of character after _
        p2 = category.find("-", p1); // Position of - or end of string
        gain_ = std::stoul(category.substr(p1, p2 - p1), nullptr);
    }

    // std::cout << "corr: " << r9Max_ << " " << etaMax_ << " " << etMax_ << " " << gain_ << std::endl;
    //so turns out the code does an includes X<=Y<=Z search for bins
    //which is what we want for run numbers
    //however then the problem is when we get a value exactly at the bin boundary
    //for the et/eta/r9 which then gives multiple bins
    //so we just decrement the maxValues ever so slightly to ensure that they are different
    //from the next bins min value
    etMax_ = std::nextafterf(etMax_, std::numeric_limits<float>::min());
    etaMax_ = std::nextafterf(etaMax_, std::numeric_limits<float>::min());
    r9Max_ = std::nextafterf(r9Max_, std::numeric_limits<float>::min());
}
bool EnergyScaleCorrection::CorrectionCategory::
    inCategory(const unsigned int runnr, const float et, const float eta, const float r9,
        const unsigned int gainSeed) const
{
    // std::cout << "check values" << std::endl;
    // std::cout << runnr << "  "<< et << "  "<< eta << "  "<< r9 << "  "<< gainSeed << std::endl;
    // std::cout << runMin_ << "  "<< runMax_ << "  "<< etMin_ << "  "<<  etMax_ << "  "<< etaMin_ << "  "<< etaMax_ << "  "<< r9Min_ << "  "<< r9Max_ << "  " << gain_ << std::endl;
    return runnr >= runMin_ && runnr <= runMax_ && et >= etMin_ && et <= etMax_ && eta >= etaMin_ && eta <= etaMax_ && r9 >= r9Min_ && r9 <= r9Max_ && (gain_ == 0 || gainSeed == gain_);
}

bool EnergyScaleCorrection::CorrectionCategory::operator<(const EnergyScaleCorrection::CorrectionCategory& b) const
{
    if (runMin_ < b.runMin_ && runMax_ < b.runMax_)
        return true;
    if (runMax_ > b.runMax_ && runMin_ > b.runMin_)
        return false;

    if (etaMin_ < b.etaMin_ && etaMax_ < b.etaMax_)
        return true;
    if (etaMax_ > b.etaMax_ && etaMin_ > b.etaMin_)
        return false;

    if (r9Min_ < b.r9Min_ && r9Max_ < b.r9Max_)
        return true;
    if (r9Max_ > b.r9Max_ && r9Min_ > b.r9Min_)
        return false;

    if (etMin_ < b.etMin_ && etMax_ < b.etMax_)
        return true;
    if (etMax_ > b.etMax_ && etMin_ > b.etMin_)
        return false;

    if (gain_ == 0 || b.gain_ == 0)
        return false; // if corrections are not categorized in gain then default gain value should always return false in order to have a match with the category
    if (gain_ < b.gain_)
        return true;
    else
        return false;
    return false;
}

std::ostream& EnergyScaleCorrection::CorrectionCategory::print(std::ostream& os) const
{
    os << runMin_ << " " << runMax_
       << "\t" << etaMin_ << " " << etaMax_
       << "\t" << r9Min_ << " " << r9Max_
       << "\t" << etMin_ << " " << etMax_
       << "\t" << gain_;
    return os;
}
