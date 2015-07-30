from math import fabs, log10, sin, cos, acos
from ROOT import TFile, TChain

def areEqual( x1, x2):
    if( x1 == 0. or x2 == 0.):
        if( fabs(x1-x2) > 1.e-5):
            return False
    else:
        if (x1 * x2 < 0.):
            return False
        else:
            if (fabs(log10(fabs(x1))-log10(fabs(x2))) > 1.e-5):
                return False
    return True

def dPhiFrom2P( Px1, Py1, Px2, Py2):
    prod = Px1*Px2 + Py1*Py2

    mod1 = sqrt(Px1*Px1+Py1*Py1)
    mod2 = sqrt(Px2*Px2+Py2*Py2)
    
    cosDPhi = prod/(mod1*mod2)
    
    return acos(cosDPhi)


def deltaR( Eta1, Phi1,	Eta2, Phi2):
    Px1 = cos(Phi1)
    Py1 = sin(Phi1)

    Px2 = cos(Phi2)
    Py2 = sin(Phi2)

    dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2)
    dEta = Eta1 - Eta2

    return sqrt(dPhi*dPhi+dEta*dEta)


class CompareSpring15:
    def __init__(self, ref, test):
        self.cref=TChain(ref.split(":")[1])
        self.ctest=TChain(test.split(":")[1])

        self.cref.Add(ref.split(":")[0])
        self.ctest.Add(test.split(":")[0])
        
    def Compare(self):
        print self.cref.GetEntries(), self.ctest.GetEntries()
        


        
