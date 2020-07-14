import ROOT
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

def findWP(fFake = "test.root", hFake = "highestDisc", sampleFake = "histsTTbar2l", rate = .1):
    f_fake = ROOT.TFile.Open(fFake)

    hdisc_fake = f_fake.Get(sampleFake + "/" + hFake)

    hdisc_fake.Print()

    fakeInt = hdisc_fake.Integral(0, hdisc_fake.GetNbinsX() + 1)
    nBins = hdisc_fake.GetNbinsX()

    closestDisc = -1
    error = 100

    for i in xrange(1, nBins + 2):
        if abs(hdisc_fake.Integral(i, nBins + 1)/fakeInt - rate) < error:
            error = abs(hdisc_fake.Integral(i, nBins + 1)/fakeInt - rate)
            closestDisc = i

#            print hdisc_fake.GetBinCenter(i), error

    r1 = hdisc_fake.Integral(closestDisc-1, nBins + 1)/fakeInt
    d1 = hdisc_fake.GetBinCenter(closestDisc-1)

    r2 = hdisc_fake.Integral(closestDisc, nBins + 1)/fakeInt
    d2 = hdisc_fake.GetBinCenter(closestDisc)

    r3 = hdisc_fake.Integral(closestDisc+1, nBins + 1)/fakeInt
    d3 = hdisc_fake.GetBinCenter(closestDisc+1)

    M = np.array([[d1*d1,d1,1],[d2*d2,d2,1],[d3,d3,1]])
    Minv = np.linalg.inv(M)

#    print M
#    print Minv
#
#    print np.dot(M,Minv)
#    print np.dot(Minv,[[r1],[r2],[r3]])

    a = np.dot(Minv,[[r1],[r2],[r3]])[0]
    b = np.dot(Minv,[[r1],[r2],[r3]])[1]
    c = np.dot(Minv,[[r1],[r2],[r3]])[2]

#    print a*d1*d1+b*d1+c-r1
#    print a*d2*d2+b*d2+c-r2
#    print a*d3*d3+b*d3+c-r3

    f_fake.Close()
    return round(abs((rate-c-a*d2*d2)/b),3)



loose = findWP(fFake = "QCD.root", hFake = "highestDisc", sampleFake = "histsTTbarNol", rate = .1)
medium = findWP(fFake = "QCD.root", hFake = "highestDisc", sampleFake = "histsTTbarNol", rate = .05)
tight = findWP(fFake = "QCD.root", hFake = "highestDisc", sampleFake = "histsTTbarNol", rate = .01)

print "loose: ",loose
print "medium: ",medium
print "tight: ",tight

looseQCD = findWP(fFake = "QCD.root", hFake = "highestDisc", sampleFake = "QCD", rate = .1)
mediumQCD = findWP(fFake = "QCD.root", hFake = "highestDisc", sampleFake = "QCD", rate = .05)
tightQCD = findWP(fFake = "QCD.root", hFake = "highestDisc", sampleFake = "QCD", rate = .01)

print "loose: ",looseQCD
print "medium: ",mediumQCD
print "tight: ",tightQCD
