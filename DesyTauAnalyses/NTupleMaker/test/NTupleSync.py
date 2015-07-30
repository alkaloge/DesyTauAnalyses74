from DesyTauAnalyses.NTupleMaker.CompareSpring15 import areEqual

print areEqual(0.,0.0001)
print areEqual(3, 3)
print areEqual(3, 3.000001)
print areEqual(3, 3.00001)

print areEqual(0.00005, -0.0000000005)
