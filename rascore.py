from RAscore import RAscore_NN

nn_scorer = RAscore_NN.RAScorerNN()

#  Imatinib mesylate
imatinib_mesylate = 'CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5.CS(=O)(=O)O'
score = nn_scorer.predict(imatinib_mesylate)
print(f"RAscore of imatinib mesylat: {score}")
