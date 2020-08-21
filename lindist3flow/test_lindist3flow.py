from lindist3flow import LinDist3Flow


dss_path = './data/compare_opendss_05node_threephase_unbalanced_oscillation_03'

slackidx = 0
Vslack = np.array([1, np.exp(1j*-120*np.pi/180), np.exp(1j*120*np.pi/180)])

LinDist3Flow(dss_path=dss_path, slackidx, Vslack)