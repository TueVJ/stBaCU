import countries

#Wind production in units of average load
def get_Wind(D):
    return D[:,0] * D[:,4] / D[:,0].mean()

#Solar production in units of average load
def get_PV(D):
    if D[:,1].sum()==0:
        return D[:,1]
    return D[:,1] * D[:,4] / D[:,1].mean()

#Load in units of average load
def get_Load(D):
    return D[:,2] * D[:,4] / D[:,2].mean()

# Mismatch for the specified gamma (penetration fraction)
# and a (Wind/Solar production)
def get_mismatch(gamma = 1., a = .6):
    D = countries.Region
    Wind = get_Wind(D = D)
    PV = get_PV(D = D)
    Load = get_Load(D = D)
    return gamma * (a * Wind + (1. - a) * PV) - Load

