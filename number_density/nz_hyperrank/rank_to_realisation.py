import numpy as np
import sys

def rank_to_realisation(chain_file, rank_file):

    ranks = np.loadtxt(rank_file)
    chain = np.genfromtxt(chain_file, names = True)

    rank_chain = chain['ranksrank_hyperparm_1']

    n_realisations = np.max(ranks)+1

    realisations = np.argwhere([ranks == np.floor(r * n_realisations) for r in rank_chain])[:,1]

    return realisations

if __name__ == "__main__":

    chain_file = sys.argv[1]
    rank_file = sys.argv[2]

    realisations = rank_to_realisation(chain_file, rank_file)

    np.savetxt('realisations_' + chain_file, realisations, fmt='%d')
