# FIXME dirty hack
import sys
sys.path.append('/home/consti/Work/software/bld/LP_MP/python')
# this assumes that lp_mp is in the python path
import lp_mp

def simple_multicut():
    uv_ids = [(0,1), (1,2), (2,3), (3,4), (4,0)]
    weights = [-1. , -1.  , -1.  , -1.  , -1   ]

    mc_opts = lp_mp.solvers.MulticutOptions()
    print mc_opts.vec()
    edges = lp_mp.solvers.multicut(uv_ids, weights, mc_opts)

    assert len(edges) == len(weights)
    print edges


if __name__ == '__main__':
    simple_multicut()
