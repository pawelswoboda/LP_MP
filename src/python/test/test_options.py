# FIXME dirty hack
import sys
sys.path.append('/home/consti/Work/software/bld/LP_MP/python')
# this assumes that lp_mp is in the python path
import lp_mp

def test_mc_opts():
    mc_opts = lp_mp.solvers.MulticutOptions()

    assert mc_opts.primalComputationInterval == 100
    assert mc_opts.standardReparametrization == "anisotropic"
    assert mc_opts.roundingReparametrization == "damped_uniform"
    assert mc_opts.tightenReparametrization == "damped_uniform"
    assert mc_opts.tighten == True
    assert mc_opts.tightenInterval == 100
    assert mc_opts.tightenIteration == 2
    assert mc_opts.tightenSlope == 0.05
    assert mc_opts.tightenConstraintsPercentage == 0.1

    print "Passed"


if __name__ == '__main__':
    test_mc_opts()
