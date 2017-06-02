from ._solvers import *

# FIXME all parameters are mor or less wild guesses...
def multicutMp( uvIds,
        costs,
        primalComputationInterval = 50,
        standardReparametrization = "anisotropic",
        roundingReparametrization = "damped_uniform",
        tightenReparametrization  = "damped_uniform",
        tighten                   = True,
        tightenInterval           = 50,
        tightenIteration          = 5,
        tightenSlope              = 0.1,
        tightenConstraintsPercentage = 0.05,
        maxIter                   = 1000,
        minDualImprovement        = 0,
        minDualImprovementInterval= 0,
        timeout                   = 0,
        nThreads                  = 1
    ):

    mc_opts = MulticutOptions()
    mc_opts.primalComputationInterval = primalComputationInterval
    mc_opts.standardReparametrization = standardReparametrization
    mc_opts.roundingReparametrization = roundingReparametrization
    mc_opts.tightenReparametrization  = tightenReparametrization
    mc_opts.tighten                   = tighten
    mc_opts.tightenInterval           = tightenInterval
    mc_opts.tightenIteration          = tightenIteration
    mc_opts.tightenSlope              = tightenSlope
    mc_opts.tightenConstraintsPercentage = tightenConstraintsPercentage
    mc_opts.maxIter                   = maxIter
    mc_opts.minDualImprovement        = minDualImprovement
    mc_opts.minDualImprovementInterval= minDualImprovementInterval
    mc_opts.timeout                   = timeout
    mc_opts.nThreads                  = nThreads

    return multicutImpl(uvIds, costs, mc_opts)
