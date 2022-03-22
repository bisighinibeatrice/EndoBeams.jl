
function init(solvertype)

    if solvertype == :MKL

        ps = MKLPardisoSolver()
        set_nprocs!(ps, Threads.nthreads())
        set_matrixtype!(ps, Pardiso.REAL_NONSYM)
        pardisoinit(ps)

        set_iparm!(ps, 1, 1) # set to allow non-default iparams
        set_iparm!(ps, 12, 2) # CSC matrix format

        return ps

    end

end


function analyze!(ps::Pardiso.MKLPardisoSolver, K, r)

    set_phase!(ps, Pardiso.ANALYSIS)
    pardiso(ps, K, r)

end

function solve!(ps::Pardiso.MKLPardisoSolver, x, K, r)

    set_phase!(ps, Pardiso.NUM_FACT_SOLVE_REFINE)
    pardiso(ps, x, K, r)

end

function release!(ps::Pardiso.MKLPardisoSolver)

    set_phase!(ps, Pardiso.RELEASE_ALL)
    pardiso(ps)

end



