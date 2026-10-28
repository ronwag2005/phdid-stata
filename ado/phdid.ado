*! version 0.1.0  phdid -- partial homogeneity in staggered DiD
*! Rohan Wagle and Parush Arora
program define phdid
    version 17.0
    gettoken sub 0 : 0, parse(" ,")

    if ("`sub'" == "") {
        di as err "nothing to do"
        di as txt "{p 4 4 2}Usage: {bf:phdid setup}, then {bf:phdid l0}, " ///
                  "{bf:phdid bayes} or {bf:phdid test}. See {help phdid}.{p_end}"
        exit 198
    }

    if ("`sub'" == "setup") {
        phdid_setup `0'
        exit
    }
    if ("`sub'" == "describe") {
        phdid_desc `0'
        exit
    }
    if ("`sub'" == "fit") {
        phdid_fit `0'
        exit
    }
    if ("`sub'" == "l0") {
        phdid_l0 `0'
        exit
    }
    if ("`sub'" == "test") {
        phdid_test `0'
        exit
    }
    if ("`sub'" == "bayes") {
        phdid_bayes `0'
        exit
    }
    if ("`sub'" == "enumerate") {
        phdid_enum `0'
        exit
    }
    if ("`sub'" == "aggregate") {
        phdid_agg `0'
        exit
    }
    if ("`sub'" == "rhat") {
        phdid_rhat `0'
        exit
    }
    if ("`sub'" == "partition") {
        phdid_point `0'
        exit
    }
    if ("`sub'" == "covcheck") {
        phdid_cov `0'
        exit
    }
    if ("`sub'" == "sensitivity") {
        phdid_sens `0'
        exit
    }
    if ("`sub'" == "graph") {
        phdid_graph `0'
        exit
    }

    di as err "unknown subcommand {bf:`sub'}"
    di as txt "{p 4 4 2}Valid subcommands are {bf:setup}, {bf:describe}, " ///
              "{bf:l0}, {bf:bayes}, {bf:fit}, {bf:test}, " ///
              "{bf:aggregate}, " ///
              "{bf:partition}, {bf:rhat}, {bf:sensitivity}, " ///
              "{bf:covcheck}, {bf:enumerate} and {bf:graph}.{p_end}"
    exit 198
end
