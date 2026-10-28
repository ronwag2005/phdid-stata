*! version 0.1.0  phdid_mata -- put the Mata library in the session
*! Rohan Wagle and Parush Arora
*
*  The numerical core of phdid is written in Mata and ships as source,
*  in src/phdid_{core,drivers,bayes,post}.mata.  It is not distributed
*  as a byte-compiled library (.mlib) because a library written by one
*  Stata cannot be read by an older one, which would silently break
*  every command for anyone running a Stata older than the build
*  machine's.  Compiling the sources instead takes about a hundredth of
*  a second and works on any Stata from 17 on.
*
*  Every phdid command calls this first.  It returns at once when the
*  functions are already compiled, so the cost is paid once per session.

program define phdid_mata
    version 17.0

    * already in memory?  nothing to do
    capture mata: phdid_abi()
    if (_rc == 0) {
        exit
    }

    * the sources compile themselves under matastrict; leave the user's
    * own Mata settings as we found them
    local ms "`c(matastrict)'"

    foreach f in core drivers bayes post {
        capture findfile phdid_`f'.mata
        if (_rc) {
            di as err "phdid: phdid_`f'.mata not found on the adopath"
            di as txt "{p 4 4 2}The installation is incomplete.  " ///
                      "Reinstall with {bf:ssc install phdid, replace}, " ///
                      "or see {help phdid}.{p_end}"
            exit 601
        }
        qui do "`r(fn)'"
    }

    if ("`ms'" != "") {
        qui mata: mata set matastrict `ms'
    }

    capture mata: phdid_abi()
    if (_rc) {
        di as err "phdid: the Mata sources did not compile"
        exit 3499
    }
end
