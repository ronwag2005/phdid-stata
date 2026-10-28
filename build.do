*! version 0.1.0  phdid: build the Mata library
*! Rohan Wagle and Parush Arora
*
*  Compiles the Mata source into src/lphdid.mlib, which is what the ado
*  files call at run time.  Run this after changing anything under src/.
*
*      do build.do
*
*  The library is byte-compiled, so an installed copy of phdid does not
*  need the .mata sources; they ship for inspection and for certification.

version 17.0
clear all
set more off

local here "`c(pwd)'"

mata: mata clear
mata: mata set matastrict on
foreach f in core drivers bayes post {
    qui do "`here'/src/phdid_`f'.mata"
}

mata: mata mlib create lphdid, dir("`here'/src") replace
mata: mata mlib add lphdid phdid_*() mm_median_phdid(), dir("`here'/src")
adopath ++ "`here'/src"
mata: mata mlib index

di as txt _n "built src/lphdid.mlib"
mata: mata describe using lphdid
