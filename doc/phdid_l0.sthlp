{smcl}
{* *! version 0.1.0  26oct2026}{...}
{vieweralsosee "phdid" "help phdid"}{...}
{vieweralsosee "phdid bayes" "help phdid_bayes"}{...}
{vieweralsosee "phdid test" "help phdid_test"}{...}
{vieweralsosee "phdid postestimation" "help phdid_postestimation"}{...}
{viewerjumpto "Syntax" "phdid_l0##syntax"}{...}
{viewerjumpto "Description" "phdid_l0##description"}{...}
{viewerjumpto "Options" "phdid_l0##options"}{...}
{viewerjumpto "Remarks" "phdid_l0##remarks"}{...}
{viewerjumpto "Examples" "phdid_l0##examples"}{...}
{viewerjumpto "Stored results" "phdid_l0##results"}{...}
{title:Title}

{phang}
{bf:phdid l0} {hline 2} The l0-penalised partial-homogeneity estimator


{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:phdid l0} [{cmd:,} {it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Selection}
{synopt:{opt sel:ect(rule)}}{opt bic} (default), {opt lambda} or {opt m}{p_end}
{synopt:{opt l:ambda(#)}}penalty, with {cmd:select(lambda)}{p_end}
{synopt:{opt m(#)}}number of groups, with {cmd:select(m)}{p_end}
{synopt:{opt nbic(#)}}sample size in the BIC penalty; default {it:K}{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt nopath}}suppress the agglomeration path table{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
{cmd:phdid l0} is {cmd:eclass}.  {helpb phdid_setup:phdid setup} must be run
first.


{marker description}{...}
{title:Description}

{pstd}
{cmd:phdid l0} computes the two-step l0-penalised estimator of Appendix B of
Arora and Wagle (2026).

{pstd}
{bf:Step 1} is a greedy agglomerative search from {it:K} singleton groups.  At
each stage it merges the pair minimising

{p 12 12 2}
{it:n_A n_B} / ({it:n_A} + {it:n_B}) ({it:tau_A} - {it:tau_B})^2,

{pstd}
where {it:n} is a precision weight and {it:tau} a group mean, giving a full
path of nested partitions from {it:K} groups down to one.

{pstd}
{bf:Step 2} re-estimates by restricted generalised least squares under the
selected partition.  This is not a refinement.  The merge criterion in step 1
treats the design as orthogonal, and only the re-estimated fit is best linear
unbiased in a panel where the cells are correlated.  Step 2 is always
performed.

{pstd}
The estimator returns one partition.  Its standard errors condition on that
partition being the true one, which is the price of a single interpretable
answer; {helpb phdid_bayes:phdid bayes} marginalises over the partition
instead.


{marker options}{...}
{title:Options}

{phang}
{opt select(rule)} chooses the point on the path.

{pmore}
{opt bic}, the default, minimises {it:deviance} + {it:m} ln({it:n}) over the
path.

{pmore}
{opt lambda} applies the Appendix B stopping rule: merge while the merge cost
is below {it:lambda} |{it:A}| |{it:B}|.  Requires {opt lambda()}.

{pmore}
{opt m} takes the partition with exactly that many groups.  Requires
{opt m()}.

{phang}
{opt lambda(#)} is the l0 penalty.  Larger values pool more.

{phang}
{opt m(#)} fixes the number of groups, between 1 and {it:K}.

{phang}
{opt nbic(#)} sets the sample size {it:n} in the BIC penalty.  The theory does
not fix it: the {it:K} first-stage estimates are the observations, while the
asymptotics justifying Schwarz's rate run in the number of units.  The default
is {it:K}, which reproduces the group counts reported in the paper.  The
option is exposed rather than buried because it changes the selected
partition.

{phang}
{opt level(#)} sets the confidence level for the reported intervals.

{phang}
{opt nopath} suppresses the table of the agglomeration path.


{marker remarks}{...}
{title:Remarks}

{dlgtab:Reading the output}

{pstd}
The per-cell table gives the flexible first-stage estimate, the group the cell
was assigned to, the grouped estimate and its standard error, and a variance
ratio: the grouped sampling variance relative to the flexible one.  A ratio
below one means pooling bought precision for that cell, and in a balanced
design the ratio is 1/|{it:C}| for a cell in a group of size |{it:C}|.  The
ratio can never exceed one.

{dlgtab:Lambda is not monotone in the number of groups}

{pstd}
The Appendix B stopping rule compares each merge cost against {it:lambda}
|{it:A}| |{it:B}|, and that pair product grows as groups absorb one another.
The effective per-merge threshold is therefore not monotone in the merge
order, and some group counts are unreachable by any single {it:lambda}.  On the
paper's first application the thresholds run 0.079, 0.686, 0.397, 1.350,
2.765, 3.571, so no value of {it:lambda} selects five groups.  To walk the path
one merge at a time instead, use {cmd:phdid sensitivity, by(m)}.  The
{cmd:pairprod} column of {cmd:e(path)} reports |{it:A}| |{it:B}| at each merge.

{dlgtab:The deviance is not a test of the selected partition}

{pstd}
The deviance reported here was minimised over the path, so referring it to a
chi-squared distribution would be invalid.  To test a partition, fix it in
advance and use {cmd:phdid fit}, which reports the same statistic as an honest
test.


{marker examples}{...}
{title:Examples}

{phang2}{cmd:. phdid l0}{p_end}
{phang2}{cmd:. phdid aggregate}{p_end}
{phang2}{cmd:. phdid graph groups}{p_end}

{pstd}A fixed number of groups, and a fixed penalty:{p_end}

{phang2}{cmd:. phdid l0, select(m) m(3)}{p_end}
{phang2}{cmd:. phdid l0, select(lambda) lambda(0.5)}{p_end}

{pstd}What the BIC sample size does to the answer:{p_end}

{phang2}{cmd:. phdid l0, nbic(2500) nopath}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:phdid l0} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(m)}}number of groups selected{p_end}
{synopt:{cmd:e(K)}}number of cells{p_end}
{synopt:{cmd:e(deviance)}}deviance of the selected fit{p_end}
{synopt:{cmd:e(bic)}}its BIC{p_end}
{synopt:{cmd:e(nbic)}}sample size used in the penalty{p_end}
{synopt:{cmd:e(att)}}overall ATT{p_end}
{synopt:{cmd:e(att_se)}}its standard error{p_end}
{synopt:{cmd:e(lambda)}}penalty, with {cmd:select(lambda)}{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:phdid_l0}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(select)}}selection rule used{p_end}
{synopt:{cmd:e(partition_str)}}the partition, written out{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}grouped cell effects{p_end}
{synopt:{cmd:e(V)}}their restricted GLS covariance{p_end}
{synopt:{cmd:e(cells)}}flexible and grouped estimates, with variance
ratios{p_end}
{synopt:{cmd:e(path)}}deviance, BIC, merge cost and pair product by
{it:m}{p_end}
{synopt:{cmd:e(phi)}}the {it:m} distinct group effects{p_end}
{synopt:{cmd:e(partition)}}group number of each cell{p_end}
{p2colreset}{...}
