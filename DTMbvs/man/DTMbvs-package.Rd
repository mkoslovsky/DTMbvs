\name{DTMbvs-package}
\alias{DTMbvs-package}
\alias{DTMbvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DTMbvs 
}
\description{
Performs Bayesian variable selection for multivariate count data with tree-like structure. Accommodates various prior distributions for inclusion indicators, including beta-binomial and Markov random field. }
\usage{
DTMbvs_R()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iterations}{integer number of MCMC samples, default = 50000}
  \item{thin}{integer MCMC by # thin, default = 10}
  \item{tree}{tree object .tre file to be sourced}
  \item{Y}{N x branch matrix of count data. Make sure branch order corresponds to edge list in tree}
  \item{X}{N x P matrix of covariates}
  \item{prior}{name of selection prior to use. Takes values "BB" beta-binomial (default), "MRF_fixed" known Markov random field, "MRF_unknown" unknown Markov random field}
  \item{alpha}{branch x 1 vector of branch intercepts initialization}
  \item{phi}{branch x P matrix of regression coefficients initialization}
  \item{zeta}{branch x P matrix of inclusion indicators initialization}
  \item{sigma2_alpha}{prior variance for alpha, default = 10}
  \item{sigma2_phi}{prior variance for phi, default = 10}
  \item{aa}{double beta-binomial hyperparameter}
  \item{bb}{double beta-binomial hyperparameter}
  \item{a_G}{double MRF hyperparameter}
  \item{b_G}{double MRF hyperparameter}
  \item{Omega}{P x P concentration matrix initialization}
  \item{G}{P x P adjacency matrix for relational graph initialization}
  \item{v0}{double variance of exclusion for the relational graph, default = 0.1}
  \item{v1}{double variance of inclusion for the relational graph, default = 10}
  \item{pie}{double prior probability of inclusion for the relational graph}
  \item{lambda}{double exponential prior hyperparameter}
  \item{subject_sim}{integer number of subjects to simulate, default = 100}
  \item{num_leaf}{integer number of leaves in tree to simulate, default = 5}
  \item{covariates_sim}{integer number of covariates to simulate, default = 50}
  \item{corr}{double correlation structure, default = 0.20}
  \item{num_branch}{integer number of branches the associated covariates are in, default = 3}
  \item{num_cov}{integer number of associated covariates to simulate, default = 5}
  \item{seed}{integer seed #, default = 1212}
  \item{eval}{boolean if true, run model on simulated data. Otherwise, simply simulate}
  \item{warmstart}{boolean if true, start model using informed initial values. }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Koslovsky, M.D., and Vannucci, M. DTMbvs: Dirichlet-tree Multinomial Regression Models with Bayesian Variable Selection for Microbiome Data - an R Package 2019+.
}
\author{
Matthew Koslovsky email: mkoslovsky12@gmail.com
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
run <- DTMbvs_R()
output <- selected( run, plotting = TRUE )

