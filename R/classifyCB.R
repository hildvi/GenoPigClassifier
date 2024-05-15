#' Function for classifying cross breed combination based on genomic data (GeneSeek 50K (Illumina) SNP chip).

#' @param DNA: matrix (n x m) with 0-1-2 coding of SNP's to be evaluated, typically data returned from ped2num()
#' @param alpha0: Deafault value is 73.58105, i.e. the free parameter \eqn{\alpha_0} used 
#' in calculation of \eqn{V(\theta)}.
#' @param Unkn_LogLike: the log likelihood for the uniform distribution for unknown breed in PLS-QDA. 
#' Defaults value is NULL, for which the value is calculated based on the m dimensional
#' space spanned bu the scores in the PLS-model for PB's.
#' @param PriorDist: "NULL" (default) for which a totally flat prior is applied, 
#' 'informative' for which the informative prior in Vinje et.al.#' is applied, 
#' or a scalar >0 and <1 defining the prior probability for class "Unknown" 
#' and a flat prior is applied to 
#' the rest of the classes, 
#' @param pred_min: Parameter to be passed to "KL_dist_CBpred_func()". Default to 10^(-10).
#' 
#' @import mvtnorm
#' 
#' @return a list bla bla
#' 
#' @author  Lars Erik Gangsei & Hilde Vinje
#' @references Vinje H, Brustad HK, Heggli A, Sevillano CA, Van Son M, 
#' Gangsei LE. Classification of breed combinations for slaughter pigs based 
#' on genotypes-modeling DNA samples of crossbreeds as fuzzy sets from purebred 
#' founders. Front Genet. 2023 Dec 4;14:1289130. doi: 10.3389/fgene.2023.1289130. 
#' PMID: 38116292; PMCID: PMC10729766.
#' 
#' @examples
#' To come


#' @export

classifyCB <- function(DNA=NULL,Predictions = NULL, Train = 'TrainP+',
                           alpha0 = 73.58105,Unkn_LogLike=NULL,PriorDist = NULL,
                          pred_min = 10^(-10))
{
  if(!is.list(DNA))
  {
  # DNA to 0-1-2 coding
  if(class(DNA)[1]=='character'){DNA <- read_ped_data(DNA)}
  if(dim(DNA)[2]==(2*dim(PigBreedPrediction_map)[1])){DNA <- ped2num(DNA)}
    message('Data loaded and prediction is ongoing')
    Scores <- predict(Mod_pls[[Train]],type='scores',newdata = DNA,
                      comps = 1:Mod_pls[[Train]]$ncomp)
    RF_pred <- predict(Classifier_RF[[Train]],newdata = DNA, type = "prob")
    PLSR_pred <- predict(Mod_pls[[Train]], newdata = DNA,
                           type = "response",comps = 1:Mod_pls[[Train]]$ncomp)
    ID <- rownames(DNA)
  }else{
    Scores <-  DNA$Scores
    RF_pred <- DNA$RF_pred
    PLSR_pred <- DNA$PLSR_pred
    ID <- rownames(Scores)
  }
    
  # Define dimentions, breeds to be evaluated and number of PLS components
  nn_new <- dim(Scores)[1]
  mm <- dim(Scores)[2]
  
  # Likelihoods for PLS-QDA
  T_likelihood <- LikelihoodParameters_func(Train,alpha0)
  
  nn_comb <- length(T_likelihood)
  
  # Log - Likelihood for uniform distribution for unknown breed
  if(is.null(Unkn_LogLike))
  {Unkn_LogLike <- -sum(log(apply(apply(Mod_pls[[Train]]$scores,2,range),2,diff)))}
  
  
  # Construct the posterior distribution
  Posterior_distribution <- vector('list',length=nn_comb+1)
  names(Posterior_distribution) <- c(names(T_likelihood),'Unknown')
  
  # Get the prior distribution
  if(is.null(PriorDist)){
    Prior_dist <- split(rep(1,nn_comb+1)/(nn_comb+1),
                        f=names(Posterior_distribution))
  }else{
    if(PriorDist == 'informative')
    {if(Train == 'TrainP+'){Prior_dist <- split(Prior$`P+`$value,f=Prior$`P+`$BreedComb)
    }else{Prior_dist <- split(Prior$`P-`$value,f=Prior$`P-`$BreedComb)}
    }else{Prior_dist <- split(c(rep(1-PriorDist,nn_comb)/nn_comb,PriorDist),
                          f=names(Posterior_distribution))
      }
  }
  
  for(bb in names(T_likelihood))
  {Posterior_distribution[[bb]] <- (dmvnorm(Scores,mean = T_likelihood[[bb]]$mu,
                                                 sigma = T_likelihood[[bb]]$cov,
                                           log=TRUE)+log(Prior_dist[[bb]]))}
  
  Posterior_distribution$Unknown <- rep(Unkn_LogLike+log(Prior_dist$Unknown),nn_new)
  
  Log_posteriors <- matrix(unlist(Posterior_distribution),
                           nn_new,nn_comb+1,byrow=FALSE)
  
  Log_posteriors <- Log_posteriors - matrix(apply(Log_posteriors,1,max),
                                            nn_new,nn_comb+1,byrow=FALSE)
  
  Posteriors <- exp(Log_posteriors)
  Posteriors <- Posteriors/matrix(rowSums(Posteriors),
                                  nn_new,nn_comb+1,byrow=FALSE)
  colnames(Posteriors) <- names(Posterior_distribution)
  
 
  
  res <- data.frame(ID = ID,
                  RF = I(RF_pred),
                  PLSR = I(PLSR_pred),
                  PLSQDA = I(PosteriorToSoft_func(Posteriors,
                             PBs = intersect(c('D','H','L','P','W'),
                                  unique(strsplit(paste(colnames(Posteriors),
                                                        collapse=''),'')[[1]])))),
                  PLSQDA_Posteriors = I(Posteriors),
                  PLSScores = I(Scores))
  
 
  return(res)
  
}
