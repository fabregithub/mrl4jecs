#' MRL LCMRL estimatation scripts
#'

#' MRL.Summary
#'
#' @param fh.labdata lab data
#' @param rnnr rnnr
#' @param biter biter
#' @param Seed seed
#' @param abbr.out abbreviation
#'
#' @export
#'

MRL.Summary <- function(fh.labdata, rnnr=1, biter=200, Seed=0261948, abbr.out=TRUE)
   {
   # ----------------Outputs------------------------------------
   #  via MRL.Stats
   #   =  a character matrix or NULL object.
   #            If a matrix:
   #              Column 1: Laboratory Name
   #              Column 2: LCMRL Calculation
   #                          numerical value if it exists else NA
   #              Column 3: LCMRL Result Flag 1, -1 or NA
   #                          1 = Valid LCMRL
   #                         -1 = LCMRL < lowest Non-Zero Spking Level
   #                         NA = Non-Valid LCMRL or LCMRL cannot be computed
   #            Length of BB.tmp is:
   #               number of labs with Valid LCMRL for original data * boot.iter
   #            Returns NULL if no Labs have a valid LCMRL for orginal data or
   #            inputs are in error.
   #  5 or more files are also output if at least 1 Lab has a valid LCMRL for original
   #          data:
   #
   #  "fh".Rdata.filehandles = filehandles of R adata files with saved objects
   #              and pdf graph objects
   #  MRL.Tables."fh" = .csv file with Tables summarizing the Bayesian Bootstrap
   #             with filename: "Analyte Name".MRL.Tables.csv
   #  "analyte.name".RData = R adata file by simulated MRL data
   #  "fh".LCMRL.Graphs.pdf = pdf files of LCRML graphs
   #  "fh".DensityPlots.pdf = pdf files for density plots of simulated LCMRL values

   # ----------------Inputs-------------------------------------
   #
   # fh.labdata = filehandle (name) with multi laboratory data in csv format
   #        header :
   #           COLUMN 1: Analyte Name
   #           COLUMN 2: Laboratory Name
   #           COLUMN 3: Spiking Level
   #           Column 4: Laboratory Measure at Spiking Level
   #           Column 5: (Optional) Measurement Units. Default = "ug/L"
   # boot.iter = number of desired boostrap iterations. Default boot.iter= 200
   # rnnr = flag for Required Non-Negative Response
   #        1 is for Required Non-Negative Response
   #        0 is for Negative Response Allowed
   # Seed = Random seed starting value for Bayesian Weisgts for consistency of results
   #        Default value set to "0261948". This allows results to be repeated.
   # abbr.out = logical flag for type of output data for MRL summary
   #            TRUE is default
   #            False for extended Output which includeds unweighted table summaries



   # Estimate MRL Summary Stats & Output Files
   MRL.Stats(fh.labdata,rnnr,biter,Seed,Abbr=abbr.out)
   # Set up density plots of Boostrapped LCMRL Estiamtes by Labs, All Labs & Predicted Lab
   infile <-paste(strtrim(fh.labdata,nchar(fh.labdata)-4),".Rdata.filehandles.csv")
   outfile <-paste(strtrim(fh.labdata,nchar(fh.labdata)-4),".DensityPlots.pdf",sep="")
   fh <- read.csv(infile,as.is=TRUE)[,1]
   numfh <- length(fh)
   holder <- NULL
   pdf(outfile,paper="letter")
   for( i in 1:numfh)
      {
      load(as.character(fh[i]))
      show(tmp.den)
      }
   dev.off()
   }


#' LCRML.Values function
#'
#' @param fh file name
#' @param LQL lower quantitation limit, default = 0.5
#' @param UQL upper quantitation limit, default = 1.5
#' @param CPR CPR
#' @param alph alpha error
#' @param bet beta error
#' @param rnnr set 1 if data consisting of only positive values, and 0 if data includes negative values
#'
#' @export
#'
LCMRL.Values <- function(fh, LQL=0.5, UQL=1.5, CPR=0.99, alph=0.05, bet=0.05, rnnr=1)
   {

   # try to read multi-lab data file in csv format

   # values for test runs

   fh.a <- strtrim(fh,nchar(fh)-4)
   outfh <- paste(fh.a,"LCMRL.values","csv",sep=".")

   dat <- try(read.csv(fh,as.is=T),silent=TRUE)
   if( class(dat)[1]=="try-error")
      {
      print("Error:")
      print(paste(".....Problem Reading File:",fh))
      print(geterrmessage())
      return()
      }
   if( length(dat[1,]) < 6)
      {
      print("Error:")
      print(".....Number of input fields is less than 6")
      print(".....Check input file for missing field")
      return()
      }
   if( length(dat[1,]) > 6)
      {
      print("Error:")
      print(".....Number of input fields is greater than 6")
      print(".....Check input file for extra fields")
      return()
      }
   AN <- dat[,1]
   uniq.AN <- unique(AN)
   num.AN <- length(uniq.AN)
   S <- dat[,3]
   M <- dat[,4]
   if ( is.numeric(S) && is.numeric(M))
            {
            } else
            {
            print( "Both the Spiking Levels and Measurents must be numeric")
            return()
            }
   #pdf(paste(fh,".LCMRL.Graphs.pdf",sep=""),height=11,width=8.5)
   A.lab <- vector("character",length(dat[,1]))
   for( i in 1:length(dat[,1])) { A.lab[i] <- paste(dat[i,1],dat[i,2],sep="--")}
   uniq.Alab <- unique(A.lab)
   numA <- length(uniq.Alab)
   Anal <- vector("character",numA)
   Mess <- DLMess <- Anal
   LCMRL <- vector("numeric",numA)
   RF <- DL <- Lc <- Ratio <- LCMRL

   for( kk in 1:numA)
      {
      subtt <- dat[A.lab==uniq.Alab[kk],]

      print("******************************")
      print(paste("**","For:",uniq.Alab[kk]) )
      print("******************************")
      print("")

         S <- subtt[,3]
         M <- subtt[,4]
         AN.t <- subtt[1,1]
         U <- try(subtt[1,6], silent=TRUE)
         if( inherits(U, "try-error") ) { U <-"ug/L" }

         tmp <-RobLCMRL(S,M,AN.t,U,lowQL=LQL,upQL=UQL,CovProbReq=CPR,alpha=alph,
                  beta=bet,ReqNonNegResponse=rnnr,grphs=0,labnum=NULL)
         Anal[kk] <- uniq.Alab[kk]
         LCMRL[kk] <- tmp[[1]]
         DL[kk] <- tmp[[2]]
         Lc[kk] <- tmp[[3]]
         Ratio[kk] <- tmp[[3]]/tmp[[2]]
         RF[kk] <- tmp[[4]]
         Mess[kk] <- tmp[[5]]
         DLMess[kk] <- tmp[[7]]
      }
   Table <- data.frame(Anal,LCMRL,DL,Lc,Ratio,RF,Mess,DLMess)
   dimnames(Table)[[2]] <- c("Analyte","LCMRL","DL","Lc","LCMRL/DL","ResultFlag","Message","DLMessage")
   write.csv(Table,outfh,row.names=FALSE)
   return("end of processing")
   }

#' LCMRL.Graphs
#'
#' @param fh file name
#' @param LQL lower quantitation limit, default = 0.5
#' @param UQL upper quantitation limit, default = 1.5
#' @param CPR CPR
#' @param alph alpha error
#' @param bet beta error
#' @param rnnr set 1 if data consisting of only positive values, and 0 if data includes negative values
#'
#' @export
#'
LCMRL.Graphs <- function(fh, LQL=0.5, UQL=1.5, CPR=0.99, alph=0.05, bet=0.05, rnnr=1)
   {
   fh.a <- strtrim(fh,nchar(fh)-4)
   outfh <- paste("LCMRL.Values.Tables.",fh,sep="")


   dat <- try(read.csv(fh,as.is=T),silent=TRUE)
   if( class(dat)[1]=="try-error")
      {
      print("Error:")
      print(paste(".....Problem Reading File:",fh))
      print(geterrmessage())
      return()
      }
   if( length(dat[1,]) < 6)
      {
      print("Error:")
      print(".....Number of input fields is less than 6")
      print(".....Check input file for missing field")
      return()
      }
   if( length(dat[1,]) > 6)
      {
      print("Error:")
      print(".....Number of input fields is greater than 6")
      print(".....Check input file for extra fields")
      return()
      }
   AN <- dat[,1]
   uniq.AN <- unique(AN)
   num.AN <- length(uniq.AN)
   S <- dat[,3]
   M <- dat[,4]
   if ( is.numeric(S) && is.numeric(M))
            {
            } else
            {
            print( "Both the Spiking Levels and Measurents must be numeric")
            return()
            }
   pdf(paste(fh.a,".LCMRL.graphs.pdf",sep=""),height=11,width=8.5)
   for( kk in 1:num.AN)
      {
      tt <- dat[dat[,1]==uniq.AN[kk],]
      LabNum <- as.character(unique(tt[,2]))
      numlabs <- length(LabNum)


      print("******************************")
      print(paste("**","For:",uniq.AN[kk]) )
      print("******************************")
      print("")

      #labflag <- vector("numeric",numlabs)
      time1 <- Sys.time()

      print(paste("-*-*-*-* Detection Values for Labs"))
      for( i in 1:numlabs)
         {
         subtt <-tt[as.character(tt[,2])==LabNum[i],]
         S <- subtt[,3]
         M <- subtt[,4]
         AN.t <- subtt[1,1]
         U <- try(subtt[1,6], silent=TRUE)
         if( inherits(U, "try-error") ) { U <-"ug/L" }
         RobLCMRL(S,M,AN.t,U,lowQL=LQL,upQL=UQL,CovProbReq=CPR,alpha=alph,
                  beta=bet,ReqNonNegResponse=rnnr,grphs=1,labnum=LabNum[i])
         }
      }
      dev.off()
      time2 <- Sys.time()
      print(paste("--------Time to do Detection Values for All Labs", difftime(time2,time1,units="hours")))
      # remove cases where robLCMRL.dev did not have and LCMRL

   return("end of processing")
   }

#' MRL.Stats
#'
#' @param fh file
#' @param rnnr rnnr
#' @param boot.iter iteration number
#' @param Seed seed
#' @param Abbr abbreviation
#'
#' @export
#'
MRL.Stats <- function(fh, rnnr=1, boot.iter=200, Seed=0261948, Abbr=TRUE)
   {

   # ----------------Outputs------------------------------------
   #  BB.tmp2 =  a character matrix or NULL object.
   #            If a matrix:
   #              Column 1: Laboratory Name
   #              Column 2: LCMRL Calculation
   #                          numerical value if it exists else NA
   #              Column 3: LCMRL Result Flag 1, -1 or NA
   #                          1 = Valid LCMRL
   #                         -1 = LCMRL < lowest Non-Zero Spking Level
   #                         NA = Non-Valid LCMRL or LCMRL cannot be computed
   #            Length of BB.tmp is:
   #               number of labs with Valid LCMRL for original data * boot.iter
   #            Returns NULL if no Labs have a valid LCMRL for orginal data or
   #            inputs are in error.
   #  3 files are also output if at least 1 Lab has a valid LCMRL for original
   #          data:
   #
   #  "fh".Rdata.filehandles = filehandles of R adata files with saved objects
   #              and pdf graph objects
   #  MRL.Tables."fh" = .csv file with Tables summarizing the Bayesian Bootstrap
   #             with filename: "Analyte Name".MRL.Tables.csv
   #  .RData = R adata file by simulated MRL data
   # ----------------Inputs-------------------------------------
   #
   # fh = filehandle (name) with multi laboratory data in csv format
   #        header :
   #           COLUMN 1: Analyte Name
   #           COLUMN 2: Laboratory Name
   #           COLUMN 3: Spiking Level
   #           Column 4: Laboratory Measure at Spiking Level
   #           Column 5: (Optional) Measurement Units. Default = "ug/L"
   # boot.iter = number of desired boostrap iterations. Default boot.iter= 200
   # Seed = Random seed starting value for Bayesian Weisgts for consistency of results
   #        Default value set to "0261948". This allows results to be repeated.
   # rnnr = flag for Required Non-Negative Response
   #        1 is for Required Non-Negative Response
   #        0 is for Negative Response Allowed
# -------------------------Begin MRL.Stats ------------------------------

   # try to read multi-lab data file in csv format

   # values for test runs

   fh.a <- strtrim(fh,nchar(fh)-4)
   outfh <- paste("MRL.Tables.",fh,sep="")

   # open output file for MRL Tables
   write(paste("MRL.Tables.",fh.a,sep=""),outfh)
   write("",outfh,append=TRUE)
   write("",outfh,append=TRUE)
   write(paste("Random Seed =",Seed),outfh,append=TRUE)
   write("",outfh,append=TRUE)

   dat <- try(read.csv(fh,as.is=T),silent=TRUE)
   if( class(dat)[1]=="try-error")
      {
      print("Error:")
      print(paste(".....Problem Reading File:",fh))
      print(geterrmessage())
      return()
      }
   if( length(dat[1,]) < 6)
      {
      print("Error:")
      print(".....Number of input fields is less than 6")
      print(".....Check input file for missing field")
      return()
      }
   if( length(dat[1,]) > 6)
      {
      print("Error:")
      print(".....Number of input fields is greater than 6")
      print(".....Check input file for extra fields")
      return()
      }
   AN <- dat[,1]
   uniq.AN <- unique(AN)
   num.AN <- length(uniq.AN)
   S <- dat[,3]
   M <- dat[,4]
   if ( is.numeric(S) && is.numeric(M))
            {
            } else
            {
            print("Error")
            print( "....Both the Spiking Levels and Measurents must be numeric")
            return()
            }
   goodflag <- 0
   for( kk in 1:num.AN)
      {
      tt <- dat[dat[,1]==uniq.AN[kk],]
      LabNum <- as.character(unique(tt[,2]))
      numlabs <- length(LabNum)
      BB.tmp <- matrix(0,boot.iter*numlabs,3)
      goodlab <- 0
      goodlab.id <- NULL
      goodlab.message <- NULL
      badlab <- 0
      badlab.id <- NULL
      badlab.message <- NULL

      print("******************************")
      print(paste("**","For:",uniq.AN[kk]) )
      print("******************************")
      print("")

      write("",outfh,append=TRUE)
      write("",outfh,append=TRUE)
      write("************************************************************************",outfh,append=TRUE)
      write(paste(uniq.AN[kk],"   MRL: Summary Tables"),outfh,append=TRUE)
      write("************************************************************************",outfh,append=TRUE)
      write("",outfh,append=TRUE)

      #labflag <- vector("numeric",numlabs)
      time1 <- Sys.time()

      print(paste("-*-*-*-* Determining which labs have usable LCMRLs"))
      for( i in 1:numlabs)
         {
         subtt <-tt[as.character(tt[,2])==LabNum[i],]
         S <- subtt[,3]
         M <- subtt[,4]
         AN.t <- subtt[1,1]
         U <- try(subtt[1,6], silent=TRUE)
         if( inherits(U, "try-error") ) { U <-"ug/L" }

         tmp <-RobLCMRL(S,M,AN.t,ReqNonNegResponse=rnnr,grphs=0,labnum=LabNum[i])

         if(tmp$LCMRLResultFlag == 1 || tmp$LCMRLResultFlag == -1)
            {
            print(paste("--------Processing Lab", LabNum[i]))

            goodlab <- goodlab + 1
            goodlab.id <- c(goodlab.id,i)
            goodlab.message <- c(goodlab.message,tmp$LCMRLResultMessage)
            BB.tmp[((goodlab-1)*boot.iter+1):(goodlab*boot.iter),]  <-
               BB.MRL(S,M,AN.t,U,biter=boot.iter,labnum=LabNum[i],Seed,rnnr)
             gc(verbose=FALSE)
            if(goodlab == 1)
               {
               LCMRL.val <- tmp[[1]]
               labflag <- tmp$LCMRLResultFlag
               } else
               {
               LCMRL.val <- c(LCMRL.val,tmp[[1]])
               labflag <- c(labflag,tmp$LCMRLResultFlag)
               }
            } else
            {
            badlab <- badlab + 1
            badlab.id <- c(badlab.id,LabNum[i])
            badlab.message <- c(badlab.message,tmp$LCMRLResultMessage)
            }
         }

      time2 <- Sys.time()
      print(paste("--------Time to do 200 iterations per Lab", difftime(time2,time1,units="hours")))
      # remove cases where robLCMRL.dev did not have and LCMRL
      BB.tmp2 <- BB.tmp[BB.tmp[,1] != "0",]

      # reformat as data.frame
      BB.tmp2 <- data.frame(BB.tmp2[,1],as.numeric(BB.tmp2[,2]),as.numeric(BB.tmp2[,3]))

      # open output file for filenames of .Rdata files with pdf graphics objects
      if( goodlab >= 2)
         {
         goodflag <- goodflag + 1
         if(goodflag == 1)
            {
            write("FileName",paste(fh.a,".Rdata.filehandles.csv"),sep="")
            }
         }

      #  Depending on number of goodlabs within Analyte determine output info
      if( goodlab == 0)
         {
         write("No Laboratory had a useable LCMRL for the original Data",outfh,append=TRUE)
         write("",outfh,append=TRUE)
         if(badlab > 0)
            {
            for( jj in 1:badlab)
               {
               write(paste("......Lab",badlab.id[jj], badlab.message),outfh,append=TRUE)
               }
            write("",outfh,append=TRUE)
            }
         }
      if( goodlab == 1 )
         {
         write("Only one laboratory had a good LCMRL for the original Data",outfh,append=TRUE)
         write("Therefore no MRL calculation was performed",outfh,append=TRUE)
         write("",outfh,append=TRUE)
         if(badlab > 0)
            {
            for( jj in 1:badlab)
               {
               write(paste("......Lab",badlab.id[jj], badlab.message[jj]),outfh,append=TRUE)
               }
            write(paste("......Lab",LabNum[goodlab.id[1]],"had LCMRL value ",LCMRL.val[1],U,"",goodlab.message),outfh,append=TRUE)
            write("",outfh,append=TRUE)
            }
         }
      if( goodlab == 2)
         {
         write("Only two labortries had a useable LCMRL for the original Data",outfh,append=TRUE)
         write("Therefore only the results for ALL LABS should be used as a preliminary MRL Estimate",outfh,append=TRUE)
         write("",outfh,append=TRUE)
         if(badlab > 0)
            {
            for( jj in 1:badlab)
               {
               write(paste("......Lab",badlab.id[jj], badlab.message[jj]),outfh,append=TRUE)
               }
            write("",outfh,append=TRUE)
            }
         tmp.den <- MRL.Output(uniq.AN[kk],U,BB.tmp2,labflag,LCMRL.val,boot.iter,outfh,fh.a,Abbr)
         Anal <- uniq.AN[kk]
         save(Anal,U,BB.tmp2,labflag,LCMRL.val,boot.iter,outfh,fh.a,tmp.den,
               goodlab,badlab,badlab.id,badlab.message,
               file=paste(uniq.AN[kk],".",fh.a,".RData",sep=""))
         write(paste(uniq.AN[kk],".",fh.a,".RData",sep=""),paste(fh.a,".Rdata.filehandles.csv"),sep="",append=T)
         }
      if( goodlab >=3 )
         {
         write("Valid MRL calculated",outfh,append=TRUE)
         write("Use results for Predicted Lab as an MRL Estimate",outfh,append=TRUE)
         write("",outfh,append=TRUE)
         if(badlab > 0)
            {
            for( jj in 1:badlab)
               {
               write(paste("......Lab",badlab.id[jj], badlab.message[jj]),outfh,append=TRUE)
               }
            write("",outfh,append=TRUE)
            }
         tmp.den <- MRL.Output(uniq.AN[kk],U,BB.tmp2,labflag,LCMRL.val,boot.iter,outfh,fh.a,Abbr)
         Anal <- uniq.AN[kk]
         save(Anal,U,BB.tmp2,labflag,LCMRL.val,boot.iter,outfh,fh.a,tmp.den,
               goodlab,badlab,badlab.id,badlab.message,
               file=paste(uniq.AN[kk],".",fh.a,".RData",sep=""))
         write(paste(uniq.AN[kk],".",fh.a,".RData",sep=""),paste(fh.a,".Rdata.filehandles.csv"),sep="",append=T)
         }
      }
   return("end of processing")
   }

#' MRL.Output
#'
#' @param An An
#' @param U U
#' @param BB.tmp2 BB.tmp2
#' @param labflag labflag
#' @param LCMRL.val LCMRL value
#' @param boot.iter number of iterations
#' @param outfh output file
#' @param fh.a file
#' @param Abbr abbreviations
#'
MRL.Output <- function(An, U, BB.tmp2, labflag, LCMRL.val, boot.iter, outfh, fh.a, Abbr=TRUE)
   {

   #densityplots of individual Labs
   labs <- unique(BB.tmp2[,1])
   numlabs <- length(labs)
   minLCMRL <- 0.0
   maxLCMRL <- max(as.numeric(BB.tmp2[,2]),na.rm=TRUE)+0.1

   # Create Summary Output Table

   # Create Table 1 Summary of Bootstrap Iteration Counts
   Table1 <- matrix(0,8,numlabs+1)

   Table1[2,1:numlabs] <- rep(boot.iter,numlabs)
   Table1[2,numlabs+1] <- sum(as.numeric(Table1[2,1:numlabs]))
   for(i in 1:numlabs)
      {
      Table1[3,i] <- sum(as.numeric((BB.tmp2[,1]==labs[i]))*(BB.tmp2[,3]==1),na.rm=TRUE)
      Table1[4,i] <- sum(as.numeric(BB.tmp2[,1]==labs[i])*(BB.tmp2[,3]==-1),na.rm=TRUE)
      Table1[5,i] <- sum(as.numeric(BB.tmp2[,1]==labs[i])*(is.na(BB.tmp2[,3])))
      }
   for(i in 1:numlabs)
      {
      Table1[6,i] <- 100*Table1[3,i]/Table1[2,i]
      Table1[7,i] <- 100*Table1[4,i]/Table1[2,i]
      Table1[8,i] <- 100*Table1[5,i]/Table1[2,i]
      }
   for(i in 3:5)
      {
      Table1[i,numlabs+1] <- sum(as.numeric(Table1[i,1:numlabs]))
      }
   tot <- sum(as.numeric(Table1[2,1:numlabs]))
   Table1[6,numlabs+1] <- 100* sum(as.numeric(Table1[3,1:numlabs]))/tot
   Table1[7,numlabs+1] <- 100* sum(as.numeric(Table1[4,1:numlabs]))/tot
   Table1[8,numlabs+1] <- 100* sum(as.numeric(Table1[5,1:numlabs]))/tot

   dimnames(Table1)[[2]] <- c(as.character(labs),"All Labs")
   Table1[1,] <- rep("Valid",numlabs+1)
   Table1[1,labflag==-1] <- "<min(SL)"
   Table1[1,numlabs+1] <- ""

   write("",outfh,append=TRUE)
   write(paste("Measurement Units = ",U),outfh,append=TRUE)
   write("",outfh,append=TRUE)
   write(paste("Table 1:",An,"Iteration Count Table"),outfh,append=TRUE)
   write("",outfh,append=TRUE)
   write.csv(Table1,outfh,
         row.names=c("LCMRL Value","Num Iter","Valid LCMRL","< min(SL)",
         "NA","% Valid","% < min(SL)","% NA"), append=TRUE)
   write("",outfh,append=TRUE)

   # Create Table 2 Summary Stats for Individual Unweighted Labs, All Labs & Predicted Labs
   if(!Abbr)
   {
   BB.x <- na.omit(data.frame(BB.tmp2[,1],as.numeric(BB.tmp2[,2]),
                   stringsAsFactors=TRUE))
   oldBBx <- BB.x
   # get lambda from box cox multivariate normal
   lamb <- powerTransform(BB.x[,2] ~ BB.x[,1])$lambda
   if( lamb > 0 )
      { BB.x[,2] <- BB.x[,2]^lamb} else
      { BB.x[,2] <- log(BB.x[,2])}
   BB.labx <- split(BB.x[,2],BB.x[,1])
   BB.locvar <- sapply(BB.labx,roblocvar4,cc=6)
   #within lab variance  is the mean of the robust lab variances
   with.var <-  mean(as.numeric(BB.locvar[2,]))
   # between lab variance is variance of the robust lab means
   betw.var <- var(as.numeric(BB.locvar[1,]))
   # grand mean is the mean of the robust lab means
   x.bar <- mean(as.numeric(BB.locvar[1,]))
   k <- length(BB.labx)
   S.L <- (1+1/k)*betw.var + with.var
   xij <- sort(x.bar + sqrt(S.L)*unlist(std.z(BB.labx,BB.locvar)))
   xij <- xij[xij > 0]
   if( lamb > 0 )
      {
      xij <- xij[xij>0]
      xij <- exp(log(xij)/lamb)
      } else
      {
      xij <- exp(xij)
      }

   BB.x <- oldBBx
   Table2 <- matrix(0,22,numlabs+2)
   tlabs <-c(as.character(labs),"ALL Labs","Predicted Lab")
   dimnames(Table2)[[2]] <- tlabs
   pcts <-c(0.05,0.10,0.20,0.25,0.50,0.60,0.70,0.75,0.80,0.85,0.90,0.95,0.96,0.97,0.98,0.99)

   for(i in 1:numlabs)
      {
      Table2[2,i] <-min(BB.x[BB.x[,1]==tlabs[i],2],na.rm=T)
      Table2[19,i] <- max(BB.x[BB.x[,1]==tlabs[i],2])
      }
   for(i in 1:numlabs)
      {
      st <- sort(BB.x[BB.x[,1]==tlabs[i],2])
      Table2[3:18,i] <- quantile(st,probs=pcts,names=F)
      Table2[20,i] <- mean(st)
      Table2[21,i] <- sd(st)
      Table2[22,i] <- Table2[21,i]/Table2[20,i]
      }
   Table2[2,numlabs+1] <-min(BB.x[,2])
   Table2[19,numlabs+1] <- max(BB.x[,2])
   st <- sort(BB.x[,2])
   Table2[3:18,numlabs+1] <- quantile(st,probs=pcts,names=F)
   Table2[20,numlabs+1] <- mean(st)
   Table2[21,numlabs+1] <- sd(st)
   Table2[22,numlabs+1] <- Table2[21,numlabs+1]/Table2[20,numlabs+1]

   Table2[2,numlabs+2] <-min(xij)
   Table2[19,numlabs+2] <- max(xij)
   Table2[3:18,numlabs+2] <- quantile(xij,probs=pcts,names=F)
   Table2[20,numlabs+2] <- mean(xij)
   Table2[21,numlabs+2] <- sd(xij)
   Table2[22,numlabs+2] <- Table2[21,numlabs+2]/Table2[20,numlabs+2]
   Table2[1,] <- c(LCMRL.val,"","")
   write("Table2: Unweighted Data: Percentiles of LCMRL Bootstrap Dist.",
      outfh,append=TRUE)
   write("",outfh,append=TRUE)
   write.csv(Table2,outfh,
      row.names=c("LCMRL Value","minimum","5%","10%","20%","25%","50%","60%","70%","75%",
      "80%","85%","90%","95%","96%","97%","98%","99%","maximum",
      "mean","std.dev","cv"),append=TRUE)
   write("",outfh,append=TRUE)




   # Create Table 3 Guttman UTLs for Unweighted Data

   write("Table3: Unweighted Data: Upper Tol. Limits for 95% Probablity of p% Coverage",
            outfh,append=TRUE)
   write("",outfh,append=TRUE)

   Table3 <- matrix(0,6,numlabs+2)
   dimnames(Table3)[[2]] <-c(as.character(labs),"ALL Labs","Predicted Lab")
   #spaces <- c(""," ","   ","    ","     ","      ","       ","       ")
   #names(Table3) <- list(spaces[1:numlabs+2])
   #for(i in 1:(numlabs)) { Table3[i[i]] <- c("","","","","") }
   Table3[1,1:numlabs] <- c(LCMRL.val)
   Table3[2:6,numlabs+1] <- c(guttman.utl(st,0.95,0.70),
                  guttman.utl(st,0.95,0.75),guttman.utl(st,0.95,0.80),
                  guttman.utl(st,0.95,0.90),guttman.utl(st,0.95,0.95))
   Table3[2:6,numlabs+2] <- c(guttman.utl(xij,0.95,0.70),
                  guttman.utl(xij,0.95,0.75),guttman.utl(xij,0.95,0.80),
                  guttman.utl(xij,0.95,0.90),guttman.utl(xij,0.95,0.95))
   Table3[1,(numlabs+1):(numlabs+2)] <- c("","")
   Table3[2:6,1:numlabs] <- ""
   write.csv(Table3,outfh,append=TRUE,row.names=c("LCMRL Value","95-70 UTL","95-75 UTL",
                                                "95-80 UTL","95-90 UTL",
                                                "95-95 UTL"))
   write("",outfh,append=TRUE)
   }
   # Create Table 4  Guttman UTLs for One sided Tukey Weighted Data

   if(!Abbr)
   {
   write("Table4: Weighted Data: Upper Tol. Limits for 95% Probablity of p% Coverage",
            outfh,append=TRUE)
   write("",outfh,append=TRUE)
   } else
   {
   write("Table2: Weighted Data: Upper Tol. Limits for 95% Probablity of p% Coverage",
            outfh,append=TRUE)
   write("",outfh,append=TRUE)
   }

   BB.x <- na.omit(data.frame(BB.tmp2[,1],as.numeric(BB.tmp2[,2]),
                   stringsAsFactors=TRUE))
   Table4 <- matrix(0,6,numlabs+2)
   dimnames(Table4)[[2]] <-c(as.character(labs),"ALL Labs","Predicted Lab")
   st <- sort(BB.x[,2])
   median.st <- median(st)
   mad.st <- mad(st)
   st <- sort(st)
   tw1.lst <- wt.bw.1.sided(st,median.st,mad.st,cc=6,reltol=1e-6,maxsteps=100)
   #csum.st <- cumsum(tw1.lst[[2]])
   new.loc <- tw1.lst[[1]]
   utl70 <- wgt.guttman.utl(st,tw1.lst[[2]],.95,.70)
   utl75 <- wgt.guttman.utl(st,tw1.lst[[2]],.95,.75)
   utl80 <- wgt.guttman.utl(st,tw1.lst[[2]],.95,.80)
   utl90 <- wgt.guttman.utl(st,tw1.lst[[2]],.95,.90)
   utl95 <- wgt.guttman.utl(st,tw1.lst[[2]],.95,.95)
   Table4[1,1:numlabs] <- c(LCMRL.val)
   Table4[2:6,numlabs+1] <- c(utl70,utl75,utl80,utl90,utl95)

   # get lambda from box cox multivariate normal
   lamb <- powerTransform(BB.x[,2] ~ BB.x[,1])$lambda
   if( lamb > 0 )
      { BB.x[,2] <- BB.x[,2]^lamb} else
      { BB.x[,2] <- log(BB.x[,2])}
   med.x <- median(BB.x[,2])
   mad.x <- mad(BB.x[,2])
   BB.x.wgt.lst <- wt.bw.1.sided(BB.x[,2],med.x,mad.x,cc=6,reltol=1e-6,maxsteps=100)
   if(any(BB.x.wgt.lst[[2]] == 0 ))
      {BB.x.wgt.lst <- wt.bw.1.sided(BB.x[,2],med.x,mad.x,cc=9,reltol=1e-6,maxsteps=100)}
   BB.x <- cbind(BB.x[1],BB.x[2],BB.x.wgt.lst[[2]])
   BB.labx <- split(BB.x[,2],BB.x[,1])
   BB.labwts <- split(BB.x[,3],BB.x[,1])
   BB.labels <- split(BB.x[,1],BB.x[,1] )
   BB.locvar <- matrix(0,2,length(BB.labx))
   nlabs <-length(BB.labx)
   for(i in 1:nlabs)
      {
      tmp <- roblocvar4(BB.labx[[i]],cc=6)
      BB.locvar[1,i] <- tmp[[1]]
      BB.locvar[2,i] <- tmp[[2]]
      }
   sum.wgts <-sapply(BB.labwts,sum)
   #within lab variance is weighted average of lab robust means
   with.var <-  sum(as.numeric(BB.locvar[2,])*sum.wgts )
   # between lab variance is the weighted variance of the robust lab means
   betw.var <- w.var(as.numeric(BB.locvar[1,]),sum.wgts)
   # grand mean
   x.bar <- sum(as.numeric(BB.locvar[1,])*sum.wgts)
   k <- length(BB.labx)
   S.L <- (1+1/k)*betw.var + with.var
   BB.labx2 <- BB.labx
   for(i in 1:nlabs)
      {
      BB.labx2[[i]]<- x.bar + sqrt(S.L)*std2.z(BB.labx[[i]],BB.locvar[,i])
      }
   sxij <- sort(unlist(BB.labx2),index=TRUE)
   # check for useable values
   if(lamb > 0 )
      {
      xij <- sxij$x
      for( i in 1:length(sxij$x))
         {
         xij[i] <- ifelse(sxij$x[i]>0,exp(log(sxij$x[i])/lamb),0)
         }
      } else
      {
      xij <- exp(sxij$x)
      }
   wt.x <- unlist(BB.labwts)[sxij$ix]
   BB.labels <- unlist(BB.labels)[sxij$ix]



   utl70 <- wgt.guttman.utl(xij,wt.x,.95,.70)
   utl75 <- wgt.guttman.utl(xij,wt.x,.95,.75)
   utl80 <- wgt.guttman.utl(xij,wt.x,.95,.80)
   utl90 <- wgt.guttman.utl(xij,wt.x,.95,.90)
   utl95 <- wgt.guttman.utl(xij,wt.x,.95,.95)
   Table4[2:6,numlabs+2] <- c(utl70,utl75,utl80,utl90,utl95)
   Table4[1,(numlabs+1):(numlabs+2)] <- c("","")
   Table4[2:6,1:numlabs] <- ""
   write.csv(Table4,outfh,append=TRUE,row.names=c("LCMRL Value","95-70 UTL","95-75 UTL",
                                                "95-80 UTL","95-90 UTL",
                                                "95-95 UTL"))

   comb.dat <- cbind("All Labs",st)
   comb.dat <- rbind(comb.dat,cbind("Predicted Lab",xij))

   pos1 <- quantile(density(xij)$y,0.9)
   pos2 <- quantile(density(log10(xij))$y,0.9)
   BB.p <- na.omit(data.frame(BB.tmp2[,1],as.numeric(BB.tmp2[,2]),
                   stringsAsFactors=TRUE))
   BB.lst <- split(BB.p[,2],BB.p[,1])
   BB.pos <- sapply(BB.lst,density)
   BB.q <-sapply(BB.pos[2,],quantile,probs=0.9)
   posA <- max(BB.q)
   BB.log <- sapply(BB.pos[2,],log10)
   BB.q <- sapply(BB.log,quantile,probs=0.9)
   posB <- max(BB.q)
   #pdf(paste(An,".","DenPlot.",".pdf",sep=""),height=11,width=8.5)
      # DENSITYPLOTS of individual Labs
      dplot1 <- densityplot(~BB.p[,2]|as.factor(BB.p[,1]),
         layout=c(1,numlabs),xlab="Spiking Level",
         na.rm=TRUE,pch="|",cex=0.2,
         main =paste("Density Plots of Labs \n",paste(An,"  ",U)),
         panel=function(x)
            {
            panel.densityplot(x)
            panel.abline(v=utl75,col="red")
            panel.text(utl75,posA,labels=paste(" MRL = ",sprintf("%.2f",utl75)),col=2,pos=4,cex=0.8)
            }
         )
       dplot2 <- densityplot(~BB.p[,2]|as.factor(BB.p[,1]),
         layout=c(1,numlabs),xlab="Spiking Level",
         na.rm=TRUE,pch="|", cex=0.2, scales=list(x=list(log = 10)),
         main =paste("Density Plots of Labs \n",paste(An,"  ",U)),
         panel=function(x)
            {
            panel.densityplot(x)
            panel.abline(v=log10(utl75),col="red")
            panel.text(log10(utl75),posB,labels=paste(" MRL = ",sprintf("%.2f",utl75)),col=2,pos=4,cex=0.8)
            }
         )

      # DENSITYPLOTS of Weighted ALL LABS & PREDICTED LAB
      dplot3 <- densityplot(~as.numeric(comb.dat[,2])|as.factor(comb.dat[,1]),
         layout=c(1,2),xlab="Spiking Level", na.rm=TRUE,
         pch="|",cex=0.2, #xlim =c(min(as.numeric(comb.dat[,2])),max(as.numeric(comb.dat[,2]))),
         main =paste("Density Plots of Combined Labs \n",paste(An,"  ",U)),
         panel=function(x)
            {
            panel.densityplot(x)
            panel.abline(v=utl75,col="red")
            panel.text(utl75,pos1,labels=paste(" MRL = ",sprintf("%.2f",utl75)),col=2,pos=4,cex=0.8)
            })
      dplot4 <- densityplot(~as.numeric(comb.dat[,2])|as.factor(comb.dat[,1]),
         layout=c(1,2),xlab="Spiking Level", na.rm=TRUE,scales=list(x=list(log = 10)),
         pch="|",cex=0.2,
         main =paste("Density Plots of Combined Labs \n",paste(An,"  ",U)),
         panel=function(x)
            {
            panel.densityplot(x)
            panel.abline(v=log10(utl75),col="red")
            panel.text(log10(utl75),pos2,labels=paste(" MRL = ",sprintf("%.2f",utl75)),col=2,pos=4,cex=0.8)
            }
            )
   #dev.off()
   write("",outfh,append=TRUE)
   return(list(dplot1,dplot2,dplot3,dplot4))
   }

#' BB.MRL
#'
#' @param S S
#' @param M M
#' @param AN AB
#' @param U U
#' @param biter biter
#' @param labnum lab number
#' @param Seed seed
#' @param rnnr rnnr
#'
#' @export
#'
 BB.MRL <- function(S, M, AN, U, biter=200, labnum=NULL, Seed, rnnr=1)
   {
   #
   # This code performs a bayesian bootstrap on the spiking level and
   # measurement data generating a distribition of the LCMRL values.
   #
   # ----------------Inputs-------------------------------------
   #
   # S - vector of Spiking levels
   # M - vector of Measurents at Spiking Levels
   # AN - analyte name
   # U- measurement units. Default = 'ug/l'
   # biter - number of boostrap iterations. Default = 100
   # Seed - Random Seed. Default set at 0261948. Ensures consectutive runs give
   #        identical results
   # rnnr = flag for Required Non-Negative Response
   #        1 is for Required Non-Negative Response
   #        0 is for Negative Response Allowed
   # ----------------Outputs ------------------
   #
   # BB.mat - a matrix of LCMRL values and Result flags with dim = (biter,2 or 3)
   #     BB.mat[,1] - LCMRL values for each boostrap iteration
   #     BB.mat[,2] - LCMRLResultflag
   # ----------------Start code---------------------------------
   #

   boot.iter <- biter
   ncol <- ifelse(is.null(labnum),2,3)
   BB.1.mat <-matrix(0,boot.iter,ncol)
   dimnames(BB.1.mat)[[2]] <- list('Lab','LCMRL','ResultFlag')

   # get the Bayesian Bootstrap weights
   BB.1.wgts <- Bayes.1.Boot(S,boot.iter,Seed)

   for(i in 1:boot.iter)
      {
      tmptry <-try(RobLCMRL(S,M,AN,U,ReqNonNegResponse=rnnr,grphs=0,BW=BB.1.wgts[i,]))
      gc(verbose=FALSE)
      #print(paste("### boot.iter =",i))
      #print(tmptry)
      #if ( tmptry[2] == -1 || tmptry[2] == 1 )
      if ( class(tmptry)[1]=="try-error")
         {
         BB.1.mat[i,]<- c(labnum,NA,NA)
         }
         else
         {
         ifelse( (tmptry[2] == -1 || tmptry[2] == 1 ),
                  BB.1.mat[i,]<-c(labnum,tmptry),
                  BB.1.mat[i,]<- c(labnum,NA,NA) )
         }
      if(i%%10==0) print(paste('#########',i, "iterations completed"))
      }
   return(BB.1.mat)
   }

#' RobLCMRL
#'
#' @param SpikeConc spiked concentration
#' @param MeasConc measured concentration
#' @param AnalName analyte names
#' @param strUnits Units
#' @param lowQL lower qunatiation limit
#' @param upQL upper quantitation limit
#' @param CovProbReq covprobreq
#' @param alpha alpha
#' @param beta beta
#' @param ReqNonNegResponse required non negative response
#' @param grphs graphs
#' @param labnum lab number
#' @param BW BW
#'
#' @export
#'
 RobLCMRL <- function(SpikeConc, MeasConc, AnalName,
                     strUnits="ug/L", lowQL=0.5, upQL=1.5,
                     CovProbReq=0.99, alpha=0.05, beta=0.05, ReqNonNegResponse=1,
                     grphs=1, labnum=NULL, BW =rep(1,length(SpikeConc)))
# ----------------Outputs------------------------------------
# Returns a list with names:
#
# LCMRLval = numerical value of LCMRL
# DLval = numerical value of Hubaux-Vos detection limit
# CLval = numerical value of critical level
# R = correlation : numerical
# R2 = R-squared : numerical
# R2adj = adjusted R-squared L: numerical
# LCMRLResultFlag = result flag for LCMRL : numerical
# LCMRLResultMessage = result message for LCMRL : string
#     -4 'Aborted: Not enough non-zero Spiking Levels
#     -3 'Nonconvergence: '
#     -2 'LCMRL is above highest spiking level'
#     -1 'Lower spiking level needed to bracket the LCMRL'
#      0 - not set
#      1 'Valid LCMRL'
# DLResultFlag = result flag for HV DL  : numerical
# DLResultMessage = result message for HV DL  : string
#     -3 - 'Nonconvergence: ' output.message
#     -2 - 'PROBLEM: DL appears to be above max spiking level'
#      0 - not set
#      1 - 'Valid DL'
#      2  - 'DL calculated >= LCMRL; set DL = LCMRL'
# Var.Model = a glm object representing the Var Model
# Means.Model = a lm object representing the Means Model
# MSE.Model = a glm object representing the Condtional MSE Model
# ----------------Inputs-------------------------------------
#
# SpikeConc = vector of spiking concentrations
# MeasConc = vector of measured concentrations
# AnalName = analyte name
# strUnits = units of measurement Default = "mg/L"
# lowQL = lower quality limit in LCMRL specification  Default=0.5
# upQL = upper quality limit in LCMRL specification   Default =1.5
# CovProbReq = quality interval coverage probability requirement in LCMRL
#   specification  Default = 0.99
# alpha = probability defining critical response (yc) in HV DL  Default = 0.05
# beta = probability defining DL in HV DL: Pr(y <= yc | x = DL),  Default =0.05
# ReqNonNegResponse = Is instrument data non-negative? (1=Y, 0=N) Default = 1
# grphs = switch for graphical output 0= no output, 1= LCMRL output, Default=1
# labnum = laboratory number. Used only with MRL estimation in multi-lab-studies
# BW - normalized vector of Bayesian Weights used for MRL calculations.
#      Default is set to vector of 1's the length of SpikeConc to produce the
#      original LCMRL calculation. For the MRL the weights are generated by
#      Bayes.1.Boot() in BB.MRL()
# -------------------------Begin RobLCMRL ------------------------------
   {
   theta <- 1e-6
   maxIter <- 100
   DnuSL <- NULL
   LowLimLCMRL <- 0

   xs <- sort(SpikeConc,index=TRUE)
   x <- xs$x
   y <- MeasConc[xs$ix]
   xorg <- x    # save orgibnal x values
   yorg <- y    # save original y values
   BWorig <- BW

   ###
   # Data Conditioning
   ###

   # temporary fix, don't allow negative results if ReqNonNegResponse =1
   if( ReqNonNegResponse ==1) { y <- y*(y>=0) }


   MinNZResp <- min(y[y>0])    # find min NZ response
	ZeroSubstVal <- MinNZResp/2 # 0 response substitution value for certain situations

  # check for 0 results at NZ SL
	xnzsl <- x[x>0]  # get NZ spike values
	ynzsl <- y[x>0]  # get corresponsing measurements
	ZeroRespFlag <- any(ynzsl==0)
	if (ZeroRespFlag)
      {
		# check for fraction 0 response at NZ SL
		ZeroRespSL <- unique(xnzsl[ynzsl==0]) # get list of NZ SL with 0 response
		LowLimLCMRL <- min(x[x>max(ZeroRespSL)]) # smallest SL with no 0 responses
		LowLimDL <- LowLimLCMRL # since there should be no 0 response at DL
		nLvls <- length(ZeroRespSL) #count number of the SL
		FracNZResp <- vector("numeric",nLvls)
		for ( i in 1:nLvls )
         {
			ytmp <- ynzsl[xnzsl==ZeroRespSL[i]]
			FracNZResp[i]<- sum(ytmp>0)/length(ytmp)  # compute fraction NZ response
			if ( FracNZResp[i]>=0.5 ) # for NZ SL with 0 response fraction between 0.5 and 1
            {
				ytmp[ytmp==0] <- ZeroSubstVal; # replace 0 responses with half min NZ response
				ynzsl[xnzsl==ZeroRespSL[i]] <- ytmp
			   }
         }
		DnuSL <- ZeroRespSL[FracNZResp<0.5]
		if ( (length(unique(x)) - length(DnuSL)) < 4) # check to make sure there
                                                    # are enough remaining SL
                                                    # for valid estimates
         {
			LCMRLResultFlag <- -4;
			LCMRLResultMessage <- 'Aborted: Not enough spiking levels with nonzero
                                 results';
			DLResultFlag <- -4;
			DLResultMessage <- 'Aborted: Not enough spiking levels with nonzero
                              results';
         dat.tmp <-list(NA,NA,NA,LCMRLResultFlag,LCMRLResultMessage )
         names(dat.tmp) <- list('LCMRLval','DLval','CLval','LCMRLResultFlag',
                                 'LCMRLResultMessage')
         return(dat.tmp)
			}

		Nall0Resp <- sum(FracNZResp==0);
		if (Nall0Resp)
         {
			if (Nall0Resp > 1)
            {
				warning(' Nonzero spiking levels have all 0 responses')
			   } else
				{
            warning('One nonzero spiking level has all 0 responses')
            }
         warning(' This indicates a quality problem with the data.')
         }
      }

   n <- length(xorg) # number of original data points
   if(!is.null(DnuSL))
      {
      UseIndex <- !is.element(x, DnuSL) # find index of observations to use in models
	   x <- xorg[UseIndex] #subset x
      y <- yorg[UseIndex] # subset y
      BW <- BW[UseIndex]
      } else
      {
      x <- xorg
      y <- yorg
      BW <- BWorig
      }
   n <- length(x)
   normBW <- BW/sum(BW)
   xbar <- sum(normBW*x)
   #ssx <- n*sum(normBW*(x-xbar)^2)
   ssx <- n*sum(normBW*(x-xbar)^2)
   # plot raw data
   #plot(x,y, xlab=paste("SpikeConc"," ","(",strUnits,")"), ylab="MeasConc",
   #      main=AnalName)

   #compute replicate variances using robust w-estimator
   tmp <- RobRepVar(x,y,aPwts=BW)
   SpikeLevels <- tmp[[1]]
   RobVars <- tmp[[2]]
   RobLoc <- tmp[[3]]
   RobWt <- tmp[[4]]
   RepVarDoF <- tmp[[5]]
   nLevels <- tmp[[6]]
   remove(tmp)
   #plot(SpikeLevels,RobVars,xlab=paste("SpikeConc"," ","(",strUnits,")"),
   #     ylab="RobVars", main=AnalName)
   # assemble vector of robust weights, uses all data
   robWts <-unlist(RobWt)

   # check for 0 spiking level (method blanks) and remove for replicate
   # variance model calc
   if (SpikeLevels[1]==0 )
      {
      VarSpikeLevels <- SpikeLevels[2:nLevels]
      RobVars <- RobVars[2:nLevels]
      RepVarDoF <- RepVarDoF[2:nLevels]
      nvLevels <- nLevels-1
      } else
      {
      VarSpikeLevels <- SpikeLevels
      nvLevels <- nLevels
      }
   if (nvLevels < 4)
      {
		LCMRLResultFlag <- -4;
		LCMRLResultMessage <- 'Aborted: Not enough spiking levels with all nonzero results'
		warning(LCMRLResultMessage)
      DLResultFlag <- -4;
		DLResultMessage <- 'Aborted: Not enough spiking levels with all nonzero
                            results';
		dat.tmp <-list(NA,NA,NA,LCMRLResultFlag,LCMRLResultMessage )
         names(dat.tmp) <- list('LCMRLval','DLval','CLval','LCMRLResultFlag',
                                 'LCMRLResultMessage')
      return(dat.tmp)
      }

   #Check for 0 variances
   if (all(RobVars == 0 ) )
      {
      stop("All replicate variances are 0. Likely cause is data entry error")
      }
   if( any(RobVars == 0))
      {
      warning(paste("Warning: One or more replicate variances is 0. ",
                  "Likely causes are data entry error or variance damping (due",
                  " to insufficient digits, smoothing or thresholding"))
      VarOK <- (RobVars > 0)
      VarSpikeLevels <- VarSpikeLevels[VarOK]
      RobVars <- RobVars[VarOK]
      RepVarDoF <- RepVarDoF[VarOK]
      nvLevels <- length(VarSpikeLevels)
      }

   # estimate variance model

   Vtmp <- VarMod(VarSpikeLevels,RobVars,RepVarDoF)
   gc(verbose=FALSE)
   VarFnPar <- Vtmp

   # estimate means model
   m.mod <- IRLS(x,y,VarFnPar,robWts,ReqNonNegResponse,theta,maxIter,grphs,LN=labnum,BW)
   gc(verbose=FALSE)
   #print("########## end of IRLS")
   #print(m.mod)

   MeanFnPar <- m.mod[[1]]
   MSEFnPar <- m.mod[[2]]
   R2adj <-m.mod[[3]]
   Cp <- m.mod[[4]]
   pred <-m.mod[[5]]
   resids <-m.mod[[6]]
   residPearson <-m.mod[[7]]
   S_StdRes <-m.mod[[8]]
   cMSE <- m.mod[[9]]
   cmseSpikeLevels <- m.mod[[10]]
   Means.Model <- m.mod[[11]]

   # predicted MSE
   if( MSEFnPar$type == 'constant')
      {
      predMSE <- MSEFnPar$a*rep(1,nLevels)
      } else
      {
      predMSE <- pmax(MSEFnPar$a + MSEFnPar$b * SpikeLevels^MSEFnPar$c,
                    rep(MSEFnPar$minVar,nLevels))
      }

   # plot conditional MSE model

   # check for 0 spiking level (method blanks) and remove for LCMRL calc
   LCMRLSpikeLevels <- SpikeLevels
   if( SpikeLevels[1] == 0 )
      {
      LCMRLSpikeLevels <- SpikeLevels[2:nLevels]
      } else
      {
      LCMRLSpikeLevels <- SpikeLevels
      }

   # compute LCMRL
   LCMRL.tmp <-
       srchLCMRL(LCMRLSpikeLevels,lowQL,upQL,MeanFnPar,VarFnPar,MSEFnPar,
            CovProbReq,n,xbar,ssx,ReqNonNegResponse,LowLimLCMRL )
   LCMRLval <- LCMRL.tmp[[1]][[1]]
   CovProbMat <- LCMRL.tmp[[2]]
   LCMRLResultFlag <- LCMRL.tmp[[3]]
   LCMRLResultMessage <- LCMRL.tmp[[4]]

   # for MRL calculations return only LCMRLval & LCMRLResultflag
   if( !all(BW == 1) ) { return(c(LCMRLval,LCMRLResultFlag )) }

   # compute HV DL
   if(ZeroRespFlag)
      {
      HVDL.lst <- HubauxVos(MeanFnPar,VarFnPar,MSEFnPar,alpha,beta,
                         LCMRLval,LCMRLSpikeLevels,ReqNonNegResponse,ZeroRespFlag)
      } else
      {
      HVDL.lst <- HubauxVos(MeanFnPar,VarFnPar,MSEFnPar,alpha,beta,
                         LCMRLval,SpikeLevels,ReqNonNegResponse,ZeroRespFlag)
      }

   DLval <- HVDL.lst[[1]]
   CLval <- HVDL.lst[[2]]
   HVResultFlag <- HVDL.lst[[3]]
   HVResultMessage <- HVDL.lst[[4]]

   # GRAPHICAL OUTPUT
   if( grphs > 0  && (LCMRLResultFlag==1 || LCMRLResultFlag==-1) )
      {
      #par(ask=TRUE)
      xpts <- CovProbMat[,1]
      ypts <- CovProbMat[,2]
      npts <- length(xpts)
      dof <- MSEFnPar$DoF

      # create mean and variance vectors
      mu <- MeanFn(xpts,MeanFnPar)
      v <- VarFn(xpts,MSEFnPar)*(1+1/n+((xpts-xbar)^2)/ssx)

      # probability levels for upper and lower prediction
      lcmrlLowProb <- (1-CovProbReq)/2
      lcmrlUpProb <- 1-lcmrlLowProb
      lowPred <- vector("numeric",npts)
      upPred <- vector("numeric",npts)

      # build vectors for moment inequalty check
      std <- sqrt(v)
      M <- mu +3*std

      if(ReqNonNegResponse == 1)
         {
         for( i in 1:npts)
            {
            if( v[i] <= mu[i]*(M[i]-mu[i]) ) # if moment ineq. OK use gamma
               {
               gpars <- GammaPars(mu[i],v[i])
               lowPred[i] <- qgamma(lcmrlLowProb,gpars[[1]],scale=gpars[[2]])
               upPred[i] <- qgamma(lcmrlUpProb,gpars[[1]],scale=gpars[[2]])
               } else     # otherwise use truncted t
               {
               z <- qt(pt(-mu[i]/std[i],dof) + lcmrlLowProb*pt(-mu[i]/std[i],dof))
               lowPred[i] <- mu[i]+z*std[i]
               z <- qt(pt(-mu[i]/std[i],dof) + lcmrlUpProb*pt(-mu[i]/std[i],dof))
               upPred[i] <- mu[i]+z*std[i]
               }
            }
         } else
         {
         lowPred <- mu + qt(lcmrlLowProb,dof)*sqrt(v)
         upPred  <- mu + qt(lcmrlUpProb, dof)*sqrt(v)
         }
     maxS <- max(SpikeConc)
     MinCov <- min(CovProbMat[,2],CovProbReq)
     xlow <-  min(xpts)/2

     # plot LCMRL in terms of QC interval coverage probability

      if(is.null(labnum))
         {
         plot(xpts,ypts,main=paste(AnalName,"QC Interval Coverage Plot"),
         xlab=paste("True Concentration"," (", strUnits, ")", sep = ''),
         ylab="Probability of QC Interval Coverage",
         xlim=c(0,maxS),ylim=c(0.96*MinCov, 1.0),
         pch="", ask=TRUE)
         } else
         {
         plot(xpts,ypts,main=paste(AnalName,"Lab",
                     labnum,"\n","QC Interval Coverage Plot"),
         xlab=paste("True Concentration"," (", strUnits, ")", sep = ''),
         ylab="Probability of QC Interval Coverage",
         xlim=c(0,maxS),ylim=c(0.96*MinCov, 1.0),
         pch="", ask=TRUE)
         }
      lines(xpts,ypts,col='blue')
      lines(c(LCMRLval,LCMRLval),c(0.96*MinCov,CovProbReq),lty=2)
      lines(c(0,maxS),c(CovProbReq,CovProbReq),lty=2, col='red')
      points(LCMRLval,CovProbReq, pch='+',col='green',cex=1.5)
      if(LCMRLResultFlag==1)
         {
         legend(1.5*LCMRLval,CovProbReq,
         c(paste('LCMRL =',sprintf("%.3f",LCMRLval),strUnits),
           "Qual Lim: 50-150%", "Coverage Prob: 0.99"),
         bty='n')
         } else
         {
         op <- par(bg="grey81")
         legend(1.5*LCMRLval,CovProbReq,
         c(paste('LCMRL =',sprintf("%.3f",LCMRLval),strUnits),
           "Qual Lim: 50-150%", "Coverage Prob: 0.99" ))
         par(op)
         legend("center",legend="LCMRL is Below Lowest Non-Zero SL",text.col=2,bty='n')
         }
      # plot LCMRL, DL, response model, prediction limits and quality limits
      if(is.null(labnum))
         {
         plot(SpikeConc,MeasConc, main=paste(AnalName,"- LCMRL Plot"),
         xlab=paste("True Concentration", " (", strUnits, ")", sep = ''),
         ylab=paste("Measured Concentration", " (", strUnits, ")", sep = ''),
                  col='blue',pch="o",
         ylim=c(0,1.05*upQL*maxS),xlim=c(0,maxS))
         } else
         {
         plot(SpikeConc,MeasConc, main=paste(AnalName,"Lab",
                                 labnum,"\n","LCMRL Plot"),
         xlab=paste("True Concentration", " (", strUnits, ")", sep = ''),
         ylab=paste("Measured Concentration", " (", strUnits, ")", sep = ''),
                     col='blue',pch="o",
         ylim=c(0,1.05*upQL*maxS),xlim=c(0,maxS))
         }
      points(LCMRLval,LCMRLval, pch='+',cex=1.5, col='green')
      points(DLval,DLval, pch='x',cex=1,col='green' )
      lines(xpts,mu)
      lines(c(xlow,maxS),c(xlow,maxS)*lowQL, lty=2,col='red')
      lines(c(xlow,maxS),c(xlow,maxS)*upQL, lty=2, col='red')
      lines(xpts,lowPred,lty=3,col='blue')
      lines(xpts,upPred,lty=3,col='blue')
      topleft <- par()$usr[c(1,4)]
      if(LCMRLResultFlag==1)
         {
         legend(topleft[1],topleft[2],
         c('Data',paste('LCMRL =',sprintf("%.3f",LCMRLval),strUnits),
         paste('Hubaux-Vos DL',sprintf("%.3f",DLval),strUnits),
         '50-150% Recovery','Lower/Upper Prediction Limits'),
         col=c('blue','green','green','red','blue'),
         pch=c("o","+","x","",""),lty=c(0,0,0,2,3) )
         } else
         {
         op <- par(bg="grey81")
         legend(topleft[1],topleft[2],
         c('Data',paste('LCMRL =',sprintf("%.3f",LCMRLval),strUnits),
         paste('Hubaux-Vos DL',sprintf("%.3f",DLval),strUnits),
         '50-150% Recovery','Lower/Upper Prediction Limits'),
         col=c('blue','green','green','red','blue'),
         pch=c("o","+","x","",""),lty=c(0,0,0,2,3))
         par(op)
         legend("center",legend='LCRML is Below Lowest Non-Zero SL',text.col=2,bty='n')
         }
      par(ask=FALSE)
      }
      wts <- cbind(x,y,robWts)
      dimnames(wts)[[2]] <-  c("SL","Measurement","RobWts")
      dat.tmp <-list(LCMRLval,DLval,CLval,
                  LCMRLResultFlag,LCMRLResultMessage,HVResultFlag,
                  HVResultMessage,MeanFnPar,VarFnPar,MSEFnPar)
      names(dat.tmp) <- list('LCMRLval','DLval','CLval',
        'LCMRLResultFlag','LCMRLResultMessage','HVResultFlag','HVResultMessage',
        'MeanFnPar','VarFnPar','MSEFnPar')
      return(dat.tmp)
   }

#' RobRepVar
#'
#' @param X X
#' @param Y Y
#' @param cc cc
#' @param k k
#' @param aPwts aPwts
#'
#' @export
#'
RobRepVar <- function (X, Y, cc=9, k=1, aPwts)

# computes replicate variances at various spiking levels,
# degrees of freedom and weights for modeling
#
# ----------------Outputs-------------------------------------
#
# SpikeLevels = vector of spiking levels with 2 or more results
# RepVars = vector of computed replicate variance at each spiking level
# RepVarDoF = vector of degrees of freedom for each element of RepVars
# RepVarWt = vector of weights for modeling mean-variance function
#
# ----------------Inputs-------------------------------------
#
# X = vector of spiked concentrations
# Y = vector of measured concentrations
# c= tuning constant for biweight
# k = threshold for Huber estimator
# aPwys = normalized Bayesian Weights for MRL calculations#
# ----------------Start code---------------------------------
#
  {
  #print("#### In RobRepVar")
  tol12 <-1e-12
  # get unique values of X & counts
  ST <- cbind(as.numeric(names(table(X))),table(X))
  # only retain those for which count > 1
  SpkTmp <- ST[ST[,2]>1,c(1,2)]
  # first column is spiking level
  # second column is count; DoF is count - 1
  SL <- SpkTmp[,1]
  # number of spiking levels
  nLs <- length(SL)
  RVars <- rep(0,nLs)
  RLoc <- RVars
  RWt <- vector("list",nLs)
  RVarDoF <- vector('numeric',length(SL))
  # compute replicate variances and their associated weights
  for(i in 1:nLs)
    {
    #print(paste("SpikeLevel:",i))
    yi <- Y[ abs(X-SL[i]) < tol12]
    aPwtsi <-aPwts[abs(X-SL[i]) < tol12]
    if( var(yi) < tol12 )
      {
      ni <- length(yi)
      RLoc[i] <- yi[1]
      RVars[i] <- 0
      RWt[[i]] <- rep(1,ni)/ni
      RWt[[i]] <- RWt[[i]]*aPwtsi
      RWt[[i]] <- RWt[[i]]/sum(RWt[[i]])
      RVarDoF[i] <- ni*(1-sum(RWt[[i]]^2)) # will be (ni-1) if aPwts are all 1
      #print(c(yi,aPwtsi,ni,RLoc[i],RVars[i],RWt[[i]],RVarDoF[i]))
      } else
      {
      zz <- roblocvar4(yi,cc,k,0,aPwtsi)
      RLoc[i] <-  zz[[1]]
      RVars[i] <- zz[[2]]
      RWt[[i]] <- zz[[3]]
      RVarDoF[i] <- zz[[4]]
      #print(c(yi,aPwtsi,RLoc[i],RVars[i],RWt[[i]],RVarDoF[i]))
      }
    }
    return(list(SL, RVars, RLoc, RWt, RVarDoF, nLs))
   }

#' roblocvar4
#'
#' @param XX XX
#' @param cc cc
#' @param k k
#' @param convtype conversion type
#' @param Pwts Pwts
#'
#' @export
#'
#'
roblocvar4 <- function(XX, cc=9, k=1, convtype=0, Pwts=rep(1,length(XX)))
   {
   # roblocvar computes M-estimators of location and also returns weight vector
   #
   # ----------------Inputs-------------------------------------
   #
   # XX = numeric data vector
   # cc = parameter for Bi-weight
   # k = parameter for Huber weight
   #
   # ---------------Outputs-------------------------------------
   #
   # rloc - M-estimator of location
   # rvar - M-estimator of scale
   # rwts - vector of M-estimator of location weights at convergence
   # nw -   weighted number of observations :(nw-1) is dgrees of freedom
   # IsError = Error flag (0-1)
   # Err = Matlab error structure
   #
   # ----------------Start main code---------------------------------
   #
   # initialize variables
   #
   ##print("In roblocvar4")
   ##print(paste("ConvType:",convtype))
   ##print(paste("XX:",XX))
   ##print(paste("aPwts:",Pwts))
   tol <- 1e-4
   maxsteps <- 10
   nsteps <- 0

   # use modified Hodges-Lehmann estimator for inital location estimate
   # and Average Deviation from HL location estimate

   n <- length(XX)
   m <- n*(n-1)/2

   pair_means <- rep(0,m)
   i <- 0
   for (rowi in 1:(n-1))
     {
     for ( colj in (rowi+1):n )
       {
       i <- i+1
       pair_means[i] <- (XX[rowi] + XX[colj])/2
       }
     }

   rloc <- median(c(pair_means, median(XX))) # mod. Hodges-Lehmann est
   scal <- 1.4826*mean(abs(XX-rloc))
   ##print(paste("intial.roloc:",rloc))
   ##print(paste("intial.scal:",scal))

   # compute Huber(1.0) M-estimate of location and associated weights
   rl_old <- rloc
   delta <- 1
   while ((delta > tol) && (nsteps <= maxsteps))
     {
     ##print("..In Huber while loop")
     rwtsH <- HuberWt(XX,rloc,scal,k)
     ##print(paste("Huber Wt",rwtsH))
     rwtsH <- BlendWts(Pwts,rwtsH)
     ##print(paste("Blend Wt",rwtsH))
     rloc <- sum(rwtsH*XX)
     if(convtype == 0)
      {
      if(rl_old != 0)
         {
         delta <- try(abs((rl_old - rloc)/rl_old),silent==TRUE)
         #if(class(delta)=="try-error" ) {delta <- 0 }
         }
      }else
      {
      delta <- try(abs((rl_old - rloc)/scal),silent==TRUE)
      #if(class(delta)=="try-error" ) {delta <- 0 }
      }
     nsteps <- nsteps + 1
     rl_old <- rloc
     ##print(paste("delta:",delta,"nsteps",nsteps))
     }
   nw <- n*(1-sum(rwtsH*rwtsH))
   rvar <- (n/nw)*sum(rwtsH*(XX-rloc)^2)
   scal <- sqrt(rvar)

   #compute Biweight(9) M-estimate of location and associated weights
   rl_old <- rloc
   delta <- 1
   nsteps <- 0
   while (delta > tol && nsteps <= maxsteps)
     {
     rwtsBW <- BisquareWt(XX,rloc,scal,cc)
     rwtsBW <- BlendWts(Pwts,rwtsBW)
     rloc <- sum(rwtsBW*XX)
     if(convtype == 0)
      {
      if( rl_old != 0) {delta <- try(abs((rl_old - rloc)/rl_old),silent=TRUE)}
      if( inherits(delta, "try-error") ) {delta <- 0 }
      } else
      {
      delta <- try(abs((rl_old - rloc)/scal),silent=TRUE)
      if( inherits(delta, "try-error") ) {delta <- 0 }
      }
     if(is.nan(delta)) { delta <- 0}
     nsteps <- nsteps+1
     rl_old <- rloc
     }

   # average weights from Huber(1.0) and biweight(9) estimates
   # compute robust weighted mean and variance with the average weights
   rloc <- sum(rwtsBW*XX)
   nw <- n*(1-sum(rwtsBW*rwtsBW))
   rvar <- (n/nw)*sum(rwtsBW*(XX-rloc)^2)
   return(list(rloc,rvar,rwtsBW,nw))
   }


#' HuberWt
#'
#' @param x data
#' @param m mean
#' @param s standard deviation
#' @param k k
#'
#' @export
#'
HuberWt <- function(x, m, s, k)
   {
   u <- (x-m)/s
   nu <- length(u)
   tu <- (abs(u) <= k)   # k=1
   rwts <- tu
   for(i in 1:nu)
     {
     if(tu[i] == 0) { rwts[i] <- k/abs(u[i]) }
     }
   return(rwts/sum(rwts))
   }


#' BlendWts
#'
#' @param wt1 weight 1
#' @param wt2 weight 2
#'
#' @export
#'
BlendWts <- function(wt1, wt2)
   {
   wt <- wt1*wt2
   wt <- wt/sum(wt)
   return(wt)
   }

#' BisquareWt
#'
#' @param x data
#' @param m mean
#' @param scal scal
#' @param cc cc
#'
#' @export
#'
BisquareWt <- function(x, m, scal, cc)
   {
   u <- (x - m)/(cc*scal)   #c=9  tunig constant
   tu <- (abs(u)<=1)
   rwts <- ((1-u^2)^2)*tu
   rwts <- rwts/sum(rwts)
   return(rwts)
   }


#' x.var
#'
#' @param xx data
#' @param wts weights
#'
#' @export
#'
w.var <- function(xx, wts)
   {
   n <- length(xx)
   normwts <- wts/sum(wts)
   m <- sum(normwts*xx)
   vx <- sum(normwts * (xx - m)^2)/(1-sum(normwts*normwts))
   return(vx)
   }

#' VarFn
#'
#' @param xv xv
#' @param VFnP VFnP
#'
#' @export
#'
VarFn <- function (xv,VFnP)
  {
  epsilon <- 1
  xv <- xv*(xv > 0)
  switch(EXPR =VFnP$type,
      power = {v <-pmax((VFnP$b * xv^VFnP$c),rep(VFnP$minVar,length(xv)))},
      constant.power = {v <- VFnP$a + VFnP$b*xv^VFnP$c },
      constant = { v <- VFnP$a * rep(1,length(xv)) }
      )
  return(v)
  }

#' MeanFn
#'
#' @param xm mean
#' @param MFP MFP
#'
#' @export
#'
MeanFn <- function(xm, MFP)
  {
  switch(MFP$type,
    'linear'=
      {mu <- MFP$a + MFP$b*xm},
    'quadratic'=
      { mu <- MFP$a + MFP$b*xm + MFP$c*xm^2},
    'cubic'=
      {mu <-MFP$a + MFP$b*xm + MFP$c*xm^2 + MFP$d*xm^3})
  # do not allow a negative mean response
  minv <- max(0,MFP$a)
  mu<- pmax(rep(minv,length(mu)),mu)
  return(mu)
  }

#' VarMod
#'
#' @param SL SL
#' @param RV RV
#' @param RVDF RVDF
#'
#' @export
#'
VarMod <- function (SL, RV, RVDF)

   # VarModdev function
   #
   # computes Gamma nl model for replicate variance at various spiking
   # levels
   #
   # ----------------Inputs-------------------------------------
   #
   # SpikeLevels = vector of spiking levels with 2 or more results
   # RepVars = vector of computed replicate variance at each spiking level
   # RepVarDoF = vector of degrees of freedom for each element of RepVars
   #
   # ---------------Outputs-------------------------------------
   #
   # VFnP - structure for error variance model, type and parameters
   # ------------------------------------------------------------
   #
   {
   maxC <- 2
   VFnP <- list("",NaN,NaN,NaN,NaN,NaN)
   names(VFnP) <- list("type","a","b","c","DoF","minVar")
   nL <- length(SL)
   DoF <- sum(RVDF)

   nL <- length(SL)
   DoF <- sum(RVDF)

   # Estimate initial values for b and c using log-scale regression
   # make design matrix omiting first spiking level
   lBC<-lm(log(RV[2:nL])~log(SL[2:nL]),weights=RVDF[2:nL])

   # estimate initial values for a, b and c
   ndx <- 1:floor((nL/2)-1)
   A <- max(sum(RV[ndx]*RVDF[ndx])/sum(RVDF[ndx]),1e-8)
   B <- exp(lBC$coefficients[1])
   CC <-min(max(0, lBC$coefficients[2]),maxC)

   newCoef <- OptABC(A,B,CC,SL,RV,RVDF)
   A <- newCoef[1]; names(A) <- 'constant'
   B <- newCoef[2]; names(B) <- 'slope'
   CC <- newCoef[3]; names(CC) <- 'power'

   if( B <= 0 || CC <= 1e-2 || (B*max(SL)^CC < 0.10*A))
      {
      VFnP$type<-'constant'
      VFnP$a <- mean(RV)
      VFnP$b <- 0
      VFnP$c <- 0
      VFnP$minVar <- A
      VFnP$minVar <- VFnP$a
      VFnP$DoF <- DoF
      }else
      {
      if ( A < (1e-6)*mean(RV))
         {
         VFnP$type <- 'power'
         VFnP$a <- 0
         VFnP$b <- B
         VFnP$c <- CC
         VFnP$minVar <- mean(RV[1:2])
         VFnP$DoF <- DoF - 2
         } else
         {
         VFnP$type <- 'constant.power'
         VFnP$a <- A
         VFnP$b <- B
         VFnP$c <- CC
         VFnP$minVar <- A
         VFnP$DoF <- DoF - 3
         }
      }
   return(VFnP)
   }


#' VfnLoss
#'
#' @param x data
#' @param S S
#' @param V V
#' @param DF DF
#'
#' @export
#'
VfnLoss <-function(x, S, V, DF)
   {
   maxC <- 2
   BigLoss <- 1e12
   d <- rep(0,length(S))
   e <- d
   for(i in 1:length(S))
      {
      d[i] <- x[1] + x[2]*S[i]^x[3]
      ifelse(d[i] > 0,e[i] <- DF[i]*((V[i]-d[i])^2)/d[i], e[i] <-0)
      }
   loss <- sum(e)
   # ensure parameters in acceptable ranges
   if( x[1] < 1e-08 || x[2] < 0 || x[3] <0 || x[3] > maxC )
      {
      loss <- BigLoss
      }
   return(loss)
   }

#' OptABC
#'
#' @param A A
#' @param B B
#' @param CC CC
#' @param SL SL
#' @param RV RV
#' @param RVDF RVDF
#'
#' @export
#'
OptABC <-function (A, B, CC, SL, RV, RVDF)
   {

   # initialize variables
   maxC <- 2
   delta <- -1
   maxiter <- 5
   eps <- 1e-4

   x <-c(A,B,CC)
   rslt <- optim(x,VfnLoss,control=list(abstol=1e-16,reltol=1e-16,maxit=10000),S=SL,V=RV,DF=RVDF)
   fold <- rslt$value
   iter <- 1

   while ( (delta <= -eps ) && (iter < maxiter) )
      {
      rslt <- optim(rslt[[1]],VfnLoss,control=list(abstol=1e-16,reltol=1e-16,maxit=10000),S=SL,V=RV,DF=RVDF)
      delta <- (rslt$value - fold)/rslt$value
      fold <- rslt$value
      iter <- iter + 1
      }

   A <- rslt[[1]][1]
   if (A < 0 ) A <- 0
   B <- rslt[[1]][2]
   if (B < 0) B <- 0
   CC <- rslt[[1]][3]
   if (CC < 0 ) CC <- 0
   if (CC > maxC)  CC <- maxC
   return(c(A,B,CC))
   }

#' IRLS
#'
#' @param xx xx
#' @param yy yy
#' @param VFnP VFnP
#' @param Wts weights
#' @param RNNR RNNR
#' @param theta theta
#' @param maxIter max iteration
#' @param grphs graphs
#' @param LN NL
#' @param BW BW
#'
#' @export
#'
IRLS <- function(xx, yy, VFnP, Wts, RNNR, theta, maxIter, grphs=1, LN=NULL, BW)
   # WLS computes Iteratively Reweighted Least Squares regressions estimates
   # covB using variance and MSE function models
   #
   # ----------------Inputs-------------------------------------
   #
   # x = column vector of Spiking concentrations
   # y = column vector of Measured response concentrations
   # VFPar = structure of parameters for variance function fit in VarModdev
   # initWts = vector of weights from robust location-scale estmation
   # RNNR = require mean function to be non-negative over range
   # theta = tolerance for convergence
   # maxIter = maximum allowed number of iterations
   #
   # ----------------Outputs-------------------------------------
   #
   # MeanFnPar = structure of parameters for mean concentration function
   # MSEfnPar = structure of parameters for conditional MSE function
   # pred = vector of predictions corresponding to input vector x
   # resids = vector of raw residuals
   # residPearson = vector of Pearson residuals
   # S_StdRes = std dev of standardized (by variance function) residuals
   # IsError = error flag
   # Err = error messsage structure
   #
   # ----------------Start code---------------------------------
   #
   # initialize variables
   #
   {
   k <- 1   # threshold for Huber estimator
   cc <- 9  # initial tuing constant for biweight estimator
   converged <- 1


   #print(paste("###### start IRLS"))
   #print(paste("########### cc =",cc))
   ############

   n <- length(xx) # number of observatons
   MFnP <- list("",0,0,0,0,0,0)
   names(MFnP) <- list("type","a","b","c","d","npar","DoF")
   MSEFnP <-VFnP

   # make design matrices for linear, quadratic, cubic & quartic models
   DX <- cbind(xx)
   DX2 <- cbind(DX,xx^2)
   DX3 <- cbind(DX2,xx^3)
   DX4 <- cbind(DX3,xx^4)
   if(RNNR ==1)
      {
      xmx <- max(xx)
      xgrid <- seq(0,xmx,length.out = 100)
      }
   # compute initial 1-step BW w-estimator using WLS for linear model
   # start with inital WLS estimate

   # initial 1-step BW w-estimator using WLS
   oldMSEfnPar1 <- VFnP  # intial value for MSEfnPar
   oldB1 <- lm(yy~xx,weights = Wts)$coefficients # initial seed est. for WLS
   oldB1[is.na(oldB1)] <- 0
   #[oldB1,oldseB1,oldmse1,oldcovB1,oldp1,pred] = ...
   #  compWLS(xx,yy,VFnP,oldMSEfnPar1,oldB1,cc,ReqNonNegResponse);
   # compWLS returns a list:
   #     [[1]] with an lm.wfit object
   #     {[2]] a numeric giving the residual mean square error of the fit
   #     [[3]] a numeric p giving the number of fit parameters

   wls.tmp <- compWLS(xx,DX,yy,VFnP,oldMSEfnPar1,oldB1,cc,RNNR,BW)
   old.mean.model.1 <- wls.tmp[[1]]
   oldB1 <-wls.tmp[[4]]
   oldB1[is.na(oldB1)] <- 0
   #oldseB1 <- summary(old.mean.model.1 )$coefficients[,2]
   oldmse1 <-  wls.tmp[[2]]
   #oldcovB1 <- summary(old.mean.model.1 )$cov
   oldp1 <- wls.tmp[[3]]
   resids <- yy - old.mean.model.1$fitted	# compute raw residuals
   oldnw1 <- wls.tmp[[5]]
   #compute MSE fn
   # compMSE returns a list with MSEFnPar, cMSE estimates, cMSE Spiking Levels
   #         the MSE model object
   mse.t <- compMSE(xx,resids,cc,k,VFnP,BW)
   oldMSEfnPar1 <- mse.t[[1]]
   oldcMSE1 <- mse.t[[2]]
   cmseSL <- mse.t[[3]]

   # compute initial 1-step BW w-estimator using WLS for quadratic model
   # start with inital WLS estimate
   oldMSEfnPar2 <- oldMSEfnPar1 # use MSE fn from linear as initial value
   oldB2 <- lm(yy~DX2,weights = Wts)$coefficients # initial seed est. for WLS
   oldB2[is.na(oldB2)] <- 0
   # initial 1-step BW w-estimator using WLS
   #[oldB2,oldseB2,oldmse2,oldcovB2,oldp2,pred] =  ...
   # compWLS(x,DX2,y,VarFnPar,oldMSEfnPar2,oldB2,c,ReqNonNegResponse);
   wls.tmp <- compWLS(xx,DX2,yy,VFnP,oldMSEfnPar2,oldB2,cc,RNNR,BW)
   old.mean.model.2 <- wls.tmp[[1]]
   oldB2 <-wls.tmp[[4]]
   oldB2[is.na(oldB2)] <- 0
   #oldseB2 <- summary(old.mean.model.2 )$coefficients[,2]
   oldmse2 <-  wls.tmp[[2]]
   #oldcovB2 <- summary(old.mean.model.2)$cov
   oldp2 <- wls.tmp[[3]]
   resids <- yy - old.mean.model.2$fitted;	# compute raw residuals
   oldnw2 <- wls.tmp[[5]]
   mse.t <- compMSE(xx,resids,cc,k,oldMSEfnPar1,BW); #compute MSE fn
   oldMSEfnPar2 <- mse.t[[1]]
   oldcMSE2 <- mse.t[[2]]

   # compute initial 1-step BW w-estimator using WLS for cubic model
   # start with inital WLS estimate
   oldMSEfnPar3 <- oldMSEfnPar2 # use MSE fn from linear as initial value
   oldB3 <- lm(yy~DX3,weights = Wts)$coefficients # initial seed est. for WLS
   oldB3[is.na(oldB3)] <- 0
   # initial 1-step BW w-estimator using WLS
   #[oldB3,oldseB3,oldmse3,oldcovB3,oldp3,pred] = ...
   #  compWLS(x,DX3,y,VarFnPar,oldMSEfnPar3,oldB3,c,ReqNonNegResponse);
   wls.tmp <- compWLS(xx,DX3,yy,VFnP,oldMSEfnPar3,oldB3,cc,RNNR,BW)
   old.mean.model.3 <- wls.tmp[[1]]
   oldB3 <-wls.tmp[[4]]
   oldB3[is.na(oldB3)] <- 0
   #oldseB3 <- summary(old.mean.model.3 )$coefficients[,2]
   oldmse3 <-  wls.tmp[[2]]
   #oldcovB3 <- summary(old.mean.model.3)$cov
   oldp3 <- wls.tmp[[3]]
   resids <- yy - old.mean.model.3$fitted	# compute raw residuals
   oldnw3 <- wls.tmp[[5]]
   mse.t <- compMSE(xx,resids,cc,k,oldMSEfnPar2,BW); #compute MSE fn
   oldMSEfnPar3 <- mse.t[[1]]
   oldcMSE3 <- mse.t[[2]]

   #print("####### end of first WLS call")
   #print(t(DX3))
   #print(oldB3)

   # compute initial 1-step BW w-estimator using WLS for quartic model
   # start with inital WLS estimate
   oldMSEfnPar4 <- oldMSEfnPar3 # use MSE fn from linear as initial value
   oldB4 <- lm(yy~DX4,weights = Wts)$coefficients # initial seed est. for WLS
   oldB4[is.na(oldB4)] <- 0
   # initial 1-step BW w-estimator using WLS
   #[oldB4,oldseB4,oldmse4] = ...
   #  compWLS(x,DX4,y,VarFnPar,oldMSEfnPar4,oldB4,c,ReqNonNegResponse);
   wls.tmp <- compWLS(xx,DX4,yy,VFnP,oldMSEfnPar2,oldB4,cc,RNNR,BW)
   old.mean.model.4 <- wls.tmp[[1]]
   oldB4 <-wls.tmp[[4]]
   oldB4[is.na(oldB4)] <- 0
   #oldseB4 <- summary(old.mean.model.4 )$coefficients[,2]
   oldmse4 <-  wls.tmp[[2]]
   #oldcovB4 <- summary(old.mean.model.4)$cov
   oldp4 <- wls.tmp[[3]]
   # finished initial 1-step BW w-estimators (c=9)


   # initially try IRLS with c=9
   # if any model fails to converge, go to c=9 for all models
   cc <- 9; # initial value for Tukey's biweight (BW) parameter


   all.irls <- compAllIRLS(xx,DX,DX2,DX3,DX4,yy,VFnP,oldMSEfnPar1,oldB1,
                           oldMSEfnPar2,oldB2,oldMSEfnPar3,oldB3,
                           oldMSEfnPar4,oldB4,
                           cc,k,RNNR,theta,maxIter,grphs,LN,BW)
   converged <- all.irls[[1]]
   #print("end of compALLIRLS")
   if( converged==0 )
      # check for convergence of all models
      # if any fail to converge, revert to initial 1-step estimate
      {
      B1<-oldB1
      mse1<-oldmse1
      p1<-oldp1
      MSEfnPar1 <- oldMSEfnPar1
      cMSE1 <- oldcMSE1
      nw1 <- oldnw1
      B2<-oldB2
      mse2<-oldmse2
      p2<-oldp2
      MSEfnPar2 <- oldMSEfnPar2
      cMSE2 <- oldcMSE2
      nw2 <- oldnw2
      B3<-oldB3
      mse3<-oldmse3
      p3<-oldp3
      MSEfnPar3 <- oldMSEfnPar3
      cMSE3 <- oldcMSE3
      nw3 <- oldnw3
      mse4<-oldmse4
      mean.model.1 <- old.mean.model.1
      mean.model.2 <- old.mean.model.2
      mean.model.3 <- old.mean.model.3
      mean.model.4 <- old.mean.model.4
      #cmse.model.4 <- old.cmse.model.4

      } else
      {
      B1<-all.irls[[2]]
      mse1<-all.irls[[3]]
      p1<-all.irls[[4]]
      MSEfnPar1 <- all.irls[[5]]
      cMSE1 <- all.irls[[6]]
      nw1 <- all.irls[[23]]
      B2<-all.irls[[7]]
      mse2<-all.irls[[8]]
      p2<-all.irls[[9]]
      MSEfnPar2 <- all.irls[[10]]
      cMSE2 <- all.irls[[11]]
      nw2 <- all.irls[[24]]
      B3<-all.irls[[12]]
      mse3<-all.irls[[13]]
      p3<-all.irls[[14]]
      MSEfnPar3 <- all.irls[[15]]
      cMSE3 <- all.irls[[16]]
      nw3 <- all.irls[[25]]
      mse4<-all.irls[[18]]
      mean.model.1 <- all.irls[[19]]
      mean.model.2 <- all.irls[[20]]
      mean.model.3 <- all.irls[[21]]
      mean.model.4 <- all.irls[[22]]
      }

   # compute DoF and RSS for linear through cubic models
   d1 <- nw1 - p1
   RSS1 <- mse1*d1
   d2 <- nw2 - p2
   RSS2 <- mse2*d2
   d3 <- nw3 - p3
   RSS3 <- mse3*d3

   # compute Cp values for different models
   Cp1 <- RSS1/mse4 - (d1-p1)
   Cp2 <- RSS2/mse4 - (d2-p2)
   Cp3 <- RSS3/mse4 - (d3-p3)
   Cp <- c(Cp1,Cp2,Cp3)

   minCp <- min(Cp)
   ndx <- 1:3
   minCpndx <- ndx[Cp == minCp]

   # set values to return based on best fitting model
   switch (minCpndx,
      { # case 1 % linear
      pred <- mean.model.1$fitted
      MFnP$type <- 'linear'
      MFnP$a <- B1[1]
      MFnP$b <- B1[2]
      MFnP$npar <- p1
      MFnP$DoF <- d1
      MSEFnP <- MSEfnPar1
      cMSE <- oldcMSE1
      M.model <- mean.model.1
      },
      { #case 2 % quadratic
      pred <- mean.model.2$fitted
      MFnP$type <- 'quadratic'
      MFnP$a <- B2[1]
      MFnP$b <- B2[2]
      MFnP$c <- B2[3]
      MFnP$npar <- p2
      MFnP$DoF <- d2
      MSEFnP <- MSEfnPar2
      cMSE <- oldcMSE2
      M.model <- mean.model.2
      },
      { # case 3 % cubic
      pred <- mean.model.3$fitted
      MFnP$type <- 'cubic'
      MFnP$a <- B3[1]
      MFnP$b <- B3[2]
      MFnP$c <- B3[3]
      MFnP$d <- B3[4]
      MFnP$npar <- p3
      MFnP$DoF <- d3
      MSEFnP <- MSEfnPar3
      cMSE <- oldcMSE3
      M.model <- mean.model.3
      }
      )

   resids <- yy - pred;
   vf <- VarFn(xx,MSEFnP) # variances
   wt <- Wts/vf
   wt <- wt/sum(wt)
   wtm <- sum(yy*wt)/sum(wt)
   dev <- yy - wtm
   SSrat <- sum(wt*resids^2)/sum(wt*dev^2)
   R2adj <- 1 - SSrat*(MFnP$DoF+MFnP$npar-1)/MFnP$DoF
   StdRes <- resids/sqrt(vf)
   S_StdRes <- sd(StdRes)
   residPearson <- StdRes/S_StdRes
   #print("########end IRLS")
   return(list(MFnP,MSEFnP,R2adj,Cp,pred,resids,residPearson,S_StdRes,cMSE,
          cmseSL,M.model))
   }

#' compAllIRLS
#'
#' @param x x
#' @param DX DX
#' @param DX2 DX2
#' @param DX3 DX3
#' @param DX4 DX4
#' @param y y
#' @param VFnP VFnP
#' @param oldMSEfnPar1 oldMSEfnPar1
#' @param oldB1 oldB1
#' @param oldMSEfnPar2 oldMSEfnPar2
#' @param oldB2 oldB2
#' @param oldMSEfnPar3 oldMSEfnPar3
#' @param oldB3 oldB3
#' @param oldMSEfnPar4 oldMSEfnPar4
#' @param oldB4 oldB4
#' @param cc cc
#' @param k k
#' @param RNNR RNNR
#' @param theta theta
#' @param maxIter max iteration
#' @param grphs graphs
#' @param LN LN
#' @param BW BW
#'
#' @export
#'
#'
compAllIRLS <- function(x, DX, DX2, DX3, DX4, y, VFnP,
                        oldMSEfnPar1, oldB1, oldMSEfnPar2, oldB2,
                        oldMSEfnPar3, oldB3, oldMSEfnPar4, oldB4,
                        cc, k, RNNR, theta, maxIter, grphs, LN, BW)
   {
   converged <- 1


   # create grid to test for negative values
   if( RNNR==1 )
      {
      xmx <- max(x)
      xgrid <- seq(0,xmx,length.out=100)
      DXgrid <- cbind(rep(1,100),xgrid)
      }

   # compute IRLS for linear model
   irls.tmp <- compIRLS(x,DX,y,VFnP,oldMSEfnPar1,oldB1,cc,k,RNNR,
                        theta,maxIter,grphs,LN,BW)
   B1 <- irls.tmp[[1]]
   #seB1 <- irls.tmp[[2]]
   mse1 <- irls.tmp[[2]]
   #covB1 <- irls.tmp[[4]]
   p1 <- irls.tmp[[3]]
   MSEfnPar1 <- irls.tmp[[4]]
   cMSE1 <- irls.tmp[[5]]
   cmseSpikeLevels <- irls.tmp[[6]]
   nw1 <- irls.tmp[[9]]
   converged <- irls.tmp[[7]]
   m.model.1 <- irls.tmp[[8]]

   if (converged == 0) { return(list(converged)) }

   # compute IRLS for quadratic model
   irls.tmp <- compIRLS(x,DX2,y,VFnP,oldMSEfnPar2,oldB2,cc,k,RNNR,
                       theta,maxIter,grphs,LN,BW)
   B2 <- irls.tmp[[1]]
   #seB2 <- irls.tmp[[2]]
   mse2 <- irls.tmp[[2]]
   #covB2 <- irls.tmp[[4]]
   p2 <- irls.tmp[[3]]
   MSEfnPar2 <- irls.tmp[[4]]
   cMSE2 <- irls.tmp[[5]]
   cmseSpikeLevels <- irls.tmp[[6]]
   nw2 <- irls.tmp[[9]]
   converged <- irls.tmp[[7]]
   m.model.2 <- irls.tmp[[8]]

   if (converged == 0) {return(list(converged)) }

   irls.tmp <- compIRLS(x,DX3,y,VFnP,oldMSEfnPar3,oldB3,cc,k,RNNR,
                       theta,maxIter,grphs,LN,BW)
   B3 <- irls.tmp[[1]]
   #seB3 <- irls.tmp[[2]]
   mse3 <- irls.tmp[[2]]
   #covB3 <- irls.tmp[[4]]
   p3 <- irls.tmp[[3]]
   MSEfnPar3 <- irls.tmp[[4]]
   cMSE3 <- irls.tmp[[5]]
   cmseSpikeLevels <- irls.tmp[[6]]
   nw3 <- irls.tmp[[9]]
   converged <- irls.tmp[[7]]
   m.model.3 <- irls.tmp[[8]]

   if( converged == 0 ){  return(list(converged))  }

   # compute IRLS for quartic model
   irls.tmp <- compIRLS(x,DX4,y,VFnP,oldMSEfnPar4,oldB4,cc,k,RNNR,
                        theta,maxIter,grphs,LN,BW)
   B4 <- irls.tmp[[1]]
   #seB4 <- irls.tmp[[2]]
   mse4 <- irls.tmp[[2]]
   converged <- irls.tmp[[7]]
   m.model.4 <- irls.tmp[[8]]

   if( converged == 0 ) {return(list(converged)) }
   #print("## end compAllIRLS")
   # if convergence for all models return full list
   return(list( converged,B1,mse1,p1,MSEfnPar1,cMSE1,    #6
                B2,mse2,p2,MSEfnPar2,cMSE2,              #11
                B3,mse3,p3,MSEfnPar3,cMSE3,              #16
                B4,mse4,                                 #18
                m.model.1,m.model.2,m.model.3,m.model.4, #22
                nw1,nw2,nw3))                            #25
    }

#' compIRLS
#'
#' @param xx xx
#' @param DX DX
#' @param yy yy
#' @param VFnP VFnP
#' @param MSEFnP MSEFnP
#' @param B B
#' @param cc cc
#' @param k k
#' @param RNNR RNNR
#' @param theta theta
#' @param maxIter max iteration
#' @param grphs graph
#' @param LN LN
#' @param BW BW
#'
#' @export
#'
#
#function [B,seB,mse,covB,p,MSEfnPar,cMSE,cmseSpikeLevels,converged,...
#  IsError,Err] = compIRLS(x,DX,y,VarFnPar,MSEfnPar,B,c,k,...
#  ReqNonNegResponse,theta,maxIter)
compIRLS <- function(xx, DX, yy, VFnP, MSEFnP, B, cc, k, RNNR, theta, maxIter, grphs, LN, BW)
   {
   converged <- 1  # intialize converge to 1 meaning convergence is TRUE
   innerEpsilon <- 1
   outerEpsilon <- 1
   innerNumIter <- 0
   outerNumIter <- 0

   while( outerEpsilon > theta && outerNumIter < maxIter)
      {
      ooB <- B
      while( innerEpsilon > theta && innerNumIter < maxIter )
         {
         oB <- B
         #[B,seB,mse,covB,p,pred]
         #print("####in loop:")
         #print(t(DX))
         #print(oB)
         wls.tmp <- compWLS(xx,DX,yy,VFnP,MSEFnP,oB,cc,RNNR,BW)
         mean.model <- wls.tmp[[1]]
         B <-wls.tmp[[4]]
         B[is.na(B)] <- 0   # remove possible NAs from oldB from bad fit
         #seB <- summary(mean.model )$coefficients[,2]
         Rmse <-  wls.tmp[[2]]
         #covB <- summary(mean.model)$cov
         p <- wls.tmp[[3]]
         pred <- mean.model$fitted
         nw <- wls.tmp[[5]]
         innerEpsilon <- max(abs(oB - B))
         innerNumIter <- innerNumIter + 1
         }
      resids <- yy - pred # compute raw residuals
      #compute MSE fn
      mse.t <- compMSE(xx,resids,cc,k,MSEFnP,BW); #compute MSE fn
      MSEFnP <- mse.t[[1]]
      cMSE <- mse.t[[2]]
      cmseSL <- mse.t[[3]]
      outerEpsilon <- max(abs(ooB - B))
      outerNumIter <- outerNumIter + 1;
      }
   if(outerNumIter >= maxIter ) {converged <- 0 }
   #print("## end compIRLS")
   return(list(B,Rmse,p,MSEFnP,cMSE,cmseSL,converged,
            mean.model,nw))
   }

#' compWLS
#'
#' @param x x
#' @param DX DX
#' @param y y
#' @param VFnPar VFnPar
#' @param MSEFnP MSEFnP
#' @param oldB oldB
#' @param cc cc
#' @param RNNR RNNR
#' @param BW BW
#'
#' @export
#'
#
compWLS <-function (x, DX, y, VFnPar, MSEFnP, oldB, cc, RNNR, BW)
   {

   # get minimum nonzero response
   minNZR <- min(y[y>0])

   # compute WLS
   vf <- VarFn(x,MSEFnP)  # variances
   sf <- sqrt(vf) #standard deviations
   pred <- as.vector(cbind(rep(1,length(DX[,1])),DX) %*% oldB)  # compute predicted values
   rwts <- BisquareWt(x,pred,sf,cc)
   rwts <- BlendWts(BW,rwts)
   Wt <-  rwts/vf
   Wt <- Wt/sum(Wt)
   p <- length(DX[1,])+1
   n <- length(DX[,1])
   nw <- n*(1-sum(Wt*Wt)) + 1

   #write.csv(Wt,"weights",append=TRUE)
   # [B,seB,mse,covB] = lscov(DX,y,Wt);
   wls.tmp <- lm.wfit(cbind(rep(1,length(DX[,1])),DX),y,Wt)
   B <- wls.tmp$coefficients
   B[is.na(B)] <- 0    # remove possible NAs from oldB from bad fit
   #print("### in compWLS lmfit B")
   #print(B)

   #pred <- as.vector(cbind(rep(1,length(x)),DX) %*% B)  #compute predicted values
   resids <- y - wls.tmp$fitted
   RSS <- sum(Wt*resids^2)
   R.mse <- RSS/(nw-p)
   #print(paste("in compWLS B return =",B))
   #print(" end compWLS")
   return(list(wls.tmp,R.mse,p,B,nw))

   }

#' compMSE
#'
#' @param X x
#' @param Resids Resids
#' @param cc cc
#' @param k k
#' @param VFnP VFnP
#' @param BW BW
#'
#' @export
#'
#
#function [MSEfnPar,cMSE,cmseSpikeLevels,IsError,Err] = compMSE(X,resids,c,k)
compMSE <- function(X, Resids, cc, k, VFnP, BW)
   {
   # compute (conditional) mean squared error (MSE) by spiking level
   ## compute (conditional) mean squared error (MSE) by spiking level
   ## [cMSE,cmseSpikeLevels, nLevels, MSEcnt,IsError, Err] ...
   ##     = robCondMSE(X,resids,c,k)
   cond.mse.tmp <- robCondMSE(X,Resids,cc,k,BW)
   cMSE <- cond.mse.tmp[[1]]
   cmseSL <- cond.mse.tmp[[2]]
   nLevels <- cond.mse.tmp[[3]]
   MSEcnt <- cond.mse.tmp[[4]]
   #print("## in compMSE")

   # check for 0 spiking level (method blanks) and remove for replicate
   # variance model calc
   #if( cmseSL[1] == 0)
   #   {
   #   cmseSL <- cmseSL[2:nLevels]
   #   cMSE <- cMSE[2:nLevels]
   #   MSEcnt <- MSEcnt[2:nLevels]
   #   }

   # check for 0 values in vector of conditional MSEs
   if(all(cMSE == 0))
      {
      stop('All conditional MSEs are 0. Likely cause is data entry error')
      }
   if( any(cMSE == 0))
      {
      warning(paste('Warning: One or conditional MSEs is 0. Likely causes are:/n',
        'data entry error or variance damping (due to insufficient significant digits/n',
        'smoothing or thresholding'))
      VarOK <- cMSE > 0
      cmseSL <- cmseSL[VarOK]
      cMSE <- cMSE[VarOK]
      MSEcnt <- MSEcnt[VarOK]
      nmseLevels <- length(cmseSpikeLevels)
      }

   # DoF for conditonal MSE
   cmseDoF <- MSEcnt

   # estimate model for (conditional) mean squared error (MSE) by spiking
   # level

   mse.tmp <- MseMod(cmseSL,cMSE,cmseDoF,VFnP)
   #print("## end compMSE")
   return(list(mse.tmp,cMSE,cmseSL))
   }

#' robCondMSE
#'
#' @param X x
#' @param Resids Resids
#' @param cc cc
#' @param k k
#' @param BW BW
#'
#' @export
#'
#function [robMSE,SpikeLevels, nLevels, MSEcnt,IsError, Err] ...
#    = robCondMSE(X,resids,c,k)
robCondMSE <- function(X, Resids, cc, k, BW)
   {
   # RepVar function
   #
   # computes MSE at various spiking levels,
   # degrees of freedom and weights for modeling
   #
   # ----------------Outputs-------------------------------------
   #
   # robMSE = robust estimates of MSE by spiking level
   # SpikeLevels = vector of spiking levels with 2 or more results
   # nLevels = number of spiking levels
   # MSEcnt = number of residuals per spiking level
   # IsError = Error flag (0-1)
   # Err = Matlab error structure
   #
   # ----------------Inputs-------------------------------------
   #
   # X = vector of spiked concentrations
   # resids = vector of regression residuals from means model
   # c = tuning constant for biweight
   # k = threshold for Huber estimator
   #

   # check for 0 variance at spiking level (method blanks with all 0 results)
   # and remove for replicate variance model calc
   if(X[1]==0)
      {
      NZndx <- (X != 0)
      X <- X[NZndx]
      Resids <- Resids[NZndx]
      BW <- BW[NZndx]
      }

    # get unique values of X & counts
    STable <- cbind(as.numeric(names(table(X))),table(X))
    # first column is spiking level
    SL <- STable[,1]
    # second column is count; DoF is count
    nw <- STable[,2]
    # number of spiking levels
    nLs <- length(SL)
    robmse <- rep(0,nLs)
    # estimate MSE
    for( i in 1:nLs )
      {
      yi <- Resids[abs(X-SL[i])< 1e-12]
      BWi <- BW[abs(X-SL[i])< 1e-12]
      if (w.var(yi,BWi) <= 1e-12 )
         {
         robmse[i] <- sum(yi*BWi/sum(BWi))^2
         } else
         {
         tmp <- roblocvar4(yi,cc,k,1,BWi)
         robmse[i] <- tmp[[2]] + tmp[[1]]^2
         nw[i] <- tmp[[4]]+1
         }
      }
    return(list(robmse, SL, nLs, nw))
    }

#' MseMod
#'
#' @param SL SL
#' @param RMse RMSE
#' @param MDF MDF
#' @param MsePar MsePar
#'
#' @export
#'

MseMod <- function (SL, RMse, MDF, MsePar)
   #
   # computes Gamma nl model for replicate variance at various spiking
   # levels
   #
   # ----------------Inputs-------------------------------------
   #
   # SpikeLevels = vector of spiking levels with 2 or more results
   # RMse = vector of computed replicate variance at each spiking level
   # MDF = vector of degrees of freedom for each element of MSE
   # MsePar
   # ---------------Outputs-------------------------------------
   #
   # MseFnPar - structure for error variance model, type and parameters
   # IsError = Error flag (0-1)
   # Err = Matlab error structure
   #

   # initialize variables
   #
   {
   maxC <- 2
   nL <- length(SL)
   DoF <- sum(MDF)

   MFnP <- MsePar

   # Take initial values for a, b and c from input MsePar
   A <- max(0,MsePar$a)
   B <- MsePar$b
   CC <- min(MsePar$c,maxC)

   newCoef <- OptABC(A,B,CC,SL,RMse,MDF)

   A <- newCoef[1] ; names(A) <- 'constant'
   B <- newCoef[2] ; names(B) <- 'slope'
   CC <- newCoef[3]; names(CC) <-'power'

   if( B <= 0 || CC <= 1e-2 || (B*max(SL)^CC < 0.10*A))
      {
      MFnP$type<-'constant'
      MFnP$a <- mean(RMse)
      MFnP$b <- 0
      MFnP$c <- 0
      MFnP$minVar <- A
      MFnP$minVar <- MFnP$a
      MFnP$DoF <- DoF
      } else
      {
      if ( A < (1e-6)*mean(RMse) )
         {
         MFnP$type <- 'power'
         MFnP$a <- 0
         MFnP$b <- B
         MFnP$c <- CC
         MFnP$minVar <- mean(RMse[1:2])
         MFnP$DoF <- DoF - 2
         } else
         {
         MFnP$type <- 'constant.power'
         MFnP$a <- A
         MFnP$b <- B
         MFnP$c <- CC
         MFnP$minVar <- A
         MFnP$DoF <- DoF - 3
         }
      }
   return(MFnP)
   }

#' feq
#'
#' @param x x
#' @param y y
#' @param tol tolerance
#'
#' @export
#'
feq <- function(x, y, tol)
   {
   return(abs(x-y) <= tol)
}

#' GammaPar
#'
#' @param mu mean
#' @param v variance
#'
#' @export
#'
GammaPars <- function(mu, v)

   # GammaPars Function
   #
   # Computes the shape and scale parameters of a gamma distribution using the
   # method of moments
   #
   # ----------------Inputs-------------------------------------
   #  mu = mean
   #  v = variance
   #
   # ----------------Outputs-------------------------------------
   #  A = shape
   #  B = scale
   #
   {
   B <- v/mu
   A <- mu/B
   return(list(A,B))
}


#' covProb
#'
#' @param xc xc
#' @param LQL LQL
#' @param UQL UQL
#' @param MFP MFP
#' @param VFP VFP
#' @param MSEFP MSEFP
#' @param CPR CPR
#' @param nn nn
#' @param xxbar xxbar
#' @param ssxx ssxx
#' @param RNNR RNNR
#'
#' @export
#'
covProb <- function(xc, LQL, UQL, MFP, VFP, MSEFP, CPR, nn, xxbar, ssxx, RNNR)
   {
   if( RNNR == 1)
      {
      # assign parameter values for Gamma
      lmm <- MeanFn(xc,MFP)
      lv <- VarFn(xc,MSEFP)*(1+1/nn+((xc-xxbar)^2)/ssxx)
      gpar <-GammaPars(lmm,lv)
      A <- gpar[[1]]
      B <- gpar[[2]]
      # compute coverage probability under Gamma
      cvP <- pgamma(xc*UQL,A, scale=B) - pgamma(xc*LQL,A,scale=B)
      } else
      {
      lmm <- MeanFn(xc,MFP)
      lv <-  VarFn(xc,MSEFP)*(1+1/nn+((xc-xxbar)^2)/ssxx)
      lss <- sqrt(lv)
      DoF <- min(VFP$DoF,MSEFP$DoF)
      #compute coverage probability under Normal
      cvP <- pt((xc*UQL-lmm)/lss,DoF)-pt((xc*LQL-lmm)/lss,DoF)
      }
   return(cvP)
   }

#' GammaEstFn
#'
#' @param xs xs
#' @param LQL LQL
#' @param UQL UQL
#' @param MFP MFP
#' @param VFP VFP
#' @param MSEFP MSEFP
#' @param CPR CPR
#' @param nn nn
#' @param xxbar xxbar
#' @param ssxx ssxx
#' @param RNNR RNNR
#'
#' @export
#'
GammaEstFn <-function(xs, LQL, UQL, MFP, VFP, MSEFP, CPR, nn, xxbar, ssxx, RNNR)
   # GammaEstFn function
   #
   # Estimating function for difference between probability content of gamma
   # distribution between two limits and target coverage. Mean and variance
   # are specified by the independent variable x and with specified mean
   # and variance functions.
   #
   # ----------------Outputs-------------------------------------
   #
   # discrep = difference between probability content of interval
   #           and target coverage
   #
   # ----------------Inputs-------------------------------------
   #
   # xs = spiking level
   # MFP =
   # VFP =
   # MSEFP =
   # mV =
   # LQL =
   # UQL =
   # CovPrReq =
   # -----------------------------------------------------------
   {
   # get mean and variance
   #print("######## in GammaEstFn")
   #print(xs)
   #print(MFP)
   lmm <- MeanFn(xs,MFP)
   lv <- VarFn(xs,MSEFP)*(1+1/nn+((xs-xxbar)^2)/ssxx )
   if( RNNR==1)
      {
      # assign parameters for Gamma
      gpars <- GammaPars(lmm,lv)
      A <-gpars[[1]]
      B <- gpars[[2]]
      # compute coverage probability discrepancy
      discrep <- pgamma(xs*UQL,A,scale=B) - pgamma(xs*LQL,A,scale=B)- CPR
     } else
     {
     # assign parameters for Normal
     lss <- sqrt(lv)
     DoF <- min(VFP$DoF,MSEFP$DoF)
     #compute coverage probability discrepancy
     discrep <- pt((xs*UQL-lmm)/lss,DoF)- pt((xs*LQL-lmm)/lss,DoF)-CPR
     }
   return(discrep)
   }


#' srchLCMRL
#'
#' @param SL SL
#' @param LQL LQL
#' @param UQL UQL
#' @param MFP MFP
#' @param VFP VFP
#' @param MSEFP MSEFP
#' @param CPR CPR
#' @param nn nn
#' @param xxbar xxbar
#' @param ssxx ssxx
#' @param RNNR RNNR
#' @param LLLCMRL LLLCMRL
#'
#' @export
#'
srchLCMRL <-function (SL, LQL, UQL, MFP, VFP, MSEFP, CPR, nn, xxbar, ssxx, RNNR, LLLCMRL)
   {
   # Searches for value of x that sets GammaEstFn value to 0. This is the
   # estimated LCMRL
   #
   # ----------------Outputs-------------------------------------
   #
   # LCMRL = Lowest Concentration Minimum Reporting Limit
   # CovProbMat = matrix with spiking levels and estimated coverage
   #   probabilities for specified quality interval
   #
   # ----------------Inputs ------------------
   #
   # SL = spiking levels
   # LowQL = lower quality limit for recovery of spike
   # UpQL = upper quality limit for recovery of spike
   # CPR = desired probability for acceptable spike recovery
   # MFP = parameters of mean function for response
   # VFP = parameters of variance function for response
   #
   # ----------------Start code---------------------------------
   #
   # get endpoints of spiking interval
   #print("############## in Search LCMRL")
   ResultFlag <- 0
   xmx <- max(SL) #upper end of search interval
   if( LLLCMRL ==0)
      {
      xmn <- min(SL) #tentative lower end of search interval

      # ensure that sign of EstFn is negative at lower end of interval
      while(sign(GammaEstFn(xmn,LQL,UQL,MFP,VFP,MSEFP,CPR,nn,xxbar,ssxx,RNNR))>0)
         {
         xmn <- xmn/2
         ResultFlag <- -1
         ResultMessage <- 'LCMRL is below lowest spiking level'
         }
      } else
      {
      xmn <- LLLCMRL
      if( sign(GammaEstFn(xmn,LQL,UQL,MFP,VFP,MSEFP,CPR,nn,xxbar,ssxx,RNNR))>0)
         {
         ResultFlag <- -5
         LCMRL.lst <- list(LLLCMRL,0,0,0)
         }
      }
    xgrid <- seq(xmn,xmx,length.out=100)
    GEFgrid <- rep(0,100)

    # compute signs of GammaEstFn on grid over (xmn,xmx)

    GEFgrid <- GammaEstFn(xgrid,LQL,UQL,MFP,VFP,MSEFP,
                                    CPR,nn,xxbar,ssxx,RNNR)
    # check to see whether LCMRL is above highest spiking level
    if( all(GEFgrid<0) )
      {
      ResultFlag <- -2
      ResultMessage <- 'LCMRL is above highest spiking level'
      LCMRLval <- 0  # explicitly set to 0; test for 0 in srchDL
      # generate points for plot of QC interval coverage
      xp <- xgrid
      yp <- covProb(xp,LQL,UQL,MFP,VFP,MSEFP,CPR,nn,xxbar,ssxx,RNNR)
      CPMat <- c(xp,yp)
      return(list(LCMRLval,CPMat,ResultFlag,ResultMessage))
      }

    # hunt for change of sign in GEFgrrid
    ndxu <-match(1,sign(GEFgrid[1:100]))
    if ( ndxu > 2) { ndxl <- ndxu - 2 } else { ndxl <- ndxu - 1 }
    xl <- xgrid[ndxl]
    xu <- xgrid[ndxu]

    # check again to see whether LCMRL is above highest spiking level
    if ( !all(GEFgrid[ndxu:100] > 0))
      {
      ResultFlag <- -2
      ResultMessage <- 'LCMRL is above highest spiking level'
      LCMRLval <- 0 # explicitly set to 0; test for 0 in srchDL
      # generate points for plot of QC interval coverage
      xp <- xgrid
      yp <- covProb(xp,LQL,UQL,MFP,VFP,MSEFP,CPR,nn,xxbar,ssxx,RNNR)
      CPMat <- c(xp,yp)
      return(list(LCMRLval,CPMat,ResultFlag,ResultMessage))
      }

    #find zero of EStFn; this is the LCMRL
    if( ResultFlag != -5)
      {
      LCMRL.lst <- try(uniroot(GammaEstFn,c(xl,xu),LQL=LQL,UQL=UQL,
                           MFP=MFP,VFP=VFP,MSEFP=MSEFP,CPR=CPR,
                           nn=nn,xxbar=xxbar,ssxx=ssxx,RNNR=RNNR,tol=1e-8))

      if( class(LCMRL.lst)[1]=="try-error")
         {
         ResultFlag <- -3
         ResultMessage <- geterrmessage()
         }
      }
   switch( as.character(ResultFlag),

      '-5' = {ResultMessage <- 'LCMRL below lowest Spiking Level with all non-zero results: set equal to lowest spiking level with all non-zro results'
              CPMat <- 0},
      '-3' = {ResultMessage <- 'Nonconvergence' },
      '-1' = {ResultMessage <-'Lower spiking level needed to bracket the LCMRL'
             # generate points for plot of QC interval coverage probability
             xp <- 10^seq(log10(xmn),log10(max(SL)),length.out=50)
             yp <- covProb(xp,LQL,UQL,MFP,VFP,MSEFP,CPR,nn,xxbar,ssxx,RNNR)
             CPMat <- cbind(xp,yp)},
      '0'  = {ResultFlag <- 1
             ResultMessage <- 'Valid LCMRL'
             # generate points for plot of QC interval coverage probability
             xp <- 10^seq(log10(xmn),log10(xmx),length.out=50)
             yp <- covProb(xp,LQL,UQL,MFP,VFP,MSEFP,CPR,nn,xxbar,ssxx,RNNR)
             CPMat <- cbind(xp,yp)}
       )
   #print("############## end Search LCMRL")
   return(list(LCMRL.lst,CPMat,ResultFlag,ResultMessage))
   }

#' HubauxVos
#'
#' @param MFP MFP
#' @param VFP VFP
#' @param MSEFP MSEFP
#' @param a a
#' @param b b
#' @param lcmrlval LCMRL value
#' @param sL sL
#' @param RNNR RNNR
#' @param ZRF ZRD
#'
#' @export
#'
HubauxVos <- function (MFP, VFP, MSEFP, a, b, lcmrlval, sL, RNNR, ZRF)
   {
   # function HubauxVos
   #
   # determines regression based detection limit based on modification of
   # approach of Hubaux & Vos
   #
   # ----------------Outputs-------------------------------------
   #
   # DL = detection limit
   # ResultFlag = result flag for HV DL
   # ResultMessage = result message for HV DL
   #     -3 - 'Nonconvergence: ' output.message
   #     -2 - 'PROBLEM: DL appears to be above max spiking level'
   #      0 - not set
   #      1 - 'Valid DL'
   #      2  - 'DL calculated >= LCMRL; set DL = LCMRL'
   # IsError = error flag
   # Err = error messsage structure
   #
   # ----------------Inputs-------------------------------------
   #
   # MeanFnPar = parameter structure for mean response function
   # VarFnPar = parameter structure for variance function of response
   # alpha = probability defining critical response (yc) in HV DL
   # beta = probability defining DL in HV DL--Pr(y <= yc | x = DL) = beta
   #
   # ----------------Start code---------------------------------
   #
   #
   # initialize variables
   #
   #MFP <- MeanFnPar
   #VFP <- VarFnPar
   #MSEFP <- MSEFnPar
   # a <-alpha
   # b <- beta
   ResultFlag <-0
   ResultMessage <- 'General Error in DL determination'

   MV <- max(VFP$minVar,MSEFP$minVar)
   stdv <- sqrt(MV)
   # expected meas conc at 0 spiking level
   y0 <- MeanFn(0,MFP)
   # est'd variance at 0 spiking level
   DoF <- min(VFP$DoF,MSEFP$DoF)
   # estimate critical concentration, yc
   if(RNNR == 1)
      {
      if (y0 > 0 )  # assume truncated t if predicted mean is positive
         {
         yc <- y0+stdv*qt(pt(-y0/stdv,DoF)+(1-a)*(1-pt(-y0/stdv,DoF)),DoF)
         } else  # assume half-t distr
         {
         yc <- stdv*qt(1-a/2,DoF)
         }
      } else  # Response can be negative
      {
      yc = y0 + stdv*qt(1-a,DoF) # assume t distribution
      }
   Lc <- yc  # critical level
   srchdl.lst <- srchDL(lcmrlval,yc,sL,MFP,VFP,MSEFP,MV,b,RNNR,ZRF)
   DL <- srchdl.lst[[1]]
   fval <- srchdl.lst[[2]]
   ResultFlag<-srchdl.lst[[3]]
   ResultMessage<-srchdl.lst[[4]]
   return(list(DL,Lc,ResultFlag,ResultMessage))
   }

#' srchDL
#'
#' @param lcmrlval lcmrl value
#' @param yc yc
#' @param sL sL
#' @param MFP MFP
#' @param VFP VFP
#' @param MSEFP MSEFP
#' @param MV MV
#' @param b beta
#' @param RNNR RNNR
#' @param ZRF ZRF
#'
#' @export
#'
srchDL <- function(lcmrlval, yc, sL, MFP, VFP, MSEFP, MV, b, RNNR, ZRF)
   {
   #
   # Searches for value of x that sets DLEstFn value to 0. This is the
   # estimated HVDL
   #
   # ----------------Outputs-------------------------------------
   #
   # DLval = regression based detection limit based on approach of Hubaux & Vos
   # ResultFlag
   # ResultMessage
   # fval
   # exitflag
   # output
   # ----------------Inputs-------------------------------------
   #
   # ----------------Start code---------------------------------
   #
   # initialize variables
   ResultFlag <- 0
   nsL <- length(sL)
   # get endpoints of search interval
   mxspk <- max(sL)
   mnspk <- min(sL)

   # set upper end of search interval
   if( lcmrlval <= 0)   # idicates no valid LCMRL
      {
      xmx <- mxspk   # upper end of search interval
      if(ZRF)
         {
         xmn <- mnspk/10 # tenative lower end of search interval
         } else
         {
         xmn <- mnspk  # tenative lower end of serach interval
         }
      } else
      {
      xmx <- max(lcmrlval,mnspk) # upper end of serach interval
      if(ZRF)
         {
         xmn <- min(lcmrlval,mnspk) # tenative lower end of search interval
         } else
         {
         xmn <- min(lcmrlval,mnspk)/10 # tenative lower end of search interval
         }
      }
   # ensure that sign of EstFn is positive at lower end of interval
   if(ZRF == 0)
      {
      while( sign(DLEstFn(xmn,yc=yc,MFP=MFP,VFP=VFP,MSEFP=MSEFP,MV=MV,
                     b=b,RNNR=RNNR)) < 0 )
         {
         xmn <- xmn/2
         }
      } else
      {
      if( lcmrlval == mnspk || sign(DLEstFn(xmn,yc=yc,MFP=MFP,VFP=VFP,
                  MSEFP=MSEFP,MV=MV,b=b,RNNR=RNNR)) < 0)
         {
         DLval <- min(sL)
         ResultFlag <- -4
         ResultMessage <- 'DL unreliable because of non-zero spiking levels with 0 results'
         return(list(DLval,fval <- NA,ResultFlag,ResultMessage))
         }
      }
   # ensure that sign of EstFn is negative at upper end of interval
   while(sign(DLEstFn(xmx,yc=yc,MFP=MFP,VFP=VFP,MSEFP=MSEFP,MV=MV,
                     b=b,RNNR=RNNR)) == 1 )
         {
         xmx <- 1.2*xmx
         if (xmx > mxspk)
            {
            ResultFlag <- -2
            ResultMessage <- 'PROBLEM: DL may be above max spiking level'
            return(list(DLval<-NA,fval <- NA,ResultFlag,ResultMessage))
            }
         }

   # find zero of EstFn; this is the HV DL
   DL.lst <- try( uniroot(DLEstFn,c(xmn,xmx),yc=yc,MFP=MFP,VFP=VFP,MSEFP=MSEFP,
                          MV=MV,b=b,RNNR=RNNR,lower = xmn, upper = xmx,
                          tol = 1e-6))
    if( inherits(DL.lst, "try-error") )
      {
      ResultFlag <- -3
      ResultMessage <- geterrmessage()
      return(list(DLval<- NA,fval<-NA,Resultflag,Resultmessage))
      }

   DLval <- DL.lst[[1]]
   fval <- DL.lst[[2]]
   # is DL < LCMRL
   if( DLval >= lcmrlval )
      {
      DLval <- lcmrlval
      ResultFlag <- 2
      ResultMessage <- 'DL calculated >= LCMRL; set DL = LCMRL';
      }
   if(ResultFlag == 0) # ResultFlag not set
      {
      ResultFlag <- 1
      ResultMessage <- 'Valid DL'
      }
   return(list(DLval, fval, ResultFlag, ResultMessage))
   }


#' DLEstFn
#'
#' @param xd xd
#' @param yc yc
#' @param MFP MFP
#' @param VFP VFP
#' @param MSEFP MSEFP
#' @param MV MV
#' @param b b
#' @param RNNR RNNR
#'
#' @export
#'
DLEstFn <-function(xd, yc, MFP, VFP, MSEFP, MV, b, RNNR)
   {
   # DLEstFn function
   #
   # Estimating function for difference between probability content of gamma
   # distribution between two limits and target coverage. Mean and variance
   # are specified by the independent variable x and with specified mean
   # and variance functions.
   #
   # ----------------Outputs-------------------------------------
   #
   # discrep = est'd difference between probability content of interval
   # and target coverage
   #
   # ----------------Inputs-------------------------------------
   #
   # xd = hypothetical spiking level
   #
   # ----------------Start code---------------------------------
   #

   discrep <- -Inf

   #   assign gamma parameter values
   lmm <- MeanFn(xd,MFP)
   lv <- VarFn(xd,MSEFP)
   if( RNNR == 1 && lmm > 0 )
      {
      M <- lmm + 10*sqrt(lv)
      if( lv <= lmm*(M-lmm) )
         {
         gpars <- GammaPars(lmm,lv)
         A <- gpars[[1]]
         B <- gpars[[2]]
         # compute coverage probability discrepancy
         discrep <- pgamma(yc,A,scale=B)-b
         } else
         {
         # to avoid problems with gamma prob, use truncated t
         DoF <- MSEFP$DoF
         f <- pt(-lmm/sqrt(lv),DoF)
         discrep <- (pt((yc-lmm)/sqrt(lv),DoF)- f)/(1-f) - b
         }
      } else
      {
      DoF <- MSEFP$DoF
      discrep <- pt((yc-lmm)/sqrt(lv),DoF) - b
      }
   return(discrep)
   }


#' Bayes.1.Boot
#'
#' @param xx xx
#' @param boot.iter boot iteration
#' @param Seed seed
#'
#' @export
#'
  Bayes.1.Boot <- function (xx, boot.iter, Seed)

# computes Bayesian Bootstrap prior weights for Monte Carlo of LCMRL,
# Lc and HV-DL,
#
# ----------------Outputs-------------------------------------
#
# wgt.mat = matrix of prior weights for spiking levels and measurement pairs
#
# ----------------Inputs-------------------------------------
#
# X = vector of spiked concentrations
# boot.iter = number of bootstrap replications
#
# ----------------Start code---------------------------------
#
   {
   stab <- cbind(as.numeric(names(table(xx))),table(xx))
   sL <- stab[,1]     # spiking levels
   nsL <- length(sL)  # number of spiking levels
   nobs <- stab[,2]   # number of observations per spiking level

   W.mat <- matrix(0,boot.iter,length(xx))
   set.seed(Seed)
   W.mat <- apply(W.mat,2,rexp)
   wgt.mat <- W.mat
   W.pos <- 1
   for(i in 1:nsL)
      {
      wgt.mat[,W.pos:(W.pos+nobs[i]-1)]<- t(apply(W.mat[,W.pos:(W.pos+nobs[i]-1)],1,function(z) (length(z)*z/sum(z))))
      W.pos <- W.pos + nobs[i]
      }
   return(wgt.mat)
   }

#' LCMRL.data.rand.pert
#'
#' @param spike spiked
#' @param measured measured
#' @param r.scale.sd r scale sd
#'
#' @export
#'
# LCMRL.data.rand.pert -- function to generate random scale perturbation
# by spiking level of LCMRL data set
#
# model is y*_ij = b_i*y_ij
#
# random scale factor b_i is gamma with mean 1 and std dev r.scale.sd
# independent for each spiking level i

LCMRL.data.rand.pert<-function(spike, measured, r.scale.sd=0.03) {

  spk.lev<-unique(spike)
  n.lev<-length(spk.lev)
  n.reps<-as.vector(table(spike))

  b.s<-r.scale.sd^2
  a.s<-1/b.s
  r.scale<-rgamma(n.lev,shape=a.s,scale=b.s)
  r.sv<-rep(r.scale[1],n.reps[1])
  for(e in 2:n.lev) r.sv<-c(r.sv,rep(r.scale[e],n.reps[e]))

  return(measured*r.sv)

}

#' VarWghts
#'
#' @param xSL xSL
#' @param BW BW
#' @param RW RW
#'
#' @export
#'
   Varwgts <- function(xSL, BW, RW)
   # This function c
   #
   # ----------------Inputs-------------------------------------
   # xSL = column vector of Spiking concentrations
   # BW = Bayesian Bootstrap Weights normalized weights from robust estmation
   # RW = Robust Weights
   #
   # ----------------Outputs-------------------------------------
   #
   # VW = weights to be used for variamce modeling
   {
   stab <- cbind(as.numeric(names(table(xSL))),table(xSL))
   sL <- stab[,1]     # spiking levels
   nsL <- length(sL)  # number of spiking levels
   nobs <- stab[,2]   # number of observations per spiking level

   VW <- BW*RW
   tmp <- aggregate(VW,list(xSL),sum)
   #tot <- sum(tmp[,2])
   dof <- nobs-1
   #return(tmp[,2]*dof/tot)
   return(tmp[,2]*dof)
   }

#' std.z
#'
#' @param x x
#' @param locvar location variable
#'
#' @export
#'
std.z <- function(x, locvar)
   # This function standarizes a list of vector of values x to a list of vectors random variable with 0
   # mean and unit variance
   #  x = a list of numeric vectors containg the LCMRLs of the BB by lab
   #  locvar =  a matrix of lists with columns equal to length(x)
   #        first row being robust means of x from roblocvar4 for lab i
   #        second row being robust variances of x from roblocvar4 for lab i
   {
   n <- length(x)
   std.res <- vector("list",n)
   for( i in 1:n)
      {
      diff <- x[[i]]-as.numeric(locvar[1,])[i]
      std.res[[i]] <- diff/sqrt(as.numeric(locvar[2,])[i])
      }
   return( std.res )
   }

#' std2.z
#'
#' @param x x
#' @param locvar location variable
#'
#' @export
#'
std2.z <- function(x, locvar)
   # This function standarizes a list of vector of values x to a list of vectors random variable with 0
   # mean and unit variance
   #  x = a list of numeric vectors containg the LCMRLs of the BB by lab
   #  locvar =  a matrix of lists with columns equal to length(x)
   #        first row being robust means of x from roblocvar4 for lab i
   #        second row being robust variances of x from roblocvar4 for lab i
   {
   n <- length(x)
   std.res <- vector("numeric",n)
   for( i in 1:n)
      {
      diff <- x[[i]]-locvar[1]
      std.res[i] <- diff/sqrt(locvar[2])
      }
   return( std.res )
   }

#' guttman.utl
#'
#' @param x x
#' @param a alpha
#' @param b beta
#'
#' @export
#'
guttman.utl <- function(x, a=0.95, b=0.95)
   # This function returns the Guttman's non-parametric UTL of a vector x of
   # values for given UTL coverage of b and desired comfidence a. Both a and
   # b are in (0,1).
   # ----------------Inputs-------------------------------------
   #
   # x = numeric vector
   # a = confidence coefficent in (0,1)
   # b = required covergae in (0,1)
   # ---------------Outputs-------------------------------------
   #
   # utl =  upper tolerance limit of desired coverage
   #

   {
   n <- length(x)
   p.est <- floor(b*n)
   k <- p.est:n
   pp <- vector("numeric",length(k))
   pp <- pbeta(1-b,n-k+1,k-1)
   dif <- abs(a-pp)
   kutl <- match(min(dif),dif)
   return(as.vector(sort(x)[k[kutl]]))
   }

#' wgt.guttman.utl
#'
#' @param x x
#' @param y y
#' @param a alpha
#' @param b beta
#'
wgt.guttman.utl <- function(x, y, a=0.95, b=0.95)
   {
   #
   # This function returns a one sided upper Guttman's non-parametric UTL of
   # a vector x of values with associated weights y for given UTL coverage
   # of b and desired comfidence a. Both a and b are in (0,1).
   #
   # ----------------Inputs-------------------------------------
   #
   # x = numeric vector
   # y = a normalized vector of weights summing to 1 associated with each x
   # a = confidence coefficent
   # b = required coverage
   # ---------------Outputs-------------------------------------
   #
   # utl =  upper tolerance limit of desired coverage
   #

   if( (a<0 || a>1) || (b<0 || b>1) )
      {
      stop("Error in wgt.guttma.utl: a and b must be between 0 and 1")
      }
   xs <- sort(x,index=TRUE)
   x <- xs$x
   y <- y[xs$ix]
   csum.y <- cumsum(y)
   lenx <- length(x)
   # find guttman UTL of raw data
   r.utl <- guttman.utl(x,a,b)
   #find position of UTL in raw data
   flag <- x>=r.utl
   p2 <- match(TRUE,flag)
   p1 <- p2-1
   slp <- 1/(lenx+1)
   slp <- slp/(x[p2]-x[p1])
   cnst <- p2/(lenx+1) - slp*x[p2]
   r.utl.p <- slp*r.utl+cnst
   #find the postion r.utl.p in the weighted data
   flag <- csum.y >= r.utl.p
   p2 <- match(TRUE,flag)
   p1 <- p2-1
   slp <- 1/(csum.y[p2]-csum.y[p1])
   cnst <- p2-slp*csum.y[p2]
   w.utl.p <- slp*r.utl.p+cnst
   # find the LCMRL value at the weighted position w.utl.p
   p1 <- floor(w.utl.p)
   p2 <- p1+1
   slp <- (x[p2]-x[p1])
   cnst <- x[p2]-slp*p2
   utl <- slp*w.utl.p+cnst
   return(utl)
   }

#' wt.bw.1.sided
#'
#' @param x x
#' @param rloc rloc
#' @param rscale rscale
#' @param cc cc
#' @param reltol relocation tolerance
#' @param maxsteps max steps
#'
#' @export
#'
#'
wt.bw.1.sided<-function(x, rloc, rscale, cc=6, reltol=1e-6, maxsteps=100)
   {
   # compute weights for one-sided Biweight(c)
   # M-estimate of location

    rlold <- rloc
    delta <- 1
    nsteps <- 0
    while(delta > reltol && nsteps <= maxsteps)
      {
      u <- ifelse(x>rloc,(x - rloc)/(cc*rscale),0)
      tu <- abs(u) <= 1
      rwtsBW <- ((1-u^2)^2) * tu
      rwtsBW <- rwtsBW/sum(rwtsBW)
      rloc <- sum(rwtsBW*x)
      delta <- try(abs((rl_old - rloc)/rl_old),silent=TRUE)
      if( inherits(delta, "try-error") ) {delta <- 0 }
      nsteps <- nsteps +1
      rlold <- rloc
      }
   return(list(loc=rloc,weights=rwtsBW))
   }

