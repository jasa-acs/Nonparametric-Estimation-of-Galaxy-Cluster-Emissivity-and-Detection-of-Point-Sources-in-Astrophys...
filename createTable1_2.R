
#This code is implemented in R

#MSE TABLE
matrixM=function(){
	result=c()

	for(lt in 1:2){
		if(lt==1) typenames=c('blockssym','truef')
		else      typenames=c('f4sym','f4')
		
		for(PS in 0:1){
			if(PS==0) result=rbind(result,c( '\\emph{Without}',rep('',5*2)))
			if(PS==1) result=rbind(result,c( '\\emph{With}',rep('',5*2)))
			
			for(M in 2^(7:9)){
				trowQUT=c()
				trowSA=c()
				
				for(type in typenames){#'f4sym','f4')){
					if(type=='truef') nametype='cosmo1'
					if(type=='blockssym') nametype='blocks1'
					if(type=='f4sym') nametype='cosmoBump'
					if(type=='f4') nametype='cosmoAsym'					

					nk=2^(9-log(M)/log(2))*12

					print(c(PS,M,nk,nametype))

					a=readMat(paste('data_and_results/full_astrosim_diag_powerlaw_M',M,type,'PS',PS,'B1W1S1Erealnk',nk,'forcepos2.mat',sep=''))	

					logmseQUT=median(as.numeric(a$summary[,,1]$alllogmsefull))
					logmseQUT=format(logmseQUT*100,digits=2,scientific=FALSE)
					
					logmseQUTc=median(as.numeric(a$summary[,,1]$alllogmsec))
					logmseQUTc=format(logmseQUTc*100,digits=2,scientific=FALSE)
					
					logmseQUTm=median(as.numeric(a$summary[,,1]$alllogmsem))
					logmseQUTm=format(logmseQUTm*100,digits=2,scientific=FALSE)
					
					logmseQUTb=median(as.numeric(a$summary[,,1]$alllogmseb))
					logmseQUTb=format(logmseQUTb*100,digits=2,scientific=FALSE)
					
	 
					logmseSA=median(as.numeric(a$summary[,,1]$alllogmsedomfull))
					logmseSA=format(logmseSA*100,digits=2,scientific=FALSE)
					
					logmseSAc=median(as.numeric(a$summary[,,1]$alllogmsedomc))
					logmseSAc=format(logmseSAc*100,digits=2,scientific=FALSE)
					
					logmseSAm=median(as.numeric(a$summary[,,1]$alllogmsedomm))
					logmseSAm=format(logmseSAm*100,digits=2,scientific=FALSE)
					
					logmseSAb=median(as.numeric(a$summary[,,1]$alllogmsedomb))
					logmseSAb=format(logmseSAb*100,digits=2,scientific=FALSE)
					
					trowQUT=c(trowQUT,'',logmseQUT,logmseQUTc,logmseQUTm,logmseQUTb)
					trowSA=c(trowSA,'',paste('\\emph{',logmseSA,'}',sep=''),paste('\\emph{',logmseSAc,'}',sep=''),paste('\\emph{',logmseSAm,'}',sep=''),paste('\\emph{',logmseSAb,'}',sep=''))
					
					
					
				}
				
				trowQUT=c(M,'ALIAS',trowQUT[-1])
				trowSA=c('','\\emph{Conventional}',trowSA[-1])

				print(trowQUT)
				print(trowSA)

				result=rbind(result,trowQUT,trowSA)
			}

			rbind(result,c(rep(NA,5*2+1)))
		}
	}

	result=xtable(result,sanitize.colnames.function = identity)
	addtorow=list()
	addtorow$pos=list(0,0,0,0,0,3,5,7,10,12,14,14,14,14,14,14,17,19,21,24,26)
	addtorow$command=c(paste('&& \\multicolumn{9}{c}{MSE of the log-profile (*100)} \\\\\n'),
							'\\cline{3-11}\n',
							'&& \\multicolumn{4}{c}{\\tt cosmoBlocks} && \\multicolumn{4}{c}{\\tt cosmo1} \\\\\n',
							'\\cline{3-6} \\cline{8-11}\n',
							'M&&F&i&m&o&&F&i&m&o\\\\\n',
							'\\\\\n',
							'\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\\\\n',
							'\\\\\n',
							'\\cline{3-11}\n',
							'&& \\multicolumn{4}{c}{\\tt cosmo2} && \\multicolumn{4}{c}{\\tt cosmoAsym} \\\\\n',
							'\\cline{3-6} \\cline{8-11}\n',
							'M&&F&i&m&o&&F&i&m&o\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\\\\n'
							)
	align(result) <- c('l','l','l',rep("c", 9))
	xt=capture.output(print(result,add.to.row=addtorow,include.colnames=FALSE,include.rownames=FALSE,sanitize.colnames.function = identity))
	xt= gsub('$\\backslash$','\\',xt,fixed=TRUE)
	xt= gsub('\\{','{',xt,fixed=TRUE)
	xt= gsub('\\}','}',xt,fixed=TRUE)
	cat(xt, file = "data_and_results/table1.tex", sep = "\n")
	
}

#BIAS TABLE
matrixBIAS=function(){
	result=c()

	for(lt in 1:2){
		if(lt==1) typenames=c('blockssym','truef')
		else      typenames=c('f4sym','f4')
		
		for(PS in 0:1){
			if(PS==0) result=rbind(result,c( '\\emph{Without}',rep('',5*2)))
			if(PS==1) result=rbind(result,c( '\\emph{With}',rep('',5*2)))
			
			for(M in 2^(7:9)){
				trowQUT=c()
				trowSA=c()
				
				for(type in typenames){#'f4sym','f4')){
					if(type=='truef') nametype='cosmo1'
					if(type=='blockssym') nametype='blocks1'
					if(type=='f4sym') nametype='cosmoBump'
					if(type=='f4') nametype='cosmoAsym'					

					nk=2^(9-log(M)/log(2))*12

					print(c(PS,M,nk,nametype))

					a=readMat(paste('data_and_results/full_astrosim_diag_powerlaw_M',M,type,'PS',PS,'B1W1S1Erealnk',nk,'forcepos2.mat',sep=''))	

					logbiasQUT=median(as.numeric(a$summary[,,1]$alllogbiasfull))
					logbiasQUT=format(logbiasQUT*100,digits=2,scientific=FALSE)
					
					logbiasQUTc=median(as.numeric(a$summary[,,1]$alllogbiasc))
					logbiasQUTc=format(logbiasQUTc*100,digits=2,scientific=FALSE)
					
					logbiasQUTm=median(as.numeric(a$summary[,,1]$alllogbiasm))
					logbiasQUTm=format(logbiasQUTm*100,digits=2,scientific=FALSE)
					
					logbiasQUTb=median(as.numeric(a$summary[,,1]$alllogbiasb))
					logbiasQUTb=format(logbiasQUTb*100,digits=2,scientific=FALSE)
					
	 
					logbiasSA=median(as.numeric(a$summary[,,1]$alllogbiasdomfull))
					logbiasSA=format(logbiasSA*100,digits=2,scientific=FALSE)
					
					logbiasSAc=median(as.numeric(a$summary[,,1]$alllogbiasdomc))
					logbiasSAc=format(logbiasSAc*100,digits=2,scientific=FALSE)
					
					logbiasSAm=median(as.numeric(a$summary[,,1]$alllogbiasdomm))
					logbiasSAm=format(logbiasSAm*100,digits=2,scientific=FALSE)
					
					logbiasSAb=median(as.numeric(a$summary[,,1]$alllogbiasdomb))
					logbiasSAb=format(logbiasSAb*100,digits=2,scientific=FALSE)
					
					trowQUT=c(trowQUT,'',logbiasQUT,logbiasQUTc,logbiasQUTm,logbiasQUTb)
					trowSA=c(trowSA,'',paste('\\emph{',logbiasSA,'}',sep=''),paste('\\emph{',logbiasSAc,'}',sep=''),paste('\\emph{',logbiasSAm,'}',sep=''),paste('\\emph{',logbiasSAb,'}',sep=''))
					
					
					
				}
				
				trowQUT=c(M,'ALIAS',trowQUT[-1])
				trowSA=c('','\\emph{Conventional}',trowSA[-1])

				print(trowQUT)
				print(trowSA)

				result=rbind(result,trowQUT,trowSA)
			}

			rbind(result,c(rep(NA,5*2+1)))
		}
	}

	result=xtable(result,sanitize.colnames.function = identity)
	addtorow=list()
	addtorow$pos=list(0,0,0,0,0,3,5,7,10,12,14,14,14,14,14,14,17,19,21,24,26)
	addtorow$command=c(paste('&& \\multicolumn{9}{c}{Bias of the log-profile (*100)} \\\\\n'),
							'\\cline{3-11}\n',
							'&& \\multicolumn{4}{c}{\\tt cosmoBlocks} && \\multicolumn{4}{c}{\\tt cosmo1} \\\\\n',
							'\\cline{3-6} \\cline{8-11}\n',
							'M&&F&i&m&o&&F&i&m&o\\\\\n',
							'\\\\\n',
							'\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\\\\n',
							'\\\\\n',
							'\\cline{3-11}\n',
							'&& \\multicolumn{4}{c}{\\tt cosmo2} && \\multicolumn{4}{c}{\\tt cosmoAsym} \\\\\n',
							'\\cline{3-6} \\cline{8-11}\n',
							'M&&F&i&m&o&&F&i&m&o\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\\\\n'
							)
	align(result) <- c('l','l','l',rep("c", 9))
	xt=capture.output(print(result,add.to.row=addtorow,include.colnames=FALSE,include.rownames=FALSE,sanitize.colnames.function = identity))
	xt= gsub('$\\backslash$','\\',xt,fixed=TRUE)
	xt= gsub('\\{','{',xt,fixed=TRUE)
	xt= gsub('\\}','}',xt,fixed=TRUE)
	cat(xt, file = "data_and_results/tablebias.tex", sep = "\n")
	
}


#COVERAGE TABLE
matrixMcoverage=function(){
	result=c()

	for(lt in 1:2){
		if(lt==1) typenames=c('blockssym','truef')
		else      typenames=c('f4sym','f4')
		
		for(PS in 0:1){
			if(PS==0) result=rbind(result,c( '\\emph{Without}',rep('',5*2)))
			if(PS==1) result=rbind(result,c( '\\emph{With}',rep('',5*2)))
			
			for(M in 2^(8:9)){
				#trow=list(c(),c(),c())
			  trow=c()
			  trowdom=c()
				
				for(type in typenames){#'f4sym','f4')){
					if(type=='truef') nametype='cosmo1'
					if(type=='blockssym') nametype='blocks1'
					if(type=='f4sym') nametype='cosmoBump'
					if(type=='f4') nametype='cosmoAsym'					

					nboot=24*(3-(log(M)/log(2)-7))

					print(c(PS,M,nboot,nametype))
					
					#for(ktype in 0:2){
						
						#a=readMat(paste('data_and_results/rsummary_diag_powerlaw_',M,type,'ps',PS,'type',ktype,'nboot',nboot,'.mat',sep=''))
					  a=readMat(paste('data_and_results/rsummary_diag_powerlaw_',M,type,'ps',PS,'type',1,'nboot',nboot,'.mat',sep=''))
            
						if(type=='f4'){
  						covfull=median(as.numeric(a$summary[,,1]$coveragefull))
  						covfull=floor(covfull)
  						
  						covc=median(as.numeric(a$summary[,,1]$coveragefullc))
  						covc=floor(covc)
  						
  						covm=median(as.numeric(a$summary[,,1]$coveragefullm))
  						covm=floor(covm)
  						
  						covb=median(as.numeric(a$summary[,,1]$coveragefullb))
  						covb=floor(covb)
  						
  						covdomfull=median(as.numeric(a$summary[,,1]$coveragedomfull))
  						covdomfull=floor(covdomfull)
  						
  						covdomc=median(as.numeric(a$summary[,,1]$coveragedomfullc))
  						covdomc=floor(covdomc)
  						
  						covdomm=median(as.numeric(a$summary[,,1]$coveragedomfullm))
  						covdomm=floor(covdomm)
  						
  						covdomb=median(as.numeric(a$summary[,,1]$coveragedomfullb))
  						covdomb=floor(covdomb)
						}
  				  else{
  				    covfull=median(as.numeric(a$summary[,,1]$coverage))
  				    covfull=floor(covfull)
  				    
  				    covc=median(as.numeric(a$summary[,,1]$coveragec))
  				    covc=floor(covc)
  				    
  				    covm=median(as.numeric(a$summary[,,1]$coveragem))
  				    covm=floor(covm)
  				    
  				    covb=median(as.numeric(a$summary[,,1]$coverageb))
  				    covb=floor(covb) 
  				    
  				    covdomfull=median(as.numeric(a$summary[,,1]$coveragedom))
  				    covdomfull=floor(covdomfull)
  				    
  				    covdomc=median(as.numeric(a$summary[,,1]$coveragedomc))
  				    covdomc=floor(covdomc)
  				    
  				    covdomm=median(as.numeric(a$summary[,,1]$coveragedomm))
  				    covdomm=floor(covdomm)
  				    
  				    covdomb=median(as.numeric(a$summary[,,1]$coveragedomb))
  				    covdomb=floor(covdomb)
  				  }
						
						#trow[[ktype+1]]=c(trow[[ktype+1]],'',covfull,covc,covm,covb)
						trow=c(trow,'',covfull,covc,covm,covb)
						trowdom=c(trowdom,'',paste('\\emph{',covdomfull,'}'),paste('\\emph{',covdomc,'}'),paste('\\emph{',covdomm,'}'),paste('\\emph{',covdomb,'}'))
					#}
				}
				
				#trow[[1]]=c(M,'MIN',trow[[1]][-1])
				#trow[[2]]=c('','MEDIAN',trow[[2]][-1])
				#trow[[3]]=c('','MAX',trow[[3]][-1])
				trow=c(M,'ALIAS',trow[-1])
				trowdom=c('','\\emph{Conventional}',trowdom[-1])

				print(trow)
				print(trowdom)

				#result=rbind(result,trow[[1]],trow[[2]],trow[[3]])
				result=rbind(result,trow,trowdom)
			}

			rbind(result,c(rep(NA,5*2+1)))
		}
	}

	result=xtable(result,sanitize.colnames.function = identity)
	addtorow=list()
	#addtorow$pos=list(0,0,0,0,0,4,7,11,14,14,14,14,14,14,18,21,25)
	addtorow$pos=list(0,0,0,0,0,3,5,8,10,10,10,10,10,10,13,15,18)
	
	addtorow$command=c(paste('&& \\multicolumn{9}{c}{Coverage percentage} \\\\\n'),
							'\\cline{3-11}\n',
							'&& \\multicolumn{4}{c}{\\tt cosmoBlocks} && \\multicolumn{4}{c}{\\tt cosmo1} \\\\\n',
							'\\cline{3-6} \\cline{8-11}\n',
							'M&&F&i&m&o&&F&i&m&o\\\\\n',
							'\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\\\\n',
							'\\cline{3-11}\n',
							'&& \\multicolumn{4}{c}{\\tt cosmo2} && \\multicolumn{4}{c}{\\tt cosmoAsym} \\\\\n',
							'\\cline{3-6} \\cline{8-11}\n',
							'M&&F&i&m&o&&F&i&m&o\\\\\n',
							'\\hline\n',
							'\\\\\n',
							'\\hline\n',
							'\\\\\n'
							)
	align(result) <- c('l','l','l',rep("c", 9))
	xt=capture.output(print(result,add.to.row=addtorow,include.colnames=FALSE,include.rownames=FALSE,sanitize.colnames.function = identity))
	xt= gsub('$\\backslash$','\\',xt,fixed=TRUE)
	xt= gsub('\\{','{',xt,fixed=TRUE)
	xt= gsub('\\}','}',xt,fixed=TRUE)
	cat(xt, file = "data_and_results/table2.tex", sep = "\n")
	
}



library(xtable)
library(R.matlab)

#table 1
a=matrixM()

#table 2
b=matrixMcoverage()

#table bias (not in paper)
c=matrixBIAS()



