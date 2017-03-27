Fit.sMLSBM<-function(MLNet,MaxK,S){

#Description: This function is the generic function for fitting the sMLSBM model to data
#For Bugs and comments, please contact Natalie Stanley. stanleyn@email.unc.edu
#Last update March 25, 2017

#Inputs:
	#MLNet: A list object where each entry is the network adjacency matrix from a particular layer
	#MaxK: The maximum number of communities you expect in each layer
	#S: the number of strata (clusters of layers)

#Returns:
#A list oject with 3 items
	# $Strata: the vector of layer to strata assignments
	# $Thetas: A list object where the i-th entry gives the theta (probability matrix) from the SBM in the i-th stratum
	# $Comms: A list object where the i-th entry gives a vector of the node-to-community assignments in the i-th layer.

#dependencies: please install before.
library('mixer')
library('phyclust')
library('ttutils')
n=nrow(MLNet[[1]])

####Initialization Phase####

#Do initial layer fits
InitLayerFit<-list()
NNPiList<-list()

	for(i in 1:length(MLNet)){
		MixerRun<-mixer(MLNet[[i]],qmax=MaxK)
		GetBestModel<-getModel(MixerRun)
		InitLayerFit[[i]]<-GetBestModel$Pis
		Pis<-GetBestModel$Pis
		Taus<-GetBestModel$Taus
		InferComm<-apply(Taus,2,function(x) which(x==max(x))[1])
		NNPi<-ConvertThetaToNByN(Pis,InferComm,n)

		#identify isolated vertices
		ZeroInds<-which(rowSums(MLNet[[i]])==0)
		NNPi[ZeroInds,]<-0
		NNPi[,ZeroInds]<-0
		NNPiList[[i]]<-NNPi
	}

	#Cluster these ininital layers
	StrataAssn<-ClSBM(NNPiList,S)
	
	print('done intializing!')

#for monitoring convergence of layer to strata assns
RandInd<-.7

StoreInferComm<-list()

#Store final strata tau and comm memberships
FinalTheta<-list()
FinalComm<-list()

while(RandInd<.99){
UniqueStrata<-sort(unique(StrataAssn))
#This LayerFits list is to hold the nxn represesentation from each strata from the fixed tau
LayerFitsFixTau<-list()
#This list is to hold the nxn representation for each strata from fixed theta
LayerFitsFixTheta<-list()
#use consensus tau to compute a theta for each layer in the strata
 	for(i in 1:length(UniqueStrata)){
 		#Find relevant layers 
 		IndLayer<-which(StrataAssn==UniqueStrata[i])
 		#Pull out relevant layers
 		PullLayer<-MLNet[IndLayer]
 		#determine a k
 		RandNum<-sample(1:length(PullLayer),1)
 		MixerFit<-mixer(PullLayer[[RandNum]],qmin=2,qmax=10)
 		GetModelMixerFit<-getModel(MixerFit)
 		k=nrow(GetModelMixerFit$Pis)
 		FitMLSBM<-SBM.Variational.Fix.Both(PullLayer,n,k,length(PullLayer))
 		Tau<-FitMLSBM$Tau
 		Theta<-FitMLSBM$Theta
 		Alpha<-FitMLSBM$Alpha
 		
 #Infer communities from Tau
InferCommTau<-apply(Tau,1,function(x) which(x==max(x))[1])
ZeroIndsTau<-which(rowSums(Tau)==0)

for(j in 1:length(IndLayer)){ 

ThetaLearn<-ThetaFromFixedTau(PullLayer[[j]],Tau,k)
 
TauLearn<-TauFromFixedTheta(PullLayer[[j]],Theta,Alpha,Tau,k)

ZeroIndsLearn<-which(rowSums(TauLearn)==0)
InferCommLearnTau<-apply(TauLearn,1,function(x) which(x==max(x))[1])
IndSub2<-IndLayer[j]


Mat1<-ConvertThetaToNByN(ThetaLearn,InferCommTau,n)
Mat1[ZeroIndsTau,]<-0
Mat1[,ZeroIndsTau]<-0
LayerFitsFixTau[[IndSub2]]<-Mat1

Mat2<-ConvertThetaToNByN(Theta,InferCommLearnTau,n)
Mat2[ZeroIndsLearn,]<-0
Mat2[,ZeroIndsLearn]<-0
LayerFitsFixTheta[[IndSub2]]<-Mat2

 } #j
 FinalTheta[[i]]<-Theta
 #Zero Out final comms
 Communities<-InferCommTau
 Communities[ZeroIndsTau]<-0
 FinalComm[[i]]<-Communities
} #i

NumCol=n*(n+1)/2
NewPointCloud=matrix(0,nrow=2*length(LayerFitsFixTheta),ncol=NumCol)
for(ii in 1:length(LayerFitsFixTau)){
	TempTau<-LayerFitsFixTau[[ii]]
	TempTheta<-LayerFitsFixTheta[[ii]]
	NewPointCloud[((2*ii)-1),]<-TempTau[lower.tri(TempTau,diag=TRUE)]
	NewPointCloud[(2*ii),]<-TempTheta[lower.tri(TempTheta,diag=TRUE)]
}

#Create a matrix to store the clusters
ClusterStore<-matrix(0,nrow=(nrow(NewPointCloud)/2),ncol=2)

#K means clustering
KMeansCL<-kmeans(NewPointCloud,centers=S)$cluster

ClusterStore<-matrix(KMeansCL,ncol=2,byrow=TRUE)

StoreRelInds<-list()
UniqueRowPatterns<-unique(ClusterStore[,])

NewStrataAssn<-rep(0,nrow(ClusterStore))

for(kk in 1:nrow(UniqueRowPatterns)){
	NewInds<-which(ClusterStore[,1]==UniqueRowPatterns[kk,1] & ClusterStore[,2]==UniqueRowPatterns[kk,2])
	StoreRelInds[[kk]]<-NewInds
	NewStrataAssn[NewInds]<-kk
}	

#Compute rand index between new strata assignments and old assignments
RandInd<-RRand(NewStrataAssn,StrataAssn)$Rand

#Update Strata Assn
StrataAssn<-NewStrataAssn

} #while
Out<-vector(mode='list',length=3)
names(Out)<-c('Strata','Thetas','Comms')
Out[[1]]<-StrataAssn
Out[[2]]<-FinalTheta
Out[[3]]<-FinalComm
Out

}
#####################
##Helper functions###
######################

ConvertThetaToNByN<-function(PiMat,CommAssign,n){
	NNNewMat<-matrix(0,nrow=n,ncol=n)
	for(i in 1:n){
		for(j in 1:n){
			if(i<=j){
				Commi<-CommAssign[i]
				Commj<-CommAssign[j]
				NNNewMat[i,j]<-PiMat[Commi,Commj]
				NNNewMat[j,i]<-NNNewMat[i,j]
			}
		}
	}
	NNNewMat
}

###################################################

ClSBM<-function(FittedMLSBM,NumModel){
	#For dimensions
	Extract1<-FittedMLSBM[[1]]
	ForCol<-length(Extract1[lower.tri(Extract1,diag=TRUE)])
	#Create a data matrix to store results in
	DataMatrixSBM<-matrix(0,nrow=length(FittedMLSBM),ncol=ForCol)
	#Assign layers to rows of the data matrix
	for(i in 1:length(FittedMLSBM)){
		MatPreProcess<-FittedMLSBM[[i]]
		FinalMat<-MatPreProcess[lower.tri(MatPreProcess,diag=TRUE)]
		DataMatrixSBM[i,]<-FinalMat
	}
ResultSBM<-list()
KmeansSBMRes<-kmeans(DataMatrixSBM,centers=NumModel)$cluster
	
KmeansSBMRes

}

#################################################################
SBM.Variational.Fix.Both<-function(Network,N,K,L){
	SuperAdj=MultiLayerInit(Network,N)
    Clustering1=hclust(dist(SuperAdj))
    Clustering<-cutree(Clustering1,k=K)

    Tau<-matrix(0,nrow=N,ncol=K)
    for(i in 1:N){
        Cl<-Clustering[i]
        Tau[i,Cl]<-1
    }

	#Declare Theta#
    Theta<-matrix(0,nrow=K,ncol=K,byrow=TRUE)
    
    TauDiff<-5
    
    #Define tolerance for convergence#
    TolMat<-matrix(0.01,nrow=N,ncol=K)
    Tol<-norm(TolMat,type='F')
    
    IsoVertStore<-list()
    for(iv in 1:length(Network)){
        isonet<-Network[[iv]]
        ZeroInds<-which(rowSums(isonet)==0)
        IsoVertStore[[iv]]<-ZeroInds
    }

    Isovert<-Reduce(intersect,IsoVertStore)
    
    while(TauDiff>Tol){
        if(length(Isovert)>0){
            SubTau<-Tau[rowSums(Tau)==1,]
            Alpha<-colSums(SubTau)/(nrow(SubTau))
        }
        else{
            Alpha<-colSums(Tau)/N
        }
        
        for (i in 1:N) {
            Tau[i,] = Tau[i,] + 0.0000001
            Tau[i,] = Tau[i,]/sum(Tau[i,])
        }

        ThetaNew<-Theta
        for(K1 in 1:K){
        	for(K2 in K1:K){
        		Top<-0
        		Bottom<-0
        		for(l in 1:L){
        		    for(i in 1:N){
        			    Noti<-c(1:N)[-i]
        			    for(j in Noti){
                            Top<-Top+(Tau[i,K1]*Tau[j,K2]*Network[[l]][i,j])
        				    Bottom<-Bottom+(Tau[i,K1]*Tau[j,K2])
        		        } # j
        	        } # i
        		} # l
                if (Top == 0)
                {
                   
                }
        		Bottom<-max(Bottom,0.0000001)
        		ThetaNew[K1,K2]<-min(.99,Top/Bottom)
        		ThetaNew[K2,K1]<-ThetaNew[K1,K2]
        	} # K2
        } # K1
        
        Theta<-ThetaNew
        
        newTau<-matrix(0,nrow=N,ncol=K)
        
        for(ii in 1:N){
            Notii<-c(1:N)[-c(ii,Isovert)]
            for(q in 1:K){
        		Prod<-0
        		for(ll in 1:L){
        			for(jj in Notii){
        				for(q2 in 1:K){
        					Prod<-Prod+((Tau[jj,q2]*Network[[ll]][ii,jj]*log(Theta[q,q2]))+((1-Network[[ll]][ii,jj])*Tau[jj,q2]*log(1-Theta[q,q2])))
        				}
        			}
        		}
        		newTau[ii,q]<-log(Alpha[q])+Prod
        	}
        if(sum(newTau[ii,])==0){
        		newTau[ii,]<-0.00001
        	}

            
            a=max(newTau[ii,])
            Shift<-newTau[ii,]-a
            New<-a+log(sum(exp(Shift)))
            newTau[ii,]<-exp(newTau[ii,]-New)
        	

        } 
        TauDiff<-norm(abs(Tau-newTau),type='F')/(N)
        Tau<-newTau
    } #while
    
    FixBothOut<-vector(mode='list',length=3)
    names(FixBothOut)<-c('Theta','Tau','Alpha')
    if(length(Isovert)>0){
        Tau[Isovert,]<-0
    }
    FixBothOut[[1]]<-Theta
    FixBothOut[[2]]<-Tau
    FixBothOut[[3]]<-Alpha
    FixBothOut
}

#################################################################
ThetaFromFixedTau<-function(Network,Tau,k){
	NewTheta<-matrix(0,nrow=k,ncol=k)
	N<-nrow(Network)
	for(k1 in 1:k){
		for(k2 in 1:k){
			if(k1<=k2){
		Top<-0
    	Bottom<-0
    	for(i in 1:N){
    		Noti<-c(1:N)[-i]
    		for(j in Noti){
    			Top<-Top+(Tau[i,k1]*Tau[j,k2]*Network[i,j])
    			Bottom<-Bottom+(Tau[i,k1]*Tau[j,k2])
    		}
    	}
    	Bottom<-max(Bottom,0.0000001)
        NewTheta[k1,k2]<-min(.99,Top/Bottom)
    	NewTheta[k2,k1]<-NewTheta[k1,k2]   
	}
		}
	}
Inds<-which(NewTheta==0)
NewTheta[Inds]<-.0001
NewTheta
}
#####################################################################
TauFromFixedTheta<-function(Network,Theta,Alpha,Tau,k){
N=nrow(Network)
K=k

Isovert<-which(rowSums(Network)==0)
newTau<-matrix(0,nrow=N,ncol=K)
        
        for(ii in 1:N){
        	Notii<-c(1:N)[-c(ii,Isovert)]
        	for(q in 1:K){
        		Prod<-0
        		for(jj in Notii){
        			for(q2 in 1:K){
        			#Convert this to log domain
        				Prod<-Prod+((Tau[jj,q2]*Network[ii,jj]*log(Theta[q,q2]))+((1-Network[ii,jj])*Tau[jj,q2]*log(1-Theta[q,q2])))
        			}
        		}
        	newTau[ii,q]<-log(Alpha[q])+Prod       
        	}
        	
        	if(sum(newTau[ii,])==0){
        		newTau[ii,]<-0.00001
        	}
				a=max(newTau[ii,])
                Shift<-newTau[ii,]-a
                New<-a+log(sum(exp(Shift)))
                newTau[ii,]<-exp(newTau[ii,]-New)
        	}
        
if(length(Isovert)>0){
newTau[Isovert,]<-0
}
newTau
}

MultiLayerInit<-function(MultiLayerNetwork,n){
    
    
    NumLayer<-length(MultiLayerNetwork)
    #Create super adjacency matrix
    SuperAdjacency=matrix(0,n,n)
    for(l in 1:NumLayer){
        Mat=MultiLayerNetwork[[l]]
        SuperAdjacency=SuperAdjacency+Mat
    }
    SuperAdjacency
}
###############################################################################
