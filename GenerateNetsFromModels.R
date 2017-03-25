GenerateNetsFromModels<-function(PinVec,n,k,c,NumPerModel){
	#This function creates a list of multilayer networks 
	
	library('ttutils')
	library('mixer')
	library('phyclust')
	library('ttutils')
	NumPerComm<-n/k
	
	CommMembers=CommAssn(NumPerComm,n,k)

	MultiLayerNetwork<-list()

		for(p in 1:nrow(PinVec)){
		
		ParameterRes=ComputeParameters(c,PinVec[p],k,NumPerComm)
		
		pout<-ParameterRes[3]
		#Create layers from the model
			for(l in 1:NumPerModel){
				NewNetList<-list()
				NewNet<-PinPoutNetNew(PinVec[p],pout,n,k,CommMembers,c)
				NewNetList[[1]]<-NewNet
				#Merge with existing list
				MultiLayerNetwork<-merge(MultiLayerNetwork,NewNetList)
			}
		}
		MultiLayerNetwork
}

####################
##Helper Functions##
####################

CommAssn<-function(NumPerComm,n,k){
	#A function to generate a vector of community assignments
	#Inputs:
		#NumPerComm: Number of nodes per community
		#n: number of nodes
		#k: the number of communities. 

		TotalCommVec<-c()

		for(i in 1:k){
			TotalCommVec<-c(TotalCommVec,rep(i,NumPerComm))
		}
		TotalCommVec
}

ComputeParameters<-function(c,pin,k,NumPerComm){
	#A function to compute correct cout and pout for a fixed c (mean degree) and a choice of pin
	#Inputs:
		#c: Mean degree
		#pin: Within group probability of an edge
		#k: the number of communities
		#NumPerComm: Number of vertices per community. So, N=NumPerComm*k

	#Test to make sure that cin is smaller than c

	#Compute cin (mean in degree)
	cin=NumPerComm*pin

	#Compute cout (mean out degree)
	cout=c-cin

	#Compute pout
	pout=cout/((k-1)*NumPerComm)

	#Put these in a vector
	ParameterVec=c(cin,cout,pout,pin)
	names(ParameterVec)<-c('cin','cout','pout','pin')
	
	ParameterVec
}

PinPoutNetNew<-function(pin,pout,n,k,Comms,c){
	
	NumPerComm=n/k

	Theta=matrix(pout,nrow=k,ncol=k)
	diag(Theta)<-pin
	
	#Initialize Adjacency Matrix
	Adj<-matrix(0,nrow=n,ncol=n)
	#Fill in based on bernoullis
	for(i in 1:n){
		for(j in 1:n){
			if(i<j){
				Commi=Comms[i]
				Commj=Comms[j]
				Param=Theta[Commi,Commj]

				Adj[i,j]<-rbinom(1,1,Param)
				Adj[j,i]<-Adj[i,j]
			}
			
		}
	}
	Adj
}
