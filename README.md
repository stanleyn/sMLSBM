# sMLSBM
This code is for fitting the strata multilayer stochastic block model described in 

Clustering network layers with the strata multilayer stochastic block model
N Stanley, S Shai, D Taylor, PJ Mucha - IEEE Transactions on Network Science and engineering,2016.
http://ieeexplore.ieee.org/abstract/document/7442167/.

For bugs and comments please contact Natalie Stanley stanleyn@email.unc.edu.

---------------

Background: For a multilayer network with L layers and N nodes within each layer, the objective of the strata multilayer stochastic block model, sMLSBM is to partition the layers into strata (i.e. clusters), where members of the straum we generated from the same stochastic block model. The output of sMLSBM is a partition of layers into strata and the stochastic block model parameters (within/between probability matrices and community assignment vectors describing each stratum).

This code is implemnted in R and was tested in version 3.2.2

---------------

—————————————— 

BEFORE USING:

You will need to install the following packages to use this code.

>install.packages(‘mixer’)
>install.packages(‘ttutils’)
>install.packages(‘phyclust’)

________________

TO FIT sMLSBM:

Please use the function Fit.sMLSBM.R. We assume you are fitting this model to an unweighted undirected network with L layers and each layer has N nodes. 

The inputs to this function are:

	-MLNet: This is an L-length list where the i-th entry gives the NxN adjacency matrix for layer i. 

	The objective in running the code is partitioning these adjacency 	
	matrices into strata (clusters). 

	-MaxK: The maximum number of communities that you expect in any of the layers. 

	For example, MaxK=6 means that none of the layers have more than 6 communities.

	-S: The number of strata (i.e. clusters of layers to find). 

	For example S=3 means that
	we expect the layers have been generated from 3 different stochastic block models. 

The outputs of the function are a list object with 3 entries.

	# $Strata: the vector of layer to strata assignments
	# $Thetas: A list object where the i-th entry gives the theta (probability matrix) from the SBM in the i-th stratum
	# $Comms: A list object where the i-th entry gives a vector of the node-to-community assignments in the i-th layer. 


_______________________

Example:

We will use the function GenerateNetsFromModels.R to create a synthetic example with 3 strata. Each stratum is generated from one of 3 stochastic block models.  

Step 1: Generate the synthetic multilayer network. You can see that each entry of MLNet is a binary adjacency matrix.
>source(‘GenerateNetsFromModels.R’)
>MLNet=GenerateNetsFromModels(matrix(c(0.6,0.3,0.76),ncol=1),100,4,20,10)

Step2: Fit sMLSBM. In the inputs, we are specifying max 4 communities and 3 strata.
>source(‘Fit.sMLSBM.R’)
>Fit=Fit.sMLSBM(MLNet,4,3) 

Step 3: Check out the output

In the synthetic example we had 30 layers.

Fit$Strata is a length 30 vector with layer to strata assignments 

Fit$Thetas is a length-3 list giving the theta parameter for each stratum

Fit$Comms is a length 3 list giving node to community assignments for each stratum. 

