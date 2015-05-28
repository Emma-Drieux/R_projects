#setwd("R")
#source("tme3")
###TME sur les HMM

##données initiales : matrice des états cachés
x=c(0.99,0.05,0.01,0.95)
p=matrix(x,nrow=2,ncol=2)

## : matrice des probabilités d'émission des symboles
y=c(0.5,0.9,0.5,0.1)
e=matrix(y,nrow=2,ncol=2)

##condition initiales
pi0=c(0.99,0.01)

##fonction qui simule T jets de pièces
simul<-function(T){
	A=matrix(0,T,2)
	A[1,1]="F"
	A[1,2]=sample(c("H","T"),1,prob=c(0.5,0.5))[1]

	for( i in 2:T ){ 
		if(A[i-1,1]=="F") 	A[i,1]=sample(c("F","U"),1,prob=c(0.99,0.01))[1]
		else 			A[i,1]=sample(c("F","U"),1,prob=c(0.05,0.95))
		if(A[i,1]=="F") 	A[i,2]=sample(c("H","T"),1,prob=c(0.5,0.5))
		else 			A[i,2]=sample(c("H","T"),1,prob=c(0.9,0.1))

	}
	return (A)
}

seq=simul(2000)

##algorithme de Viterbi
viterbi<-function(p,e,pi0,a){
	n=length(a)
	A=matrix(0,n,2)
	pl=log(p)
	el=log(e)
	pil=log(pi0)
	id=matrix(0,n,1)
	sol=list()

	if(a[1]=="H") 	A[1,]=pil+t(el[,1])
	else 		A[1,]=pil+t(el[,2])
	for( i in 2:n ){ 
		for( j in 1:2 ){ 
			if(a[i]=="H"){
 				A[i,j]=el[j,1]+max(A[i-1,1]+pl[1,j],A[i-1,2]+pl[2,j])
			}

			else 	A[i,j]=t(el[j,2])+max(A[i-1,1]+pl[1,j],A[i-1,2]+pl[2,j])
		}
	}
	id[n]=which.max(A[n,])
	for( k in (n-1):1 ){
		temp=A[k,]+t(pl[,id[k+1]])
		id[k]=which.max(temp)
	}
	for( k in 1:n ){
		if(id[k]==1)	id[k]="F"
		else 		id[k]="U"
	}
	rep=matrix(0,n,2)
	rep[,1]=seq [,1]
	rep[,2]=id
	m=(rep[,1]==rep[,2])
	puissance=length(m[m=="TRUE"])
	sol$rep=rep
	sol$puissance=puissance
	return (sol)
}
res=viterbi(p,e,pi0,t(seq[,2]))

##fonction qui compte les nombres d'occurence
NM<-function(i,a){
	n=length(i)
	N=matrix(1,2,2)
	M=matrix(1,2,2)
	sol=list()

	for( k3 in 1:(n-1) ){
		if(i[k3]=="F"){
			if(i[k3+1]=="F")	N[1,1]=N[1,1]+1
			else 			N[1,2]=N[1,2]+1
			if(a[k3]=="H")		M[1,1]=M[1,1]+1
			else 			M[1,2]=M[1,2]+1
		}
		else if(i[k3]=="U"){
			if(i[k3+1]=="F") 	N[2,1]=N[2,1]+1
			else 			N[2,2]=N[2,2]+1
			if(a[k3]=="H")  	M[2,1]=M[2,1]+1
			else 			M[2,2]=M[2,2]+1
		}				

	}
	if(i[n]=="F"){
		if(a[n]=="H")	M[1,1]=M[1,1]+1
		else 		M[1,2]=M[1,2]+1
	}
	else if(i[n]=="U"){
		if(a[n]=="H")  	M[2,1]=M[2,1]+1
		else 		M[2,2]=M[2,2]+1
	}
	sol$N=N
	sol$M=M
	return(sol)
}

nm=NM(seq[,1],seq[,2])

##fonction qui estime e et p à partir des nombres d'ocurrence
pe<-function(nm){
	sol=list()
	e2=matrix(0,2,2)
	p2=matrix(0,2,2)

	for( i in 1:2 ){
		for( j in 1:2 ){
			e2[i,j]=nm$M[i,j]/rowSums(nm$M)[i]
			p2[i,j]=nm$N[i,j]/rowSums(nm$N)[i]
		}
	}

	sol$p2=p2
	sol$e2=e2
	return(sol)
}

pe1=pe(nm)


##fonction qui calcule le log odd ratio de deux couples s=(p,e)
logOddRatio<-function(a,ii,s1,s2){
	n=length(a)
	ps=0

	for( k in 1:n ){
		if(a[k]=="H") 	j=1
		else  		j=2
		if(ii[k]=="F") 	i=1
		else 		i=2
		if(k>1 && ii[k-1]=="F") 	i0=1
		else 				i0=2
		if(k>1)	ps=ps+log(s1$p2[i0,i]/s2$p2[i0,i])
		ps=ps+log(s1$e2[i,j]/s2$e2[i,j])	

	}

	return (ps)
}

##algorithme Viterbi training
vt<-function(a,pi0){

	max=100
	sol2=list()
	n1=runif(2)
	e2=matrix(c(n1,1-n1),2,2)
	n2=runif(2)
	p2=matrix(c(n2,1-n2),2,2)
	k=1

	while( k<max){
		res=viterbi(p2,e2,pi0,a)
		nm=NM(res$rep[,2],t(a))
		sol=pe(nm)

		if(length(sol2)==0){
			sol2$e2=sol$e2
			sol2$p2=sol$p2
		}
		else if(length(sol2)>0 && logOddRatio(a,res$rep[,1],sol,sol2)>0){ 
			sol2$e2=sol$e2
			sol2$p2=sol$p2
		}

		n1=runif(2)
		e2=matrix(c(n1,1-n1),2,2)
		n2=runif(2)
		p2=matrix(c(n2,1-n2),2,2)
		k=k+1
	}
	return(sol2)
}

vt1=vt(t(seq[,2]),pi0)













