#Question 1-3.
N=30
tfinal=1000000		#instants t=1,...,10exp(6)
T=1.5
s=sample(c(-1,+1),N,replace=TRUE)	#genere de facon aleatoire N valeurs egales a -1 ou +1;;;; donnees de depart
M=array(1:tfinal)	#genere un vecteur liste allant de 1 a tfinal
for (t in 1:tfinal){
	i=sample(1:N,1)	#selectionner une position aleatoire
	pplus= exp(1/(N*T)*sum(s[-i])) /(exp(-1/(N*T)*sum(s[-i]))+exp(1/(N*T)*sum(s[-i])))
	y=sample(c(-1,+1),1,prob=(c(1-pplus,pplus)))
	s[i]=y
	M[t]= mean(s)
	
	}
#print(M)
#plot(M)
#data <- M # création d'un objet "data" contenant les valeurs 
#jpeg("mongraph.jpg")
#plot(data)
#dev.off() # permet d'ouvrir le fichier .jpg sans fermer R

#Question 4.
hist(M,10)
data <- M # création d'un objet "data" contenant les valeurs
jpeg("mongraph.jpg")
hist(data)
dev.off() # permet d'ouvrir le fichier .jpg sans fermer R
