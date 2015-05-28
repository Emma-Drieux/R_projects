#Question-1-

c1=c("H","T")
c2=c("F","U")

trans=matrix(0,nrow=2001,ncol=2)
trans[1,1]=sample(c2,1,prob=c(0.99,0.01)) # condition initale

for (i in 1:2000){
	if (trans[i,1]=="F"){
		trans[i+1,1]=sample(c2,1,prob=c(0.99,0.01))
		trans[i,2]=sample(c1,1,prob=c(0.5,0.5))}
	if (trans[i,1]=="U"){
		trans[i+1,1]=sample(c2,1,prob=c(0.05,0.95))
		trans[i,2]=sample(c1,1,prob=c(0.9,0.1))}
}
# dernière emission
if (trans[2001,1]=="F"){
	trans[2001,2]=sample(c1,1,prob=c(0.5,0.5))}
if (trans[2001,1]=="U"){
	trans[2001,2]=sample(c1,1,prob=c(0.9,0.1))}

#print (trans)

#Question-2-

emissions=t(trans[,2])
pi=log(c(0.99,0.01))
x=c(0.5,0.9,0.5,0.1)
e=matrix(data=log(x),nrow=2,ncol=2)

###
#			    F	U
#			    0.5	0.9   H
#la matrice e a la forme e=
#			    0.5 0.1   T



MatrixDel=matrix(0,nrow=2001,ncol=2)	#matrice des delta de chaque position (delta est un vecteur de forme delta(F,U))
#initialisation	XXXXXXXXX
if (emissions[1]=="H"){
	delta=pi+e[,1]}		#delta est un vecteur
if (emissions[1]=="T"){
	delta=pi+e[,2]}
#XXXXXXXXXXXXXXXXXXXXXXXX
MatrixDel[1,1]=delta[1]
MatrixDel[1,2]=delta[2]
#RecursionXXXXXXXXXXXXXXX
P=log(c(0.99,0.05,0.01,0.95))
for (i in 1:2001){
	Addi=P+matrix((delta),2,2)	#remplit par colonne---en se servant deux fois du vecteur delta
	Maximums=c(max(Addi[,1]),max(Addi[,2]))
	if (emissions[i]=="H"){
		delta=t(e[,1])+Maximums
	}
	if (emissions[i]=="T"){
		delta=t(e[,2])+Maximums
	}
	MatrixDel[i,1]=delta[1]
	MatrixDel[i,2]=delta[2]
}
#XXXXXXXXXXXXXXXXXXXXXXXX
#print (MatrixDel)

#BacktrakingXXXXXXXXXXXXX
transB=matrix(NA,nrow=2001,ncol=2)
for(i in 1:2001){
	transB[i,2]=trans[i,2]}
	#initialisation: valeur pour T
if (MatrixDel[2001,1]>MatrixDel[2001,2]){
	transB[2001,1]="F"}
if (MatrixDel[2001,1]<=MatrixDel[2001,2]){
	transB[2001,1]="U"}

	#algo: remonter jusqu'a la position 1
for (i in 1:2000){
	if (transB[(2002-i),1]=="F"){
		cA=MatrixDel[(2001-i),1]+P[1]
		cB=MatrixDel[(2001-i),2]+P[2]
		if (cA>cB){
			transB[(2001-i),1]="F"	
		}
	else{
			transB[(2001-i),1]="U"	
		}
	}
	if (transB[(2002-i),1]=="U"){
		cA=MatrixDel[(2001-i),1]+P[3]
		cB=MatrixDel[(2001-i),2]+P[4]
		if (cA<cB){
			transB[(2001-i),1]="U"
	}
	else{
			transB[(2001-i),1]="F"	
		}
}}

print(transB)
plot(trans[,1]=="F")
lines(transB[,1]=="F",col="red")


data <-trans[,1]=="F" # création d'un objet "data" contenant les valeurs 
jpeg("mongraph.jpg")
plot(data)
lines(transB[,1]=="F",col="red")
dev.off() # permet d'ouvrir le fichier .jpg sans fermer R
#XXXXXXXXXXXXXXXXXXXXXXXX

