##############################################################################################
#Fonction 1 
Prim=function(A,X=rep(1:nrow(A))){
  
  # choix d'un sommet au hasard
  Xp <- sample(X,1) 
  
  # lignes et colonnes des poids positifs
  pweight=which(A>0,arr.ind=T)
  
  # matrice U  = arêtes j-k de poids minimal + poids associé
  U <- matrix(0,length(X)-1,3)
  
  # compteur qui va servir à modifier la matrice U 
  cpt = 1 
  
  while (all(intersect(X,Xp)==X)==FALSE){ 
    # génère des warnings a cause des tailles différentes entre Xp et X,
    # ne pas en tenir compte
    
    Xb <- setdiff(X,Xp)
    
    
    # sommets J,K et P les poids de la matrice d'adjacence
    J = K = P = c() 
    
    # ind =  sommets qui sont dans Xp et les sommets qui sont dans Xb
    ind=which((pweight[,1]%in%Xp)&(pweight[,2]%in%Xb)) 
    
    for (i in ind){
      j <- pweight[i,1] # sommet j qui vérifie j appartient à Xp
      k <- pweight[i,2] # sommet k qui vérifie k appartient à Xb
      J <- c(J,j)       # ensemble des sommets j 
      K <- c(K,k)       # ensemble des sommets k
      P <- c(P,A[j,k])  # ensemble des poids
      poids_min=min(P)  # on récupère le poids minimal
      # on prend le 1er indice du poids minimal trouvé parmis les poids P
      ind_poids_min= which.min(P) 
    }
    Xp <- c(Xp,K[ind_poids_min])                            
    U[cpt,] <- c(J[ind_poids_min],K[ind_poids_min],poids_min)
    cpt=cpt+1
  }
  return(U)
}

##############################################################################################
#Fonction2

Bellman=function(X,A,s){
  k=c()
  # pi = vecteur qui représente les étiquettes de chaque sommet
  pi=rep(Inf,length(X))
  pi[s]=0
  
  # Tant qu'une étiquette du vecteur pi change
  while(setequal(k,pi)==FALSE){
    k=pi # on mémorise le dernier vecteur pi connu afin de pouvoir 
    # le comparer au prochain qui aura été modifié ou non
    
    for (i in setdiff(X,s)){
      # mise à jour des étiquettes
      pi[i]=min(pi[i],(pi[which(A[,i]!=0)]+c(A[which(A[,i]!=0),i]))) 
    }
  }
  return(pi)
}

##############################################################################################
#Fonction3

Ford_Fulkerson = function(X,A,s,p,phi,v_phi){
  
  
  C1 <- (A-phi)>0 # 1ère condition
  C2 <- t(phi)>0 # 2ème condition
  
  C=C1|C2 # matrice de booléens qui vérifie C1 ou C2
  
  # on commence par mettre dans S la source s 
  S = c(s) 
  
  # on associe à une liste m_s les triplés qui caractérisent la source
  m_s = list(0,Inf,"+")
  
  # on crée une liste m qui va contenir les listes des marquages permanents 
  m=list()
  m[[s]]=m_s # ajout de m_s dans m en position s
  
  
  Sb = setdiff(X,S) # ensemble des sommets X qui ne sont pas dans S
  
  # indices des lignes et colonnes pour lesquels les valeurs 
  # de la matrice C[S,Sb] sont à TRUE
  ind=which(matrix(C[S,Sb],nrow=length(S),ncol = length(Sb)),arr.ind = TRUE) 
  
  S_Sb=cbind(S[ind[,1]],Sb[ind[,2]]) 
  # ind[,1] : 1ère colonne de ind correspond aux sommets de S
  # ind[,2] : 2ème colonne de ind, correspond aux sommets de Sb
  
  # Tant que la matrice S_Sb n'est pas vide
  while (all(is.na(S_Sb))!=TRUE) {
    
    # on prend le premier couple (i,j) de la matrice S_Sb qui vérifient C
    # (on aurait pu prendre n'importe quel couple de cette matrice)
    i <- S_Sb[1,1]
    j <- S_Sb[1,2]
    
    alpha_i <- m[[i]][[2]]
    
    if ((A[i,j]-phi[i,j])>0){
      alpha_j <- min(alpha_i,(A[i,j]-phi[i,j]))
      m[[j]] <- list(i,alpha_j,'+')
    }
    
    else if (phi[j,i]>0){
      alpha_j <- min(alpha_i,phi[j,i])
      m[[j]] <- list(i,alpha_j,'-')
    }
    
    # on met a jour S, Sb et S_Sb avant de rentrer à nouveau dans la boucle
    S=c(S,j)
    Sb = setdiff(X,S)
    ind=which(matrix(C[S,Sb],nrow=length(S),ncol=length(Sb)),arr.ind = TRUE)
    S_Sb=cbind(S[ind[,1]],Sb[ind[,2]])
    
    # Si j correspond au puit
    if (j==p){ 
      v_phi=v_phi+alpha_j
      break
    }
  }
  
  # si le puit se trouve dans S..., lignes exécutées si il y a eu un break 
  if (p %in% S){
    # tant que j est différent de la source
    while(j!=s){ 
      
      if (m[[j]][3]=='+'){
        phi[unlist(m[[j]][1]),j] <- phi[unlist(m[[j]][1]),j]+alpha_j
      }
      
      else if (m[[j]][3]=='-'){
        phi[j,unlist(m[[j]][1])] <- phi[j,(unlist(m[[j]][1]))]-alpha_j
      }
      
      j<- unlist(m[[j]][1])
    }
    # Aller en 1 du pseudo-code -> on rappelle la fonction avec X,A,s et p
    # qui restent inchangés et les nouvelles valeurs de phi et v_phi
    return (Ford_Fulkerson(X,A,s,p,phi,v_phi))
  }
  else{
    # Fin du programme, on retourne le flot maximal
    return(v_phi)
  }
}
##############################################################################################