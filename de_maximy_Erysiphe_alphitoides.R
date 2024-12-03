rm(list=objects()); graphics.off()


library(cobiclust)
library(ggplot2)
library(purrr)
library(dplyr)

### ---- Données

counts <- read.table(file = "https://mia-paris.pages.mia.inra.fr/formation_abondance_reseau/tutoriels/PLN_TP/Data/counts.tsv")
metadata <- read.table(file = "https://mia-paris.pages.mia.inra.fr/formation_abondance_reseau/tutoriels/PLN_TP/Data/metadata.tsv")

### ---- Préparation et quelques résultats sur la donnée

counts <- as(t(counts), "matrix") # on transforme le dataframe en matrice et on transpose car on veut les OTUs sur les lignes et les feuilles sur les colonnes

OTU = rownames(counts)
Feuilles = colnames(counts)
cat("Il y a", length(OTU), "espèces biologiques différentes dans la donnée et ", length(Feuilles), "échantillons de feuilles dans la donnée.")

entite_bio =  unlist(lapply(strsplit(OTU, split = "_"), FUN = function(y) y[[1]]))
fungi = which(entite_bio %in% c("f","E"))
bact = which(entite_bio == "b")
Erysiphe_alphitoides = which(entite_bio == "E")
cat("Parmi les différentes OTUs de la matrice de comptage, il y a ", length(bact), "bactéries \n
    et ", length(fungi), "espèces de fungi, dont", length(Erysiphe_alphitoides), "champignons infectueux pour l'oïdium du chêne")

# Comparaison de profondeur de séquençage entre les deux différentes entités biologiques (bactéries / fungi)

prof_fungi = colSums(counts[fungi,]) 
prof_bact = colSums(counts[bact,])

sequencing_depth = data.frame(
  Sample = rep(colnames(counts), 2),
  Depth = c(prof_fungi, prof_bact),
  Group = rep(c("Fungi", "Bactéries"), each = ncol(counts))
)


sequencing_depth$Depth = log10(sequencing_depth$Depth)

ggplot(sequencing_depth, aes(x = Group, y = Depth, fill = Group)) +
  geom_boxplot(outlier.color = "grey", outlier.size = 2) +
  labs(
    title = "Profondeur de séquençage (échelle logarithmique)",
    x = "Groupe biologique",
    y = "Log10(Profondeur de séquençage)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Centrer le titre
    legend.position = "none"
  )



# Calcul de la matrice nu_j d'offsets 

calculate_offset_perentite_bio <- function(counts = counts) {
  entite_bio = unlist(lapply(strsplit(rownames(counts), split = "_"), FUN = function(y) y[[1]]))
  
  fungi = which(entite_bio %in% c("f","E")) # c'est un champignon si son nom dans la matrice de comptage est de la forme "f_nomchammpignon" ou "E_champignon" -> lignes correspondants à des fungi
  TC_fungi = colSums(counts[fungi,])/mean(colSums(counts[fungi,])) # somme des valeurs de comptage pour chaque champignon / moyenne de comptage pour tous les champignons
  
  bact = which(entite_bio=="b") # C'est une bactérie lorsque son nom est de la forme "b_OTU" -> ligne correspondant à des bactéries
  TC_bact = colSums(counts[bact,])/mean(colSums(counts[bact,])) # somme des valeurs de comptage pour chaque bactérie / moyenne de comptage pour toutes les bactéries 
  
  # Matrix of nu_j
  indic = c(rep(1, length(which(entite_bio %in% c("f","E")))), rep(0,length(which(entite_bio == "b")))) # vecteur qui contient 1 si la ligne correspond à un fungi et 0 si la ligne correspond à une bactérie
  indic2 = as.matrix(cbind(indic, abs(indic-1))) # on transforme le vecteur en matrice à 2 colonnes (première colonne = champignon : 1 si champi 0 sinon et deuxième colonne = bactérie : 1 si bactérie 0 sinon)
  mat_nuj = indic2%*%matrix(ncol = ncol(counts), nrow = 2, c(TC_fungi, TC_bact), byrow = TRUE)
  # crée une matrice de taille 2 * ncol(counts) au début qui contient l'offset des fungi dans sa première ligne et l'offset des bactéries dans sa seconde
  # on effectue le produit matricielle avec la matrice indic2 pour obtenir une matrice nrow x ncol dans laquelle chaque ligne est remplie avec les offsets correspondant au groupe biologique de cette ligne
  return(mat_nuj)
}

mat_nuj <- calculate_offset_perentite_bio(counts)

length(which(counts==0))*100/length(counts) # pourcentage de données nulles 
cat("Le jeu de données présente ",length(which(counts==0))*100/length(counts), "% de données nulle" )

### Fitting several Poisson-Gamma Latent Block Model

Kmax <- 6
Gmax <- 6
arg2 <- rep(rep(1:Kmax), Gmax)
arg3 <- rep(1:Gmax, each = Kmax)

res_akg <-  pmap(list(arg2, arg3), ~ cobiclust(counts, ..1, ..2, nu_j = mat_nuj, akg = TRUE))

### Selection of the best Poisson-Gamma Latent Block Model

arg1 <- 1:(Kmax*Gmax)
# Calculation of the variational BIC and variational ICL criteria
crit_akg <- pmap(list(arg1, arg2, arg3), ~ selection_criteria(res_akg[[..1]], K = ..2, G = ..3), )
crit_akg <- data.frame(Reduce(rbind, crit_akg)) 
# Selection of the best model
vICL <- crit_akg %>% dplyr::filter(vICL == max(vICL)) 
vICL
ou_best <- apply( vICL[,c(1,6:7)], 1, FUN=function(x) intersect(which(arg2==x[2]),which(arg3==x[3])))
best_akg <- res_akg[[ou_best]]

# We directly fit the best model for this example
best_akg <- cobiclust(counts, K = 5, G = 3, nu_j = mat_nuj, a = NULL, akg = TRUE)

#Some informations
#Groups in rows and column

Z <- best_akg$classification$rowclass
W <- best_akg$classification$colclass
table(Z)
table(W)

# Description of the groups
sapply(1:max(W),FUN=function(x) table(colnames(counts)[W==x]))
sapply(1:max(Z),FUN=function(x) rownames(counts)[Z==x])

table(data.frame(metadata$tree,W))

### Models parameters

# alpha_kg
print(best_akg$parameters$alpha)


# a
print(best_akg$parameters$a)

### Some useful graphics
### Shannon diversity Index according to  𝑊

don <- counts
p_i <- apply(don, 2, FUN = function(x) x / sum(x))
H <- apply(p_i, 2, FUN = function(x) - sum(x * log2(x), na.rm = TRUE))
df <- data.frame(Leaf = 1:ncol(don), H = H, W = factor(round(best_akg$info$t_jg %*% (1:max(W)))))
levels(df$W) <- paste("W",levels(df$W),sep="")
hwplot <- ggplot(data = df, aes(x = W, y = H))
hwplot <- hwplot + geom_boxplot() + labs(x = "")
print(hwplot)

rm(p_i, df)

### Data reordered according to the biclustering
image(log(counts+0.5), main = "Raw data on log-scale")
image(log(counts[order(Z), order(W)]+0.5), main = "Reordered data on log-scale")


### Boxplot of level of infection in log-scale for the groups of leaves

plot.infection <- ggplot(data = data.frame(pmInfection = log2(metadata$pmInfection), W = as.factor(paste("W",W,sep=""))), aes(y = pmInfection, x = W))
plot.infection <- plot.infection + geom_boxplot() + labs(x = "", y = "Level of infection") + theme_bw()
print(plot.infection)


### Boxplot of the mu chapeau i within each group in row 𝑍𝑘

df1 <- data.frame(Microorg = rownames(counts), Z = round(best_akg$info$s_ik %*% (1:max(Z))))
df2 <- data.frame(Microorg = rownames(counts), Mu = best_akg$parameters$mu_i)
df <- merge(df1, df2, by = "Microorg")
tmp <- sapply(1:max(Z), FUN = function(x) mean(best_akg$parameters$mu_i[Z == x]))
df$Z <- factor(df$Z)
levels(df$Z) <- paste("Z", levels(df$Z), sep = "")
muzplot <- ggplot(data = df, aes(x = Z, y = log(Mu)))
muzplot <- muzplot + geom_boxplot() + labs(x = "Z") 
muzplot <- muzplot + theme_bw() +  labs(x = "", y = expression(mu[i])) 
print(muzplot)

rm(df1, df2, df, tmp)

### Plot of the vICL criterion

crit_akg$G <- as.factor(crit_akg$G)
plot.vICL_K <- ggplot(data = crit_akg, aes(y = vICL, x = K, group = G, colour = G))

plot.vICL_K <- plot.vICL_K + geom_point()  + geom_line() +  theme_bw() + theme(legend.position = "bottom") 
print(plot.vICL_K)

