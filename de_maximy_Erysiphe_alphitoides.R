rm(list=objects()); graphics.off()


library(cobiclust)
library(ggplot2)
library(purrr)
library(dplyr)

### ---- Données

counts = read.table(file = "https://mia-paris.pages.mia.inra.fr/formation_abondance_reseau/tutoriels/PLN_TP/Data/counts.tsv")
metadata = read.table(file = "https://mia-paris.pages.mia.inra.fr/formation_abondance_reseau/tutoriels/PLN_TP/Data/metadata.tsv")

### ---- Préparation et quelques résultats sur la donnée

counts = as(t(counts), "matrix") # on transforme le dataframe en matrice et on transpose car on veut les OTUs sur les lignes et les feuilles sur les colonnes

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
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )



# Calcul de la matrice nu_j d'offsets 

calculate_offset = function(counts = counts) {
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

mat_nuj = calculate_offset(counts)

length(which(counts==0))*100/length(counts) # pourcentage de données nulles 
cat("Le jeu de données présente ",length(which(counts==0))*100/length(counts), "% de données nulle" )

### ---- Essai de plusieurs Poisson-Gamma LBM

Kmax = 6 # nombre de cluster pour les OTUs max
Gmax = 6 # nombre de cluster max pour les feuilles
arg2 = rep(rep(1:Kmax), Gmax) 
arg3 = rep(1:Gmax, each = Kmax) # pour former toutes les paires (K,G) de combinaisons possibles pour K,G <= Kmax, Gmax resp


# Cas hétérogène akg = TRUE : le paramètre de forme / de dispersion est différent pour chaque groupe

res_akg =  pmap(list(arg2, arg3), ~ cobiclust(counts, ..1, ..2, nu_j = mat_nuj, akg = TRUE))
# tous les cobiclusts pour la matrice de comptage, pour tous les K et G, avec mat_nuj comme matrice d'offset
# akg = TRUE 

# Cas homogène akg = FALSE: le paramètre de forme est le même pour chaque groupe

res_a = pmap(list(arg2, arg3), ~ cobiclust(counts, ..1, ..2, nu_j = mat_nuj, akg = FALSE)) # résultat sous la forme d'une liste à Kmax*Gmax = 36 élements contenant le résultat du cobiclust, la valeur de K et la valeur de G


### ---- Séléction du meilleure modèle

arg1 = 1:(Kmax*Gmax)

# Calcul des critères BIC et ICL variationnels
# selection_criteria = fonction du package donnant le BIC , le vICL, la pénalité associée à K et G ,
# la lower bound et l'entropie

crit_akg = pmap(list(arg1, arg2, arg3), ~ selection_criteria(res_akg[[..1]], K = ..2, G = ..3), )
crit_akg = data.frame(Reduce(rbind, crit_akg)) 

crit_a = pmap(list(arg1, arg2, arg3), ~ selection_criteria(res_a[[..1]], K = ..2, G = ..3), )
crit_a = data.frame(Reduce(rbind, crit_a))


# Séléction

vICL_kg = crit_akg[which.max(crit_akg$vICL), ]
vICL_kg
ou_best_kg = apply( vICL_kg[,c(1,6:7)], 1, FUN=function(x) intersect(which(arg2==x[2]),which(arg3==x[3]))) # on trouve où est atteint le meilleure vICL
best_akg = res_akg[[ou_best_kg]]

cat("Dans le cas hétérogène le meilleur vICL est atteint pour", best_akg$K, "clusters d'OTUs et", best_akg$G, "clusters de feuilles")

ggplot(crit_akg, aes(x = K, y = G, fill = vICL)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap du vICL pour différentes combinaisons (K, G)",
       x = "Nombre de clusters (OTUs, K)",
       y = "Nombre de clusters (Feuilles, G)") +
  theme_minimal()


vICL = crit_a[which.max(crit_a$vICL),]
vICL
ou_best = apply( vICL[,c(1,6:7)], 1, FUN=function(x) intersect(which(arg2==x[2]),which(arg3==x[3])))
best_a = res_a[[ou_best]]

BIC = crit_a[which.min(crit_a$BIC)]

cat("Dans le cas homogène le meilleur vICL est atteint pour", best_a$K, "clusters d'OTUs et", best_a$G, "clusters de feuilles")

ggplot(crit_a, aes(x = K, y = G, fill = vICL)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap du vICL pour différentes combinaisons (K, G)",
       x = "Nombre de clusters (OTUs, K)",
       y = "Nombre de clusters (Feuilles, G)") +
  theme_minimal()


# Graphique vICL en fonction du nombre de cluster

crit_a <- crit_a %>% mutate(Model = "a")
crit_akg <- crit_akg %>% mutate(Model = "akg")
crit_combined <- bind_rows(crit_a, crit_akg) # on fusionne les deux tableaux de critère en ajoutant une colonne pour préciser "a" s'il s'agit du modèle homogène ou "akg" sinon

# Créer les graphiques : combinés, modèle hétérogène et homogène
ggplot(crit_combined, aes(x = K, y = vICL, color = factor(G), shape = Model, group = interaction(G, Model))) +
  geom_line() +  # Tracer les lignes
  geom_point(size = 3) +  # Ajouter des points
  labs(
    title = "Comparaison des critères vICL pour homogène et hétérogène",
    x = "Nombre de clusters pour les OTUs (K)",
    y = "vICL",
    color = "Nombre de clusters (G)",
    shape = "Modèle"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot(crit_akg, aes(x = K, y = vICL, color = factor(G), shape = Model, group = interaction(G, Model))) +
  geom_line() +  # Tracer les lignes
  geom_point(size = 3) +  # Ajouter des points
  labs(
    title = "Comparaison des critères vICL pour hétèrogène",
    x = "Nombre de clusters pour les OTUs (K)",
    y = "vICL",
    color = "Nombre de clusters (G)",
    shape = "Modèle"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggplot(crit_a, aes(x = K, y = vICL, color = factor(G), shape = Model, group = interaction(G, Model))) +
  geom_line() +  # Tracer les lignes
  geom_point(size = 3) +  # Ajouter des points
  labs(
    title = "Comparaison des critères vICL pour homogène",
    x = "Nombre de clusters pour les OTUs (K)",
    y = "vICL",
    color = "Nombre de clusters (G)",
    shape = "Modèle"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# c'est un peu brouillon les points sont très proches 

### --- Pour le modèle hétérogène choisi :

# Groupes de lignes et colonnes

Z_kg = best_akg$classification$rowclass # clusters associés aux lignes / OTUs
W_kg = best_akg$classification$colclass # clusters associés aux colonnes / feuilles 
table(Z_kg) # taille de chaque cluster
table(W_kg) 

# Descriptions des groupes

sapply(1:max(W_kg),FUN=function(x) table(colnames(counts)[W_kg==x])) # description des clusters de bactéries / champignons
sapply(1:max(Z_kg),FUN=function(x) rownames(counts)[Z_kg==x]) # description des clusters de feuille -> le 1 correspond à celui avec le champignon E_alphitoides

table(data.frame(metadata$tree,W_kg)) # croisement entre le type d'arbre (résistant, intermédiaire, susceptible) et les clusters

### --- Pour le modèle homogène choisi : 

# Groupes de lignes et colonnes

Z = best_a$classification$rowclass
W = best_a$classification$colclass
table(Z)
table(W)

# Description des groupes 

sapply(1:max(W),FUN = function(x) table(colnames(counts)[W==x]))
sapply(1:max(Z),FUN = function(x) rownames(counts)[Z==x]) # cluster 1 avec champignon

table(data.frame(metadata$tree, W))


### Paramètres du modèle

# alpha_kg
print(best_akg$parameters$alpha)
print(best_a$parameters$alpha)
print(min(best_a$parameters$alpha))
print(max(best_a$parameters$alpha))


# a
print(best_akg$parameters$a)
print(mean(best_akg$parameters$a))
print(best_a$parameters$a)

### Graphiques

### Indice de diversité de Shannon par rapport à W


# heterogène
don = counts
p_i = apply(don, 2, FUN = function(x) x / sum(x)) # proportion relatives de chaque OTU dans chaque échantillon
H = apply(p_i, 2, FUN = function(x) - sum(x * log2(x), na.rm = TRUE)) # calcul de l'indice de Shannon pour chaque échantillon
#-> élevé = diversité élevée (proportions équilibrées d'espèce) entre 1 et 5
df_kg = data.frame(Leaf = 1:ncol(don), H = H, W = factor(round(best_akg$info$t_jg %*% (1:max(W_kg))))) # t_jg = probabilité d'appartenance des feuilles aux clusters d'OTUs
levels(df_kg$W) = paste("W",levels(df_kg$W),sep="")
hwplot = ggplot(data = df_kg, aes(x = W, y = H))
hwplot = ggplot(data = df_kg, aes(x = W, y = H)) +
  geom_boxplot() +
  labs(
    x = "",  
    y = "Indice de diversité de Shannon",
    title = "Distribution de la diversité de Shannon par cluster (W) : cas hétérogène"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")   
  )

print(hwplot)

rm(p_i, df_kg)


#homogène 
df = data.frame(Leaf = 1:ncol(don), H = H, W = factor(round(best_a$info$t_jg %*% (1:max(W)))))
levels(df$W) = paste("W",levels(df$W),sep="")
hwplot = ggplot(data = df, aes(x = W, y = H))
hwplot = ggplot(data = df, aes(x = W, y = H)) +
  geom_boxplot() +
  labs(
    x = "",  # Pas de label pour l'axe X
    y = "Indice de diversité de Shannon",
    title = "Distribution de la diversité de Shannon par cluster (W) : cas homogène"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")   
  )

print(hwplot)

rm(df)


### Niveau d'infection par groupe de feuilles

# hétérogène

plot.infection <- ggplot(
  data = data.frame(
    pmInfection = log2(metadata$pmInfection), 
    W = as.factor(paste("W", W_kg, sep = ""))
  ), 
  aes(y = pmInfection, x = W)
) +
  geom_boxplot() +
  labs(
    title = "Distribution du niveau d'infection par cluster (W) : cas hétérogène",  
    x = "", 
    y = "Level of infection"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")   
  )

print(plot.infection)

# homogène

plot.infection = ggplot(
  data = data.frame( pmInfection = log2(metadata$pmInfection),  W = as.factor(paste("W", W, sep = "")) ), 
  aes(y = pmInfection, x = W)) +
  geom_boxplot() +
  labs( title = "Distribution du niveau d'infection par cluster (W) : cas homogène",  
    x = "", 
    y = "Level of infection"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")   
  )

print(plot.infection)


### mu chapeau i dans chaque cluster ligne Zk

# cas hétérogène
df1 = data.frame(Microorg = rownames(counts), Z = round(best_akg$info$s_ik %*% (1:max(Z_kg))))
df2 = data.frame(Microorg = rownames(counts), Mu = best_akg$parameters$mu_i)
df = merge(df1, df2, by = "Microorg")
tmp = sapply(1:max(Z_kg), FUN = function(x) mean(best_akg$parameters$mu_i[Z == x]))
df$Z = factor(df$Z)
levels(df$Z) = paste("Z", levels(df$Z), sep = "")

muzplot <- ggplot(data = df, aes(x = Z, y = log(Mu))) +
  geom_boxplot() +
  labs(
    title = "Distribution des valeurs log(Mu) par cluster (Z) : cas hétérogène",
    x = "", 
    y = expression(mu[i])
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")  
  )

print(muzplot)


rm(df1, df2, df, tmp)

# cas homogène
df1 = data.frame(Microorg = rownames(counts), Z = round(best_akg$info$s_ik %*% (1:max(Z))))
df2 = data.frame(Microorg = rownames(counts), Mu = best_a$parameters$mu_i)
df = merge(df1, df2, by = "Microorg")
tmp = sapply(1:max(Z), FUN = function(x) mean(best_a$parameters$mu_i[Z == x]))
df$Z = factor(df$Z)
levels(df$Z) = paste("Z", levels(df$Z), sep = "")

muzplot <- ggplot(data = df, aes(x = Z, y = log(Mu))) +
  geom_boxplot() +
  labs(
    title = "Distribution des valeurs log(Mu) par cluster (Z) : cas homogène",
    x = "", 
    y = expression(mu[i])
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")  
  )

print(muzplot)

rm(df1, df2, df, tmp)

### log ratio des interactions entre les cluster et le cluster le plus contaminé

alpha = best_a$parameters$alpha 


W1 = 1  
W2 = 2 
W3 = 3  
W4 = 4

log_ratios_W3_W1 <- log(alpha[, W3]) - log(alpha[, W1])
log_ratios_W3_W2 <- log(alpha[, W3]) - log(alpha[, W2])
log_ratios_W3_W4 = log(alpha[,W3]) - log(alpha[, W4])

log_ratios_df <- data.frame(
  OTU_Cluster = 1:nrow(alpha),  
  Log_Ratio_W3_W1 = log_ratios_W3_W1,
  Log_Ratio_W3_W2 = log_ratios_W3_W2,
  Log_Ratio_W3_W4 = log_ratios_W3_W4
)

log_ratios_plot <- ggplot(log_ratios_df, aes(x = OTU_Cluster)) +
  geom_line(aes(y = Log_Ratio_W3_W1, color = "W3 vs W1"), size = 1) +
  geom_line(aes(y = Log_Ratio_W3_W2, color = "W3 vs W2"), size = 1) +
  geom_line(aes(y = Log_Ratio_W3_W4, color = "W3 vs W4"), size = 1) +
  labs(
    title = "Log-ratios des termes d'interaction par cluster d'OTUs",
    x = "Cluster d'OTUs (K)",
    y = "Log-ratio",
    color = "Comparaison"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("W3 vs W1" = "blue", "W3 vs W2" = "red", "W3 vs W4" = "orange"))

print(log_ratios_plot)




