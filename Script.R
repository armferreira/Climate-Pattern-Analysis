library(dplyr)
library(mice)
library(lsr)
library(confintr)
library(ggcorrplot)
library(psych)
library(factoextra)
library(gpairs)
library(cluster)
library(mclust)
library(ggplot2)
library(sf)
library(rnaturalearth)

############### 1. Dados ###############

# Ler ficheiro csv com nomes
clima <- read.csv("G:/O meu disco/MS Data Science/Pattern Recognition/TrabalhoRP/Datasets/clima.csv", 
                  row.names = 1,
                  stringsAsFactors = T)

############### 1.1 Análise Valores Omissos ############### 

# Analisar valores omissos
md.pattern(clima, rotate.names = T)
na_percentages <- sapply(clima, function(x) sum(is.na(x)) / nrow(clima) * 100)
na_percentages_nonzero <- na_percentages[na_percentages > 0]
na_table <- data.frame(Column = names(na_percentages_nonzero), NA_Percentage = na_percentages_nonzero)

na_percentages_rows <- apply(clima, 1, function(x) sum(is.na(x)) / length(x) * 100)
na_table_rows <- data.frame(Row = seq_along(na_percentages_rows), NA_Percentage = na_percentages_rows)
na_table_rows_nonzero <- na_table_rows[na_table_rows$NA_Percentage > 0, ]

# Remover variavel 'air_transport' por %NAs elevada (~11%)
clima <- select(clima, -air_transport)

############### 1.2 Matriz de correlação ###############

# Função para calcular matriz de correlação entre diferentes tipos de variáveis
calculate_correlation_matrix <- function(data) {
  if (!is.data.frame(data)) {    stop("Input is not a dataframe")  }
  
  matrix <- create_empty_matrix(data)
  variable_combinations <- combn(colnames(data), 2, simplify = FALSE)
  
  for (i in 1:length(variable_combinations)){
    var1 <- variable_combinations[[i]][1]
    var2 <- variable_combinations[[i]][2]
    
    if (is.factor(data[[var1]]) & is.factor(data[[var2]])) {
      c2 <- chisq.test(data[[var1]], data[[var2]])
      cramer_v <- cramersv(c2)
      corr <- round(cramer_v, 4) 
    }
    if (is.numeric(data[[var1]]) & is.numeric(data[[var2]])) {
      corr <- round(cor(data[[var1]], data[[var2]], method = 'pearson', use = "pairwise.complete.obs"), 4)
    }
    if (is.numeric(data[[var1]]) & is.factor(data[[var2]])) {
      model <- aov(data[[var1]] ~ data[[var2]], data = data)      
      eta_result <- etaSquared(model)
      corr <- round(eta_result[1], 4)
    } 
    if (is.factor(data[[var1]]) & is.numeric(data[[var2]])) {
      model <- aov(data[[var2]] ~ data[[var1]], data = data)      
      eta_result <- etaSquared(model)
      corr <- round(eta_result[1], 4)
    }
    row_index1 <- which(rownames(matrix) == var1)
    col_index1 <- which(colnames(matrix) == var2)
    matrix[row_index1, col_index1] <- corr
    
    row_index2 <- which(rownames(matrix) == var2)
    col_index2 <- which(colnames(matrix) == var1)
    matrix[row_index2, col_index2] <- corr
  }
  
  for (i in 1:nrow(matrix)) {
    for (j in 1:nrow(matrix)) {
      row <- rownames(matrix)[i]
      col <- colnames(matrix)[j]
      if (row == col) {
        matrix[i,j] <- 1.00
      }
    }
  }
  return (matrix)
}

create_empty_matrix <- function(dataframe) {
  if (!is.data.frame(dataframe)) {
    stop("Input is not a dataframe")
  }
  num_cols <- ncol(dataframe)
  empty_matrix <- matrix(NA, nrow = num_cols, ncol = num_cols)
  colnames(empty_matrix) <- colnames(dataframe)
  rownames(empty_matrix) <- colnames(dataframe)
  
  return(empty_matrix)
}

# Usar funções para criar matriz de correlação
cor <- round(calculate_correlation_matrix(clima),3)
pred_matrix <- cor
pred_matrix[abs(pred_matrix) > 0.3] <- 1
pred_matrix[abs(pred_matrix) <= 0.3] <- 0
diag(pred_matrix) <- 0

############### 1.3 Imputação com MICE ###############

# Imputar omissos usando MICE
imp <- mice(clima,  m = 5, maxit = 10, method = 'pmm', pred=pred_matrix, seed = 123)

# Ver valores imputados
imp$imp
imp$method

# Avaliar distribuição valores imputados vs observados
stripplot(imp, access_clean~.imp, pch=20, cex=2)
densityplot(imp, ~ access_clean)
stripplot(imp, water_total~.imp, pch=20, cex=2)
densityplot(imp, ~ water_total)
stripplot(imp, fert_cons~.imp, pch=20, cex=2)
# densityplot(imp, ~ fert_cons)
stripplot(imp, fish_total~.imp, pch=20, cex=2)
densityplot(imp, ~ fish_total)
stripplot(imp, renew_water~.imp, pch=20, cex=2)
densityplot(imp, ~ renew_water)
stripplot(imp, sdi~.imp, pch=20, cex=2)
densityplot(imp, ~ sdi)
stripplot(imp, ndgain~.imp, pch=20, cex=2)
densityplot(imp, ~ ndgain)

# Avaliar convergência
plot(imp)

# Extrair o dataset completo
clima_comp <- complete(imp)

############### 2. Análise Componentes Principais ###############

# Excluir variáveis de PROFILE do dataset
clima_input <- clima_comp[,1:21]

# Scatterplot entre variáveis métricas
#pairs(clima_input, pch = 1, lower.panel = NULL)

# Matriz de correlação
corr <- cor(clima_input)
p.mat <- cor_pmat(clima_input)

ggcorrplot(
  corr,
  hc.order = TRUE,
  type = "lower",
  insig = "blank",
  lab = TRUE,
  p.mat = p.mat,
  ggtheme = ggplot2::theme_gray,
)

# Teste de Bartlet
cortest.bartlett(corr)

# Teste KMO
KMO(corr)
## Remover varíaveis com MSA < 0.5
# clima_input <- subset(clima_input, select = -c(food_prod, crop_prod, pop_density, livestock_prod))

# Standardização
clima_scale <- scale(clima_input)

#### PCA
pc10 <- principal(clima_scale, nfactors=10, rotate="none", scores=TRUE)

# Critério de Kaiser - Variâncias dos componentes principais
round(pc10$values, 3)

# Screeplot
plot(pc10$values, type = "b", main = "Scree plot", xlab = "Número CP", ylab = "Eigenvalue")

# Explained Variance
pc10$loadings

# Extrair 5 componentes
pc5 <- principal(clima_scale, nfactors=5, rotate="none", scores=TRUE)

# Fazer a rotação para intepretação
pc5r <- principal(clima_scale, nfactors=5, rotate="varimax")
pc5r$loadings

# Visualizar comunalidades
round(pc5$communality, 2)

# Guardar dataset com novas variáveis
clima_pca <- data.frame(matrix(ncol = 0, nrow = nrow(clima)))
rownames(clima_pca) <- rownames(clima)
clima_pca$pc1 <- pc5$scores[,1]
clima_pca$pc2 <- pc5$scores[,2]
clima_pca$pc3 <- pc5$scores[,3]
clima_pca$pc4 <- pc5$scores[,4]
clima_pca$pc5 <- pc5$scores[,5]

# Biplots
pca <- prcomp(clima_scale)
biplot1 <- fviz_pca_biplot(pca,
                col.ind = "steelblue",
                col.var = "darkred",
                repel = TRUE,
                label = "var",
                labelsize = 4,
                alpha.var = 0.5)
biplot2 <- fviz_pca_biplot(pca,
                geom = "text",
                invisible = "var",
                habillage = clima$continent,
                label = "ind",
                labelsize = 4,
                alpha.ind = "cos2",
                select.ind = list(contrib = 60))

pca_biplot1 <- biplot1 +
  xlab("CP1 (25%)") +
  ylab("CP2 (19%)") +
  labs(title = NULL)
pca_biplot1

pca_biplot2 <- biplot2 +
  xlab("CP1 (25%)") +
  ylab("CP2 (19%)") +
  labs(title = NULL)
pca_biplot2

############### 3. Clustering ############### 

# Matriz de gráficos pareados
gpairs(clima_pca)

######## 3.1. Hierarchical clustering ########
clima_dist <- dist(clima_pca)

# Single-Linkage
hclust_clima <- hclust(clima_dist,method='single')
plot(hclust_clima,label=,hang=-1)

# Complete-Linkage
hclust_compclima <- hclust(clima_dist)
plot(hclust_compclima,label=,hang=-1)

# Ward.D2
hclust_w <- hclust(clima_dist,method='ward.D2')
plot(hclust_w,label=,hang=-1)

groups.k5 <- cutree(hclust_w, k=5)
rect.hclust(hclust_w, k=5, border="darkred")

# Perfilagem com PCAs
aggregate(clima_pca,list(groups.k5), mean)

# Silhueta de Validação
fviz_silhouette(silhouette(groups.k5,clima_dist))

######## 3.2.K-means ########
kmeans.k6 <- kmeans(clima_pca, 6, nstart=100)
attributes(kmeans.k6)

# Tabela de contingencia
table(groups.k5,kmeans.k6$cluster)

# Perfilagem com PCAs
aggregate(clima_pca,list(kmeans.k6$cluster), mean)

# Silhueta de Validação
fviz_silhouette(silhouette(kmeans.k6$cluster,clima_dist))

# Critério "cotovelo
wssplot <- function(xx, nc=15, seed=1234){
  wss <- (nrow(xx)-1)*sum(apply(xx,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(xx, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Número de clusters", ylab="Soma quadrados intra-cluster")
}
wssplot(clima_pca, nc=10)

# Grafico Dispersão + clustering com pontos
plots <- list()
counter <- 1
for(i in 1:4) {
  for(j in (i+1):5) {
    current_vars <- paste0("pc", c(i, j))
    
    p <- fviz_cluster(kmeans.k6, data = clima_pca, stand = FALSE, geom = "point", show.clust.cent = F,
                      ellipse = F, ellipse.type = "convex", choose.vars = current_vars,
                      ggtheme = theme_bw(), xlab = paste0("CP", i),
                      ylab = paste0("CP", j),
                      main = "", ellipse.alpha = 0.1) + theme(legend.position = "none")
    plots[[counter]] <- p
    counter <- counter + 1
  }
}
grid_plot <- do.call(gridExtra::grid.arrange, c(plots, ncol = 4))
grid_plot

# Grafico Dispersão + clustering com texto
fviz_cluster(kmeans.k6, data = clima_pca, stand = F, geom = c("text"), ellipse = F, 
             choose.vars = c("pc1", "pc2"), ggtheme = theme_bw(), xlab = "CP1 (25%)",
             ylab = "CP2 (19%)", main = "", labelsize = 9)

######## 3.3 PAM clustering ######## 
pam.k6 <- pam(clima_pca, 6)

# Tabela de contingencia
table(pam.k6$clustering, kmeans.k6$cluster)

# Perfilagem com PCAs
aggregate(clima_pca,list(pam.k6$clustering), mean)

# Silhueta de Validação
fviz_silhouette(silhouette(pam.k6$cluster,clima_dist))

######## 3.4 Gaussian Mixture Models (GMM) ######## 

# Selecção parâmetros de acordo com BIC
BIC <- mclustBIC(clima_pca)
summary(BIC)

# Estimação de modelo
gmm <- Mclust(clima_pca, x = BIC)
summary(gmm, parameters = TRUE)

# Matriz de classificação, incerteza e densidade
plot(gmm, what = "classification")
plot(gmm, what = "uncertainty")
plot(gmm, what = "density")

# Tabela de contingencia
table(gmm$classification, kmeans.k6$cluster)

# Perfilagem com PCAs
aggregate(clima_pca,list(gmm$classification), mean)

######## 4. Perfilagem ########

# Perfilagem geográfica
world <- ne_countries(scale = "medium", returnclass = "sf")
clusters <- data.frame(gmm$classification)
clusters$country.name <- row.names(clusters)
world_clusters <- merge(world, clusters, by.x = "iso_a3", by.y = "country.name", all.x = TRUE)
colors <- c("#2b2d42", "#ffb703", "lightblue", "lightgreen", "lightcoral")
world_clusters <- world_clusters %>%
  filter(iso_a3 != 'ATA')

ggplot(data = world_clusters) +
  geom_sf(aes(fill = factor(gmm.classification)),  colour = "white") +
  scale_fill_manual(values = colors, name = "Cluster") +
  theme_classic()

# Perfilagem por sdi e ndgain

## Juntar dataframes num só
clima_merged <- merge(clima_pca, clima[,22:25], by = "row.names")
clima_merged <- merge(clima_merged, clusters, by.x = "Row.names", by.y = "country.name", all.x = TRUE)
row.names(clima_merged) <- clima_merged$Row.names
clima_merged <- clima_merged[, -1]
clima_merged$gmm.classification <- as.factor(clima_merged$gmm.classification)

## Estatisticas descritivas
descriptive_stats <- clima_merged %>%
  group_by(gmm.classification) %>%
  summarise(
    mean_sdi = mean(sdi, na.rm = TRUE),
    median_sdi = median(sdi, na.rm = TRUE),
    sd_sdi = sd(sdi, na.rm = TRUE),
    min_sdi = min(sdi, na.rm = TRUE),
    max_sdi = max(sdi, na.rm = TRUE),
    variance_sdi = var(sdi, na.rm = TRUE),
    mean_ndgain = mean(ndgain, na.rm = TRUE),
    median_ndgain = median(ndgain, na.rm = TRUE),
    sd_ndgain = sd(ndgain, na.rm = TRUE),
    min_ndgain = min(ndgain, na.rm = TRUE),
    max_ndgain = max(ndgain, na.rm = TRUE),
    variance_ndgain = var(ndgain, na.rm = TRUE)
  )

## Testes ANOVA
anova_sdi <- aov(sdi ~ gmm.classification, data = clima_merged)
summary(anova_sdi)

anova_ndgain <- aov(ndgain ~ gmm.classification, data = clima_merged)
summary(anova_ndgain)

## Curvas densidade
ggplot(clima_merged, aes(x = sdi, colour = factor(gmm.classification))) + 
  geom_density(size = 1) + 
  labs(x = "Sustainable Development Index (SDI)",
       y = 'Densidade',
       colour = "Cluster") +
  theme_minimal()

ggplot(clima_merged, aes(x = ndgain, colour = factor(gmm.classification))) + 
  geom_density(size = 1) + 
  labs(x = "ND-GAIN",
       y = 'Densidade',
       colour = "Cluster") +
  theme_minimal()

# Perfilagem por paris
## Tabela contingência
table(clima_merged$gmm.classification, clima_merged$paris)

## Teste exato Fischer
fisher.test(table(clima_merged$gmm.classification, clima_merged$paris))

## Gráfico barras paris
ggplot(clima_merged, aes(x = gmm.classification, fill = paris)) +
  geom_bar(position = "dodge") +
  labs(x = "Cluster",
       y = "Freq. Absoluta",
       fill = "Paris") +
       theme_minimal()

# Tabela países por cluster GMM
clusters_countries <- subset(world_clusters, select = c(name, gmm.classification))