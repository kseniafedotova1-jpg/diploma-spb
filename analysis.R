# =============================================================================
# Анализ инвестиционной привлекательности регионов РФ
# =============================================================================

library(readxl); library(tidyverse); library(dplyr)
library(pastecs); library(psych); library(corrplot); library(writexl)
library(clustertend); library(hopkins); library(factoextra)
library(cluster); library(fpc); library(dbscan); library(NbClust); library(MASS)
library(sf); library(ggplot2); library(tmap); library(tmaptools)
library(RColorBrewer); library(leaflet); library(summarytools)

# =============================================================================
# 1. ЗАГРУЗКА ДАННЫХ
# =============================================================================
df <- read_excel("data/spb_database.xlsx", sheet = "кластерный 2023")
df <- mutate_at(df, vars(2:24), as.numeric)
df1 <- df[, 2:24]
df1 <- mutate_at(df1, vars(1:23), as.numeric)

# =============================================================================
# 2. ОПИСАТЕЛЬНАЯ СТАТИСТИКА И КОРРЕЛЯЦИЯ
# =============================================================================
df1stat <- stat.desc(df1)
cordf1 <- round(cor(df1), 2)

# =============================================================================
# 3. ПРОВЕРКА УСЛОВИЙ ФАКТОРНОГО АНАЛИЗА
# =============================================================================
KMO(df1)
cortest.bartlett(cor(df1), n = nrow(df1))
summary(df1)
apply(df1, 2, sd)

# =============================================================================
# 4. МЕТОД ГЛАВНЫХ КОМПОНЕНТ (PCA)
# =============================================================================
pca <- prcomp(df1, scale. = TRUE)
screeplot(pca, type = "lines", main = "Scree Plot")
summary(pca)

ggplot(NULL, aes(names(summary(pca)$importance[2, 1:9]),
                 summary(pca)$importance[2, 1:9])) +
  geom_point(size = 2) + geom_line(aes(group = 1)) +
  xlab("Principal components") + ylab("Variance Proportion")

PCA1 <- pca$rotation[, 1:3]

pz_pca <- tibble(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  region = df$region
)

ggplot(pz_pca, aes(x = PC1, y = PC2, label = region)) +
  geom_point() + geom_text(size = 3, vjust = -0.5) + theme_minimal() +
  labs(title = "Распределение регионов по обобщающим факторам (PC1 и PC2)",
       x = "PC1", y = "PC2")

# =============================================================================
# 5. КЛАСТЕРНЫЙ АНАЛИЗ
# =============================================================================

# Исключаем регионы-выбросы
outliers <- c("г. Москва", "Московская область", "г. Санкт-Петербург",
              "Ямало-Ненецкий автономный округ",
              "Ненецкий автономный округ", "Чукотский автономный округ")

df_clean   <- pz_pca[!pz_pca$region %in% outliers, ]
df_claster <- df_clean[, 1:3]

corddf_claster <- round(cor(df_claster), 2)
cluststat <- format(stat.desc(df_claster), scientific = FALSE)
hopkins(df_claster)

# Определение оптимального числа кластеров
fviz_nbclust(df_claster, kmeans, method = "wss")
fviz_nbclust(df_claster, kmeans, method = "silhouette")
fviz_nbclust(df_claster, kmeans, method = "gap_stat")

res <- NbClust(df_claster, min.nc = 2, max.nc = 9, method = "kmeans")
res$Best.nc

barplot(table(res$Best.nc[1, ]),
        xlab = "Число кластеров", ylab = "Количество индексов",
        col = "skyblue", main = "Выбор числа кластеров (NbClust)")

dist.dat  <- dist(df_claster)
clust.dat <- hclust(dist.dat, "ward.D")
plot(clust.dat)

# K-means (2 кластера)
kmeansclasters <- kmeans(df_claster, 2)
df_clean$clust <- kmeansclasters$cluster
centers        <- as.data.frame(kmeansclasters$centers)

ggplot(df_clean, aes(x = PC1, y = PC2, color = factor(clust))) +
  geom_point(size = 2) +
  geom_point(data = centers, aes(x = PC1, y = PC2), color = "black", size = 4, shape = 8) +
  labs(title = "Кластеры и их центры (PC1 vs PC2)", color = "Кластер") + theme_light()

ggplot(df_clean, aes(x = PC2, y = PC3, color = factor(clust))) +
  geom_point(size = 2) +
  geom_point(data = centers, aes(x = PC2, y = PC3), color = "black", size = 4, shape = 8) +
  labs(title = "Кластеры и их центры (PC2 vs PC3)", color = "Кластер") + theme_light()

ggplot(df_clean, aes(x = PC1, y = PC3, color = factor(clust))) +
  geom_point(size = 2) +
  geom_point(data = centers, aes(x = PC1, y = PC3), color = "black", size = 4, shape = 8) +
  labs(title = "Кластеры и их центры (PC1 vs PC3)", color = "Кластер") + theme_light()

c1 <- filter(df_clean, clust == 1)
c2 <- filter(df_clean, clust == 2)
a1 <- format(stat.desc(c1), scientific = FALSE)
a2 <- format(stat.desc(c2), scientific = FALSE)

# =============================================================================
# 6. ОЦЕНКА КАЧЕСТВА КЛАСТЕРИЗАЦИИ
# =============================================================================
sil      <- silhouette(as.numeric(df_clean$clust), dist.dat)
mean_sil <- mean(sil[, 3])
cat("Средняя ширина силуэта:", round(mean_sil, 3), "\n")

plot(sil, col = 2:(max(as.numeric(df_clean$clust)) + 1), border = NA)
abline(v = round(mean_sil, 3), col = "blue", lty = 2, lwd = 2)

stats <- cluster.stats(dist.dat, df_clean$clust)
cat("Среднее внутрикластерное:", stats$average.within, "\n")
cat("Среднее межкластерное:   ", stats$average.between, "\n")
cat("Индекс Данна:            ", stats$dunn, "\n")
