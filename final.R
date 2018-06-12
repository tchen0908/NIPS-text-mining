##### Final project K-means & EM #####
load("C:/Users/nohah/Desktop/Sessions/NCCU_STAT/106-2/統計計算與模擬/Bayes/final/workspace_final.RData")

########################################
##### 清理term #####
## build new matrix
## asgin the old to new
## find those haven't be used
## add to the new matrix

##### 分析 #####
## 看每個年份有多少篇論文  <OK>
## 清理完之後呢 轉成TF-IDF <OK> 
## 因為還是sparse matrix 做完SVD分解可以降低維度
## 找出頻率高的前五百個字詞作為feature
## 做分群先把同一年的放一起

##### 視覺化 #####
## 文字雲  <OK>
## 不同年份的wordcloud  
########################################
### packages ###
library(NLP)
library(tm)
library(rvest)
library(dplyr)
library(SnowballC)
library(tidytext)
library(wordcloud2)
library(magrittr)
library(textstem)
library(plyr)
library(ggplot2)
library(cluster)

NIPS <- read.csv("NIPS_1987-2015.csv")




### different years term freq -----------------------------------
cnt <- NULL
for(i in 1987:2015){
  a <- which(substr(colnames(NIPS.cleaned),2,5) == as.factor(i)) %>% length()
  cnt <- c(cnt,a)
}
cnt
year.count <- data.frame(1987:2015,cnt)
colnames(year.count) <- c("year","count")
attach(year.count)

df <- data.frame(year = year.count$year, count = year.count$count)
ggplot(df, aes(year, count)) +
  geom_col()+
  theme(axis.text.x = element_text(face = "bold",vjust = 0, angle = 45))





### wordcloud --------------------------------------------------
# build term frequency data frame
(freq <- rowSums(as.matrix(NIPS[,-1])) %>% as.matrix())
freq.frame <- data.frame(NIPS[,1],freq) 
colnames(freq.frame) <- c("word","freq")
freq.frame  


# sort the freq descendingly
word.count <- freq.frame %>% arrange(desc(freq))

# plot wordcloud 
wordcloud2(freq.frame)
wordcloud2(word.count[1:100,],size = .3,color = "random-dark",shape = "star")

# plot wordcloud for 1987-1995/1996-2005/2006-2015
which(substr(colnames(NIPS[,]),2,5) == "2015")

NIPS.87_95 <- NIPS[,2:1137]
NIPS.96_05 <- NIPS[,1138:2741]
NIPS.06_15 <- NIPS[,2742:5812]

tf.wordcloud <- function(doc){
  freq <- rowSums(as.matrix(doc)) %>% as.matrix()
  freq.frame <- data.frame(NIPS[,1],freq)
  colnames(freq.frame) <- c("word","freq")
  wordcloud2(freq.frame[1:500,], size = .5,color = "random-dark",shape = "star") 
}

tf.wordcloud(NIPS.87_95)
tf.wordcloud(NIPS.96_05)
tf.wordcloud(NIPS.06_15)





### TF-IDF -----------------------------------------------------
NIPS.fortfidf <- NIPS
NIPS.fortfidf[,1] %<>% as.character()
rownames(NIPS.fortfidf) <- NIPS[,1]
NIPS.fortfidf <- NIPS.fortfidf[,-1]
NIPS.fortfidf %<>% as.matrix()
colnames(NIPS.fortfidf)

# 轉成TF-IDF矩陣
tdm.tfidf <- as.TermDocumentMatrix( NIPS.fortfidf,weighting = weightTfIdf)
tdm.tfidf
tdm.tf <- as.TermDocumentMatrix(NIPS.fortfidf, weighting = weightTf)
tdm.tf

inspect(tdm.tfidf[5:10,5:10])

# total sum of TF-IDF or TF
freq <- rowSums(as.matrix(tdm.tfidf))

# Plot those frequencies ordered.
plot(sort(freq, decreasing = T),col="darkorange",main="Word TF-IDF frequencies", xlab="TF-IDF-based rank", ylab = "TF-IDF")

# see the most frequent terms
tail(sort(freq),n=10)

# Show most frequent terms and their frequencies in a bar plot.
high.freq=tail(sort(freq),n=10)
hfp.df=as.data.frame(sort(high.freq))
hfp.df$names <- rownames(hfp.df) 

ggplot(hfp.df, aes(reorder(names,high.freq), high.freq)) +
  geom_bar(stat="identity") + coord_flip() + 
  xlab("Terms") + ylab("Term Frequency") +
  ggtitle("TF Top 10")

# 看出至少10次的term
findFreqTerms(tdm, 10)

# find correlated terms
findAssocs(tdm.tfidf, "network", 0.5)
findAssocs(tdm.tfidf, "spike", 0.5)
findAssocs(tdm.tfidf, "units", 0.5)

# remove 40 percent of sparse
inspect(removeSparseTerms(tdm, 0.4))






### cleaning ---------------------------------------------------
# select top 2000 terms by TF-IDF
tfidf.sum <- rowSums(as.matrix(tdm.tfidf)) %>% as.matrix
tfidf.frame <- data.frame(NIPS[,1],tfidf.sum) 
colnames(tfidf.frame) <- c("word","tfidf")
tfidf.frame

top2000.tfidf <- tfidf.frame %>% arrange(desc(tfidf))
top2000.tfidf[1900:2001,]  # the 2001th TF-IDF is 1.833669, so choose those more than 1.8337

NIPS.tfidfselected <- NIPS[which(tfidf.frame$tfidf > 1.8337),]

NIPS.test <- NIPS.tfidfselected
NIPS.test[,1] %<>% as.character() 
NIPS.test[,1] %<>%  lemmatize_words()

a <- c("raiply","rapid","badly")
lemmatize_words(a)
# combine the rows that have same terms
NIPS.cleaned <- ddply(NIPS.test,"X",numcolwise(sum))  # 要跑很久小心
NIPS.cleaned <- NIPS.cleaned[-c(1,2),]

# new TF-IDF matrix
NIPS.fortfidf2 <- NIPS.cleaned
NIPS.fortfidf2[,1] %<>% as.character()
rownames(NIPS.fortfidf2) <- NIPS.cleaned[,1]
NIPS.fortfidf2 <- NIPS.fortfidf2[,-1]
NIPS.fortfidf2 %<>% as.matrix()

# 轉成TF-IDF矩陣
tdm.tfidf_clean <- as.TermDocumentMatrix( NIPS.fortfidf2,weighting = weightTfIdf)
tdm.tfidf_clean
tdm.tf_clean <- as.TermDocumentMatrix(NIPS.fortfidf2, weighting = weightTf)
tdm.tf_clean



### new TF-IDF
# total sum of TF-IDF or TF
freq <- rowSums(as.matrix(tdm.tfidf_clean))

# Plot those frequencies ordered.
plot(sort(freq, decreasing = T),col="darkorange",main="Word TF-IDF frequencies", xlab="TF-IDF-based rank", ylab = "TF-IDF")

# see the most frequent terms
tail(sort(freq),n=10)

# Show most frequent terms and their frequencies in a bar plot.
high.freq=tail(sort(freq),n=10)
hfp.df=as.data.frame(sort(high.freq))
hfp.df$names <- rownames(hfp.df) 

ggplot(hfp.df, aes(reorder(names,high.freq), high.freq)) +
  geom_bar(stat="identity") + coord_flip() + 
  xlab("Terms") + ylab("TF-IDF") +
  ggtitle("Cleaned Data TF-IDF Top 10")

### wordcloud for new TF-IDF
# build TF-IDF data frame
(tfidf <- rowSums(as.matrix(tdm.tfidf_clean)) %>% as.matrix())
tfidf.frame <- data.frame(rownames(a),tfidf) 
colnames(tfidf.frame) <- c("word","tf-idf")

# sort the freq descendingly
word.count <- tfidf.frame %>% arrange(desc(tfidf))

# plot wordcloud by TF-IDF of top 100
wordcloud2(tfidf.frame)
wordcloud2(word.count[1:100,],size = .3,color = "random-dark",shape = "star")






### clustering -------------------------------------------------
## cosine distance matrix
library(slam)
cosine_dist_mat <- crossprod_simple_triplet_matrix(tdm.tfidf_clean)/(sqrt(col_sums(tdm.tfidf_clean^2) %*% t(col_sums(tdm.tfidf_clean^2))))

a <- cosine_dist_mat
a[which(is.nan(a))] <- 0

# inspect the sorted correlation of documents
aa <- a[do.call(order, as.data.frame(a)),] 
aa[5792:5811,1:2]  # see top 20 documents related to X1987_1



## HC clustering
a <- removeSparseTerms(tdm.tfidf_clean, sparse = 0.3)
aa <- as.matrix(a)
distMatrix <- dist(scale(aa))
a <- hclust(distMatrix, method = "ward.D")
plot(a, cex=0.9,hang=-1,main="Term Cluster Dendrogram")
rect.hclust(a,k=3)

b <- removeSparseTerms(t(tdm.tfidf_clean), sparse = 0.3)
bb <- as.matrix(b)
distMatrix <- dist(scale(bb))
b <- hclust(distMatrix, method = "ward.D")
plot(b, cex=1,hang=-1,main="Document Cluster Dendrogram")
rect.hclust(b,k=5)

## K-means
kmeans.tfidfclean <- kmeans(aa,3,nstart = 20)
kmeans.tfidfclean$cluster <- as.factor(kmeans.tfidfclean$cluster)
ggplot(aa, aes(Petal.Length, Petal.Width, color = iris$cluster)) + geom_point()

