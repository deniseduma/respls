twofour <- read.table("./percentMutation_24Genes.txt",header=TRUE)
names(twofour)[1] <- c('Subject')
twofour <- tbl_df(twofour)
train.data <- readRDS("./traindata.rds")

X <- as.matrix(twofour[,2:24])
b <- svd(X)
plot(b$u[,1],b$u[,2])
plot(X%*%b$v[,1],X%*%b$v[,2])
# perfrom pca on 24 genes
twofour.pca <- prcomp(as.matrix(twofour[,2:24]))
pairs(twofour.pca$x[,1:5])

# perform pca on 24 genes for training data
twofour.train <- inner_join(train.data,twofour)
twofour.train.pca <- prcomp(as.matrix(twofour.train[,4195:4217]))
cols <- character(nrow(twofour.train))
cols[] <- "black"
cols[twofour.train$PTGENDER == 'Male'] <- "blue"
pairs(twofour.train.pca$x[,1:5],col=cols)
