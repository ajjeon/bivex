library(RColorBrewer)
library(gplots)

getwd()
read.table("output/test.txt")

readLines("output/test.txt")

file <- readLines("output/test.txt")
file <- gsub("[[:punct:]]", "", file)
file <- gsub("^ ", "", file)
tc <- textConnection(file)
processed <- read.table(tc, sep=" ", na.string="---")
close(tc)

df <- as.data.frame(t(matrix(as.vector(unlist(processed)),nrow=100)))

image(matrix(as.vector(unlist(processed)),nrow=100),col = c("white","azure3","azure4","black"))


ggplot(data = df, aes(x = metrics, y = teams)) +
  geom_tile(aes(fill = performance)) 
