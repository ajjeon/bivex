getwd()
file <- readLines("output/test.txt")
file <- gsub("[[:punct:]]", "", file)
file <- gsub("^ ", "", file)
tc <- textConnection(file)
processed <- read.table(tc, sep=" ", na.string="---")
close(tc)

df <- as.data.frame(matrix(as.vector(unlist(processed)),nrow=100))

image(matrix(as.vector(unlist(processed)),nrow=100),
      col = c("white","azure3","azure4","black"),
      xlab="nucleosomes",
      ylab="iterations",
      axes=FALSE)

axis(1, at=seq(0,1, length=11), labels=seq(0,100,10))
axis(2, at=seq(0,1, length=11), labels=seq(0,5000,500))
 

x <- read.table("output/test.txt")
