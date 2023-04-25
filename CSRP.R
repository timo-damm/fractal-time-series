#=================================
# COMPLEX SCIENCE RESEARCH PROJECT
#=================================
#uses the a temporal network dataset on ocntacts in a workplace
#found at http://www.sociopatterns.org/datasets/test/
#---- general preparation ----

library(dplyr)
library(sna)
library(foreign)
library(nonlinearTseries)
library(ggplot2)
library(pracma)

#----data import, inspection and cleaning-----
data <- read.table("/home/dot/PAOS_ResMa/complex_systems_track/CSRP/data/tij_InVS15.dat")
metadata <- read.table("/home/dot/PAOS_ResMa/complex_systems_track/CSRP/data/metadata_InVS15.txt")

colnames(data)<- c("t","i","j") #where t = time, i and j are IDs of the people in contact
sum(is.na(data)) #no missing data
range(data$t) #range of time
range(data$i) #range of indiviuals
range(data$j) #range of individuals

colnames(metadata) <- c("i","di") #where i = ID of individual and d = department

#combining the interaction information with information on the departments
data <- merge(data,metadata, by = "i")
colnames(metadata) <- c("j","dj")
data <- merge(data,metadata, by = "j")
list_of_departments <- unique(data$di)

#removing out of office hours from the data
unique_t <- sort(as.data.frame(unique(data$t)))
unique_t$n <- seq.int(nrow(unique_t))
colnames(unique_t) <- c("t","ID")
data<-merge(data,unique_t, by = "t")

#create data on interaction frequencies (whole organisation)
frequency_whole <- data %>%
  count(ID) %>%
  sort(ID, decreasing = F)

#create data on interaction frequencies (within departments)
data_SFLE <- subset(data, di == "SFLE" & dj == "SFLE")

data_DMI <- subset(data, di == "DMI" & dj == "DMI")

data_SCOM <- subset(data, di == "SCOM" & dj == "SCOM")

data_SRH <- subset(data, di == "SRH" & dj == "SRH")

data_DCAR <- subset(data, di == "DCAR" & dj == "DCAR")

data_DST <- subset(data, di == "DST" & dj == "DST")

data_DSE <- subset(data, di == "DSE" & dj == "DSE")

data_DMCT <- subset(data, di == "DMCT" & dj == "DMCT")

data_DISQ <- subset(data, di == "DISQ" & dj == "DISQ")

data_SSI <- subset(data, di == "SSI" & dj == "SSI")

data_SDOC <- subset(data, di == "SDOC" & dj == "SDOC")

data_DG <- subset(data, di == "DG" & dj == "DG")


frequency_SFLE <- data_SFLE %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_DMI <- data_DMI %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_SCOM <- data_SCOM %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_SRH <- data_SRH %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_DCAR <- data_DCAR %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_DST <- data_DST %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_DSE <- data_DSE %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_DMCT <- data_DMCT %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_DISQ <- data_DISQ %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_SSI <- data_SSI %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_SDOC <- data_SDOC %>%
  count(ID) %>%
  sort(ID, decreasing = F)

frequency_DG <- data_DG %>%
  count(ID) %>%
  sort(ID, decreasing = F)


#----network analysis InVS-----
#preparing the network data
data_net <- subset(data, select = c(i,j)) #selecting only the individual IDs
data_net <- data_net %>% distinct() #removing duplicate rows, to avoid multiplex network
data_net <- as.network(data_net)

#descriptive statistics
network.size(data_net) #number of nodes in the network

network.density(data_net) #density of the network

sum(degree(data_net))/network.size(data_net) #average degree in the network

max(degree(data_net)) #maximum degree in the network

data_graph <- intergraph::asIgraph(data_net)
igraph::transitivity(data_graph) #global clustering coefficient in the network

max(igraph::clusters(data_graph, mode = "weak")$csize) #size of the largest component
#there is only one component, so the percentage is inevitably 100%

#visualisation
gplot(data_net, usearrows = FALSE, displaylabels = FALSE, 
      displayisolates = F, vertex.col = "darkgray", vertex.cex = 1.25, 
      edge.col = "black") +
  title(main = "InVS contact network")

#----network analysis DSE-----
#preparing the network data
DSE_net <- subset(data_DSE, select = c(i,j)) #selecting only the individual IDs
DSE_net <- DSE_net %>% distinct() #removing duplicate rows, to avoid multiplex network
DSE_net <- as.network(DSE_net)

#descriptive statistics
network.size(DSE_net) #number of nodes in the network

network.density(DSE_net) #density of the network

sum(degree(DSE_net))/network.size(DSE_net) #average degree in the network

max(degree(DSE_net)) #maximum degree in the network

DSE_graph <- intergraph::asIgraph(DSE_net)
igraph::transitivity(DSE_graph) #global clustering coefficient in the network

max(igraph::clusters(DSE_graph, mode = "weak")$csize) #size of the largest component
#there is only one component, so the percentage is inevitably 100%

#visualisation
gplot(DSE_net, usearrows = FALSE, displaylabels = FALSE, 
      displayisolates = F, vertex.col = "darkgray", vertex.cex = 1.25, 
      edge.col = "black") +
  title(main = "DSE department contact network")

# ------detrended fluctuation analysis (DFA)------
#DFA for whole organisation
#timeseries_data <- subset(data, select = -c(di,dj,t))

dfa_whole <- dfa(time.series = frequency_whole$n, window.size.range = c(30,18488), npoints = 200)
hurst_whole <- hurstexp(frequency_whole$n, d = 30)

#DFA intra-department
#DMI department
dfa_DMI <- dfa(time.series = frequency_DMI$n, window.size.range = c(30,18488), npoints = 200)
hurst_DMI <- hurstexp(frequency_DMI$n, d = 30)

#DISQ department
dfa_DISQ <- dfa(time.series = frequency_DISQ$n, window.size.range = c(30,18488), npoints = 200)
hurst_DISQ <- hurstexp(frequency_DISQ$n, d = 30)

#DSE department
dfa_DSE <- dfa(time.series = frequency_DSE$n, window.size.range = c(30,18488), npoints = 200)
hurst_DSE <- hurstexp(frequency_DSE$n, d = 30)

#----visualisation----
#whole organisation
ggplot(data=frequency_whole, aes(x=ID, y=n)) +
  geom_line()

#SFLE department
ggplot(data= frequency_SFLE, aes(x=ID, y=n)) +
  geom_line()

#DMI department
ggplot(data= frequency_DMI, aes(x=ID, y=n)) +
  geom_line()

#SCOM department
ggplot(data= frequency_SCOM, aes(x=ID, y=n)) +
  geom_line()

#DST department
ggplot(data= frequency_DST, aes(x=ID, y=n)) +
  geom_line()

#SRH dpeartment
ggplot(data= frequency_SRH, aes(x=ID, y=n)) +
  geom_line()
