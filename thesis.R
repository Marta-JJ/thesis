# Data analysis, MPhil thesis

# load packages
p <- c("ggplot2", "ggfortify", "purrr", "broom",
       "magrittr", "lattice", "tidyr", "lubridate", "dplyr")
lapply(p, require, character.only = TRUE)

# load all data
stb <- read.csv("/Users/marta/Desktop/Rscripts/Rthesis/stbgla.csv", header = T, na.strings = "*")
str(stb)
head(stb, n = 50)

# give Leaf column names Leaf 1-2 instead of 1-2
stb$Leaf <- factor(stb$Leaf)
levels(stb$Leaf) <- c("Leaf 1", "Leaf 2")

# change column name from Additive to QoI
colnames(stb)[9] <- "QoI"

# change avgGLA into numeric
stb$avgGLA <- as.numeric(as.character(stb$avgGLA))
View(stb)

# own theme
stb_theme <- theme_classic() +
  theme(text = element_text(color = "red"),
        panel.border = element_rect(linetype = "solid", fill = NA), # fill needs to be NA, otherwise it covers the panel
        axis.title = element_text(size = 13), 
        plot.title = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 10),
        panel.grid.major = element_line(color = "gray", size = 0.1),
        panel.grid.minor = element_line(color = "gray", size = 0.1),
        strip.background = element_rect(fill = "gray"),
        strip.text = element_text(size = 11, face = "bold"),
        legend.text = element_text(colour = "black"))

# histogram of values
stb %>% ggplot(aes(x = avgSTB)) +
  geom_histogram(binwidth = 15, alpha = 0.6) +
  facet_wrap(~Leaf) + 
  stb_theme + ylab("Count") + xlab("Average STB")

stb %>% ggplot(aes(x = avgGLA)) +
  geom_histogram(binwidth = 15, alpha = 0.6) +
  facet_wrap(~Leaf) + 
  stb_theme + ylab("Count") + xlab("Average GLA")

(sum <- stb %>% group_by(Leaf) %>% summarize(medianSTB = median(avgSTB, na.rm  = T),
                                             IQRSTB = IQR(avgSTB, na.rm = T)))


# yield data
yno <- read.csv("/Users/marta/Desktop/Rscripts/Rthesis/yieldnoout.csv", header = T, na.strings = "*")
str(yno)

# generic labs
yno_xlab <- "Imtrex (L/ha)"
yno_ylab <- "Mean yield (t/ha)"

# calculate means, sds and ses to add to plots
ynosum <- yno %>% group_by(Imtrex, Experiment) %>% 
  summarise(meanyield = mean(yield, na.rm = TRUE), sdyield = sd(yield, na.rm = TRUE), 
            seyield = sd(yield, na.rm=  TRUE)/sqrt(sum(n()))) %>%
  ggplot(ynosum, mapping = aes(x = Imtrex, y = meanyield, colour = Experiment, group = Experiment)) + 
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = meanyield - seyield, ymax = meanyield + seyield), width = 0.1) + 
  stb_theme + xlab(yno_xlab) + ylab(yno_ylab)

yno$Imtrex
# change labels
yno_lab <- c("0", "0.5", "1", "1.5", "2", "0", "0.5", "1", "1.5", "2", "0", "0.5", "1", "1.5", "2")

# filter out Oak Park and create a 3-way interaction boxplot
yno %>% filter(Site != "Oak Park") %>% 
  ggplot(aes(x = interaction(Imtrex, Additive), y = yield)) + 
  geom_boxplot(aes(colour = Additive), outlier.shape = 1, outlier.size = 2) + 
  facet_grid(Experiment~ .) + 
  scale_x_discrete(labels = yno_lab) + 
  xlab(yno_xlab) + ylab(yno_ylab) +
  guides(fill=guide_legend(title=NULL)) + stb_theme + 
  theme(legend.position = "top")

# load ratio qPCR data
qpcr <- read.csv("/Users/marta/Desktop/Rscripts/Rthesis/ratioqPCR.csv", header = T, na.strings = "*")
qpcr$Imtrex <- as.factor(qpcr$Imtrex)
qpcr$Trtmnt <- as.factor(qpcr$Trtmnt)
str(qpcr)


# calculate leaf 1 and leaf 2 pooled together
grandsums <- stb %>% group_by(Experiment, Trtmnt, Rep) %>% 
  summarise(meanSTB = mean(avgSTB, na.rm = TRUE))

grandsums1 <- stb %>% group_by(Experiment, Trtmnt, Rep) %>% 
  summarise(meanGLA = mean(avgGLA, na.rm = TRUE))

# merge both columns according to exp, trtmnt, rep use c() for a string of values
m1 <- merge(grandsums, grandsums1, by = c("Experiment", "Trtmnt", "Rep"))

# merge grandsums and create correlation plots for ratio and STB
merge(grandsums, qpcr, by = c("Experiment", "Trtmnt", "Rep")) %>%
ggplot(aes(x = meanSTB, y = RatFV, colour = Experiment)) + 
  geom_point(size = 2, shape = 1) + 
  geom_smooth(method = "lm") + 
  facet_wrap( ~ Experiment) + 
  xlab("Average STB (%)") + 
  ylab("ratio of resistant to sensitive") + 
  stb_theme + theme(legend.position = "none")


# filter data only for July
stb %>% filter(Month == "July") %>%
  ggplot(aes(x = interaction(Imtrex, QoI), y = avgGLA, colour = QoI)) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 2) + 
  facet_grid(Experiment ~ Leaf) + 
  stb_theme + theme(legend.position="top") + 
  scale_x_discrete(labels = c("0", "0.5", "1", "1.5", "2", "0", "0.5", "1", "1.5", "2", "0", "0.5", "1", "1.5", "2")) + 
  ylab("Average GLA (%)") + xlab("Imtrex (L/ha)") + scale_y_log10(breaks = 10 ^ (0:3)) + 
  guides(fill=guide_legend(title=NULL))

# filter data only for June
stb %>% filter(Month == "June") %>%
  ggplot(aes(x = interaction(Imtrex, QoI), y = avgGLA, colour = QoI)) + 
  geom_boxplot(outlier.shape = 1, outlier.size = 2) + 
  facet_grid(Experiment ~ Leaf) + 
  stb_theme + theme(legend.position="top") + 
  scale_x_discrete(labels = c("0", "0.5", "1", "1.5", "2", "0", "0.5", "1", "1.5", "2", "0", "0.5", "1", "1.5", "2")) +
  scale_y_log10() +
  ylab("Average GLA (%)") + xlab("Imtrex (L/ha)") + 
  guides(fill=guide_legend(title=NULL))


# glasshouse data set
glass <- read.csv("/Users/marta/Desktop/Rscripts/Rthesis/glass1.csv", header = T, na.strings = "*")

# tidy up plot
glass$pot <- as.factor(glass$pot)
glass$block <- as.factor(glass$block)
glass$Iso <- factor(glass$Iso, c("32","2.4.9","317", "325", "330", 
                                 "359", "413", "442", "459", "9"))
glass$mod <- factor(glass$mod, c("none","half"))
glass$imx <- factor(glass$imx, c("none", "one-eighth"))

# give columns new names
names(glass) <- c("pot", "block", "avdis", "Isolate", "Treatment", "Modem", "Imtrex")


# fit a simple model
lm <- lm(avdis ~ Isolate * Modem * Imtrex + block, data = glass)
autoplot(lm, smooth.colour = NA)
anova(lm)
summary(lm)

# check distribution of residuals
res <- resid(lm)
hist(res)

# remove insig. figures from the model
lm <- lm(avdis ~ Isolate * Modem * Imtrex + block - Isolate:Modem:Imtrex - 
           Modem:Imtrex - Isolate:Modem, data = glass)

install.packages("phia")
library(phia)

in.means <- interactionMeans(lm)
in.means
plot(in.means)
testInteractions(lm, fixed="Iso", across="imx") # imtrex has a significant effect on Iso317

# plot residuals from the model
res <- resid(lm)
hist(res)
qqnorm(res)

sumsim <- glass %>% group_by(Isolate, Imtrex) %>% 
  summarise(mean = mean(avdis, na.rm=TRUE), sd = sd(avdis, na.rm=TRUE), 
            se = sd(avdis, na.rm=TRUE)/sqrt(sum(n())))

# interaction plot for Isolate:Imtrex
ggplot(sumsim, aes(y = mean, x = Isolate, group = Imtrex, colour = Imtrex)) + geom_point(shape = 22) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) + 
  scale_x_discrete(labels = c("sen", "res", "317", "325", "330", 
                            "359", "413", "442", "459", "9")) + xlab("Isolate") + 
  ylab("Average disease (%)") + stb_theme + theme(legend.position="top")


microec50 <- read.csv("/Users/marta/Desktop/Rscripts/Rthesis/Microtiter/microec50.csv",
                      header = T, na.strings="*")

# change to factors
str(microec50)
microec50$iso <- as.factor(microec50$iso)
microec50$conc <- as.factor(microec50$conc)

levels(microec50$fung2) <- c("azoxystrobin", "pyraclostrobin")
colnames(microec50)[4] <- "QoI"

# histogram before data transformation 
ggplot(microec50, aes(x = ec50)) + 
  geom_histogram(binwidth = 0.5, alpha = 0.6) + 
  facet_wrap(~QoI, ncol = 1) + stb_theme + xlab("EC50") + ylab("Count")

## histogram after data transformation
microec50 %>% mutate(ec50, logec50 = log10(ec50)) %>%
ggplot(aes(x = logec50)) +
  geom_histogram(binwidth = 0.5, alpha = 0.6) +
  facet_wrap(~QoI, ncol = 1) + stb_theme + xlab("EC50") + ylab("Count")

# boxplots for each concentration
microec50 %>% mutate(ec50, logec50 = log10(ec50)) %>%
ggplot(aes(x = conc, y = logec50)) + 
  geom_boxplot(outlier.size = 1.5, outlier.shape = 1) + 
  facet_wrap(~med) + stb_theme + xlab("Concentration") + ylab("LogEC50")


microec50 %>% mutate(ec50, logec50 = log10(ec50)) %>% group_by(conc, QoI) %>% 
  summarise(mean = mean(logec50, na.rm = TRUE), 
            sd = sd(logec50, na.rm = TRUE), 
            se = sd(logec50, na.rm = TRUE)/sqrt(n())) %>%
ggplot(aes(x = conc, y = mean, group = QoI, colour = QoI)) + 
  geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width=0.1) + 
  xlab("Concentration of QoI") + ylab("Mean log10 of EC50 of fluxapyroxad") + 
  stb_theme + theme(legend.position = c(0.85, 0.9), 
                    legend.background = element_rect(color = "black"))


# remove sensitive strain
microec50 %>% filter(iso != "101") %>% group_by(conc, QoI) %>% mutate(ec50, logec50 = log10(ec50)) %>% 
  group_by(conc, QoI) %>% summarise(mean = mean(logec50, na.rm = TRUE), sd = sd(logec50, na.rm = TRUE), se = sd(logec50, na.rm = TRUE)/sqrt(n())) %>%
ggplot(aes(x = conc, y = mean, group = QoI, colour = QoI)) + 
  geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) + 
  xlab("Concentration of QoI") + ylab("Mean log10 of EC50 of fluxapyroxad") + stb_theme +
  theme(legend.position = c(0.35, 0.3), 
        legend.background = element_rect(color = "black"))

# interaction plot
microec50 %>% group_by(QoI, iso) %>% 
  summarise(meanec50 = mean(ec50, na.rm = TRUE), sd = sd(ec50, na.rm = TRUE), 
            se = sd(ec50, na.rm = TRUE)/sqrt(sum(n()))) %>%
ggplot(aes(x = QoI, y = meanec50, group = iso, colour = iso)) + 
  geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = meanec50 - se, ymax = meanec50 + se), width=0.1) + 
  xlab("QoI") + ylab("mean EC50 of fluxapyroxad") + 
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6)) + 
  stb_theme + theme(legend.position = "top")





