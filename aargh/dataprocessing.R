##---
##title: "Data processing for nationam BAM analyses"
##author: "Peter Solymos"
##date: "Apr 18, 2016"
##output:
##  pdf_document:
##    toc: true
##    toc_depth: 2
##---

### Preliminaries start here

## Define root folder where data are stored
ROOT <- "c:/bam/May2015"
ROOT2 <- "e:/peter/bam/Apr2016"

## Load required packages
library(mefa4)
library(RODBC)
library(maptools)
library(QPAD)

## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

### Pulling in tables

## Define MS Access database connection
con <- odbcConnectAccess2007(file.path(ROOT, "BAM_BayneAccess_BAMBBScore.accdb"))
con2 <- odbcConnectAccess2007(file.path(ROOT2, "BBS_V4draft_2016_aprilForErin.accdb"))

#### Species lookup table

TAX <- sqlFetch(con, "dbo_DD_BAM_AOU52_SPECIES_TBL")
TAX$SSMA_TimeStamp <- NULL

#### SS level stuff for BAM+BBS combined

## This table has time zone, BCR, jurisdiction, XY
SS01 <- sqlFetch(con, "dbo_BBSBAM_V4_XYTBL_ATTRIBUTES1")
SS01 <- nonDuplicated(SS01, SS, TRUE)
SS01$COUNTRY <- ifelse(SS01$JURSALPHA %in% c("AB","BC","MB","NB",
    "NL","NS","NT","NU","ON","PEI","QC","SK","YK"), "CAN", "USA")
save(SS01, file=file.path(ROOT2, "out", "SS-regions-and-xy.Rdata"))

## Tree proportions
SS02 <- sqlFetch(con, "dbo_TREE_BBSBAM_V4_tbl")
SS02 <- nonDuplicated(SS02, SS, TRUE)
SS02 <- SS02[rownames(SS01),]

## TREE update 2016-12-01
tr <- read.csv("e:/peter/bam/Apr2016/tree_update/BAM_BBS_2015_XY_AVHRRTREE_nov30_2016.csv")
SS02$TREE_WRONG <- SS02$TREE
SS02$TREE <- tr$VCFAVHRR[match(SS02$SS, tr$SS)]

SS02$TREE[SS02$TREE > 100] <- NA
SS02$TREE[SS02$TREE < 0] <- NA
SS02$TREE <- SS02$TREE / 100
SS02$TREE3 <- factor(NA, levels=c("Open", "Sparse", "Dense"))
SS02$TREE3[SS02$TREE < 0.25] <- "Open"
SS02$TREE3[SS02$TREE >= 0.25 & SS02$TREE < 0.60] <- "Sparse"
SS02$TREE3[SS02$TREE >= 0.60] <- "Dense"

## Point level land cover
SS03 <- sqlFetch(con, "dbo_BAMBBS_LANDCOVER_PTS")
SS03 <- nonDuplicated(SS03, SS, TRUE)
SS03 <- SS03[rownames(SS01),]

## comment out all LCC05 and EOSD
## Reclass LCC05
ltlcc <- read.csv("~/repos/bamanalytics/lookup/lcc05.csv")
SS03$LCC05 <- SS03$LCC05_PT
SS03$LCC05_PT[SS03$LCC05_PT < 1 | SS03$LCC05_PT > 39] <- NA
SS03$LCC05_PT[SS01$COUNTRY == "USA"] <- NA
SS03$HAB_LCC1 <- ltlcc$BAMLCC05V2_label1[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC2 <- ltlcc$BAMLCC05V2_label2[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)]
SS03$HAB_LCC1 <- relevel(SS03$HAB_LCC1, "ConifDense")
SS03$HAB_LCC2 <- relevel(SS03$HAB_LCC2, "Conif")

if (FALSE) {
## for Nicole
SS_Nicole <- data.frame(
    PCODE=SS01$PCODE,
    SS=SS01$SS,
    X=SS01$X_GEONAD83,
    Y=SS01$Y_GEONAD83,
    JURS=SS01$JURSALPHA,
    COUNTRY=SS01$COUNTRY,
    BCR=as.factor(SS01$BCR),
    SS03[,c("LCC05","HAB_LCC1","HAB_LCC2")])

## Reclass EOSD
lteosd <- read.csv("~/repos/bamanalytics/lookup/eosd.csv")
levels(SS03$EOSD_PT) <- sub(",", "", levels(SS03$EOSD_PT))
SS03$EOSD_PT <- as.integer(as.character(SS03$EOSD_PT))
SS03$EOSD_PT[SS03$EOSD_PT < 1] <- NA
SS03$HAB_EOSD1 <- lteosd$Reclass_label1[match(SS03$EOSD_PT, lteosd$Value)]
SS03$HAB_EOSD2 <- lteosd$Reclass_label2[match(SS03$EOSD_PT, lteosd$Value)]
SS03$HAB_EOSD1 <- relevel(SS03$HAB_EOSD1, "ConifDense")
SS03$HAB_EOSD2 <- relevel(SS03$HAB_EOSD2, "Conif")
}

## Reclass NALCMS
ltnalc <- read.csv("~/repos/bamanalytics/lookup/nalcms.csv")
SS03$HAB_NALC2 <- ltnalc$Label[match(SS03$NALCMS_PT, ltnalc$Value)]
tmp <- as.character(interaction(SS03$HAB_NALC2, SS02$TREE3, sep="", drop=TRUE))
SS03$HAB_NALC1 <- as.character(SS03$HAB_NALC2)
ii <- SS03$HAB_NALC1 %in% c("Conif", "Decid", "Mixed", "Wet")
SS03$HAB_NALC1[ii] <- tmp[ii]
SS03$HAB_NALC1 <- as.factor(SS03$HAB_NALC1)
SS03$HAB_NALC1 <- relevel(SS03$HAB_NALC1, "ConifDense")
SS03$HAB_NALC2 <- relevel(SS03$HAB_NALC2, "Conif")

## NALC is used for QPADv3
if (FALSE) {
## LCC for offsets
SS03$LCC_OFF1 <- as.factor(ltlcc$qpad_num[match(SS03$LCC05_PT, ltlcc$lcc05v1_2)])
SS03$LCC_OFF2 <- factor(5, 1:5)
SS03$LCC_OFF2[is.na(SS03$HAB_NALC1)] <- NA
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("DecidSparse")] <- "4"
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("ConifSparse","MixedSparse")] <- "3"
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("DecidDense")] <- "2"
SS03$LCC_OFF2[SS03$HAB_NALC1 %in% c("ConifDense","MixedDense")] <- "1"
SS03$LCC_combo <- SS03$LCC_OFF1
SS03$LCC_combo[is.na(SS03$LCC_OFF1)] <- SS03$LCC_OFF2[is.na(SS03$LCC_OFF1)]
}

## Grid ID 4x4 km
SS_grid <- read.csv(file.path(ROOT, "BAMBBS_Gridcode.csv"))
rownames(SS_grid) <- SS_grid$SS
compare_sets(rownames(SS01), rownames(SS_grid))
SS_grid <- SS_grid[rownames(SS01),"gridcode",drop=FALSE]
levels(SS_grid$gridcode) <- gsub(",", "", levels(SS_grid$gridcode))

## Road: dist to, class, #lanes, surface
SS_road <- sqlFetch(con, "dbo_BAMBBS_2015_NearDistanceRoadJoin1000M")
rownames(SS_road) <- SS_road$SS
compare_sets(rownames(SS01), rownames(SS_road))
SS_road <- SS_road[rownames(SS01),]
SS_road$d2road <- SS_road[["Distance to Road"]]
table(SS_road$ROADCLASS, SS_road$d2road > 0, useNA="a")
table(SS_road$NBRLANES, SS_road$d2road > 0, useNA="a")
table(SS_road$PAVSTATUS, SS_road$d2road > 0, useNA="a")
table(SS_road$ROADCLASS, SS_road$NBRLANES, useNA="a")
SS_road <- SS_road[,c("d2road","ROADCLASS","NBRLANES","PAVSTATUS")]
## need to exclude # lanes >2

## Fire: year, size
SS_fire <- sqlFetch(con, "dbo_BBSBAM_2015_FIRE")
rownames(SS_fire) <- SS_fire$SS
compare_sets(rownames(SS01), rownames(SS_fire))
SS_fire <- SS_fire[rownames(SS01),]
SS_fire <- SS_fire[,c("Year","SIZE_HA")]
colnames(SS_fire) <- c("YearFire","FIRE_HA")

## Terrain: slope, twi, elev
SS_terr <- sqlFetch(con, "dbo_BBSBAM_2015_TERRAIN90")
rownames(SS_terr) <- SS_terr$SS
compare_sets(rownames(SS01), rownames(SS_terr))
SS_terr <- SS_terr[rownames(SS01),]
t(table(is.na(SS_terr$cti90), SS01$PCODE)) # mostly affects BBS in some states
SS_terr <- SS_terr[,c("slp90","cti90","elv90")]

## Climate variables from DS and NALCMS 4x4 level
SS_clim <- sqlFetch(con, "dbo_BBSBAM_2015__CLIMLU")
rownames(SS_clim) <- SS_clim$SS
compare_sets(rownames(SS01), rownames(SS_clim))
SS_clim <- SS_clim[rownames(SS01),]
tmp <- as.matrix(SS_clim[,grepl("NALCMS05_", colnames(SS_clim))])
SS_clim <- SS_clim[,!grepl("NALCMS05_", colnames(SS_clim))]
colnames(tmp) <- gsub("NALCMS05_", "", colnames(tmp))
Col <- as.character(ltnalc$Label)[match(colnames(tmp), as.character(ltnalc$Value))]
Col[is.na(Col)] <- "Water"
## 4x4 stuff is not done
#SS_NALC4x4 <- data.frame(groupSums(tmp, 2, Col, na.rm=TRUE))
#colnames(SS_NALC4x4) <- paste0("GRID4_NALC_", colnames(SS_NALC4x4))
SS_clim$NALCMS05 <- NULL
SS_clim$PCODE <- NULL
SS_clim$SS <- NULL

## 4x4 stuff is not done
if (FALSE) {
## LCC05 4x4 level
SS_LCC4x4 <- sqlFetch(con, "dbo_BBSBAM_V4_LCC05CND_4X4SUMM")
SS_LCC4x4 <- nonDuplicated(SS_LCC4x4, SS, TRUE)
#rownames(SS_LCC4x4) <- SS_LCC4x4$SS
compare_sets(rownames(SS01), rownames(SS_LCC4x4))
SS_LCC4x4 <- SS_LCC4x4[rownames(SS01),]
SS_LCC4x4$SS <- NULL
SS_LCC4x4$gridcode <- NULL
SS_LCC4x4$LCCVVSUM <- NULL
SS_LCC4x4 <- as.matrix(SS_LCC4x4)
colnames(SS_LCC4x4) <- gsub("LCCVV", "", colnames(SS_LCC4x4))
#Col <- as.character(ltlcc$BAMLCC05V2_label1)[match(colnames(SS_LCC4x4),
#    as.character(ltlcc$lcc05v1_2))]
Col <- as.character(ltlcc$BAMLCC05V2_label2)[match(colnames(SS_LCC4x4),
    as.character(ltlcc$lcc05v1_2))]
Col[is.na(Col)] <- "BARREN"
SS_LCC4x4 <- data.frame(groupSums(SS_LCC4x4, 2, Col, na.rm=TRUE))
SS_LCC4x4[is.na(SS03$HAB_LCC2),] <- NA
colnames(SS_LCC4x4) <- paste0("GRID4_LCC_", colnames(SS_LCC4x4))

## EOSD 4x4 level
SS_EOSD4x4 <- sqlFetch(con, "dbo_BBSBAM_V4_EOSD_4X4SUMM")
SS_EOSD4x4$upsize_ts <- NULL
rownames(SS_EOSD4x4) <- SS_EOSD4x4$SS
compare_sets(rownames(SS01), rownames(SS_EOSD4x4))
SS_EOSD4x4 <- SS_EOSD4x4[match(rownames(SS01), rownames(SS_EOSD4x4)),]
rownames(SS_EOSD4x4) <- SS01$SS

SS_EOSD4x4 <- as.matrix(SS_EOSD4x4[,grepl("eosdVV", colnames(SS_EOSD4x4))])
colnames(SS_EOSD4x4) <- gsub("eosdVV", "", colnames(SS_EOSD4x4))

#Col <- as.character(lteosd$Reclass_label1)[match(colnames(SS_EOSD4x4),
#    as.character(lteosd$Value))]
Col <- as.character(lteosd$Reclass_label2)[match(colnames(SS_EOSD4x4),
    as.character(lteosd$Value))]
Col[is.na(Col)] <- "BARREN"
SS_EOSD4x4 <- data.frame(groupSums(SS_EOSD4x4, 2, Col, na.rm=TRUE))
SS_EOSD4x4[is.na(SS03$HAB_EOSD2),] <- NA
colnames(SS_EOSD4x4) <- paste0("GRID4_EOSD_", colnames(SS_EOSD4x4))
}

## HEIGHT (Simard)
SS_height <- sqlFetch(con, "dbo_Height")
SS_height <- nonDuplicated(SS_height, SS, TRUE)
compare_sets(rownames(SS01), rownames(SS_height))
SS_height <- SS_height[rownames(SS01),]
SS_height <- SS_height[,"HEIGHTSIMARD",drop=FALSE]

if (FALSE) {
## Nature Serve range: 3 spp (Can clipped range used 0/1)
SS_nserv <- sqlFetch(con, "dbo_BBSBAM_SARSPPLOCATIONSRange")
SS_nserv <- nonDuplicated(SS_nserv, SS, TRUE)
compare_sets(rownames(SS01), rownames(SS_nserv))
SS_nserv <- SS_nserv[rownames(SS01),]
SS_nserv <- SS_nserv[,c("CAWAINOUT","OSFLINOUT","CONIINOUT")]
}

## GFW yearly loss intersections and 1st year of loss
#SS_gfw <- sqlFetch(con, "dbo_BAMBBS_GFWLossYear")
SS_gfw <- read.csv(file.path(ROOT, "GFWLossYear.csv"))
SS_gfw <- nonDuplicated(SS_gfw, SS, TRUE)
compare_sets(rownames(SS01), rownames(SS_gfw))
levels(SS_gfw$YearLoss) <- gsub(",", "", levels(SS_gfw$YearLoss))
SS_gfw$YearLoss <- as.integer(as.character(SS_gfw$YearLoss))
SS_gfw <- SS_gfw[rownames(SS01),]
SS_gfw <- SS_gfw[,"YearLoss",drop=FALSE]

## Pasher disturbance
SS_pash <- read.csv(file.path(ROOT, "bambbs2015beadandpasher.csv"))
SS_pash <- nonDuplicated(SS_pash, SS, TRUE)
compare_sets(rownames(SS01), rownames(SS_pash))
SS_pash <- SS_pash[rownames(SS01),]
SS_pash <- SS_pash[,c("BEADTotalL","BEADtotPol")]

## Local spring

SS_sprng <- read.csv("e:/peter/bam/May2015/NRCAN_SG_001_BAMBBS2015_71_13.csv")
SS_sprng <- SS_sprng[,c("SS","RASTERVALU")]
SS_sprng$SPRNG <- SS_sprng$RASTERVALU
levels(SS_sprng$SPRNG) <- gsub(",", "", levels(SS_sprng$SPRNG))
SS_sprng$SPRNG <- as.numeric(as.character(SS_sprng$SPRNG))
SS_sprng$SPRNG[SS_sprng$SPRNG < 0] <- NA
rownames(SS_sprng) <- SS_sprng$SS
SS_sprng <- SS_sprng[rownames(SS01),]

## Put together the main SS level object
SS <- data.frame(
    PCODE=SS01$PCODE,
    SS=SS01$SS,
    X=SS01$X_GEONAD83,
    Y=SS01$Y_GEONAD83,
    Xcl=SS01$X_CLCC,
    Ycl=SS01$Y_CLCC,
    JURS=SS01$JURSALPHA,
    COUNTRY=SS01$COUNTRY,
    TZONE=SS01$TZONE_CODE,
    BOREALLOC=SS01$BOREALLOC,
    BCR=as.factor(SS01$BCR),
    TREE=SS02$TREE,
    TREE3=SS02$TREE3,
    SPRNG=SS_sprng$SPRNG,
    #LCC05_PT=SS03$LCC05_PT, # -- FOR NICOLE
    SS03[,c("LCC05","HAB_LCC1","HAB_LCC2")],
    SS03[,c("HAB_NALC2", "HAB_NALC1")],
    SS_grid,
    #SS_nserv,
    SS_road,
    SS_terr,
    SS_fire,
    SS_clim,
    SS_pash,
    SS_gfw,
    SS_height)
    #SS_NALC4x4,
    #SS_LCC4x4,
    #SS_EOSD4x4)

#### Project summary table

## This table needed local tweaks to be operable
#PCODE <- sqlFetch(con, "dbo_National_Proj_Summary_V4_2015")
#PCODE$SSMA_TimeStamp <- NULL
PCODE <- read.csv(file.path(ROOT,"proj.csv"))
levels(PCODE$Maxdist) <- tolower(levels(PCODE$Maxdist))
levels(PCODE$Maxdist)[levels(PCODE$Maxdist)=="unlimited"] <- "Inf"
PCODE$Maxdist[PCODE$Maxdist=="unknown"] <- NA
PCODE$Maxdist <- droplevels(PCODE$Maxdist)
PCODE$Maxdist <- as.numeric(as.character(PCODE$Maxdist))
PCODE$Maxdur <- pmin(PCODE$MaxDuration, 10)

#### Survey level fields

## Pull in PKEY tables
pkbam <- sqlFetch(con, "dbo_National_PKEY_V4_2015")
pkbam$SSMA_TimeStamp <- NULL
#pkbbs <- sqlFetch(con, "dbo_PKEY_BBS_V3_2015")
pkbbs <- sqlFetch(con2, "PKEY_BBS_V4_2016")

## What columns to retain
PKCOLS <- c("PKEY","SS","PCODE","METHOD","SITE","STN","ROUND",
    "YEAR","MONTH","DAY","HOUR","MIN","PART","MAXDUR","MAXDIS")

colnames(pkbam) <- toupper(colnames(pkbam))
pkbam$MAXDUR <- PCODE$Maxdur[match(pkbam$METHOD, PCODE$Method)]
pkbam$MAXDIS <- PCODE$Maxdist[match(pkbam$METHOD, PCODE$Method)]
pkbam$PART <- 1L
pkbam$MONTH <- pkbam$MM
pkbam$DAY <- pkbam$DD
pkbam$HOUR <- pkbam$HR
pkbam$YEAR <- pkbam$YYYY
levels(pkbam$ROUND) <- sub("[[:alpha:]]+$", "", levels(pkbam$ROUND))
pkbam$ROUND <- as.integer(as.character(pkbam$ROUND))
pkbam <- pkbam[,PKCOLS]

colnames(pkbbs) <- toupper(colnames(pkbbs))
pkbbs$MAXDUR <- 3
pkbbs$MAXDIS <- Inf
pkbbs$PART <- 2L
pkbbs$METHOD <- as.factor("BBS:9999")
pkbbs$SITE <- as.factor(pkbbs$SITE)
pkbbs$YEAR <- pkbbs$YYYY
pkbbs$MONTH <- pkbbs$MM
pkbbs$DAY <- pkbbs$DD
pkbbs$HOUR <- pkbbs$HR
pkbbs <- pkbbs[,PKCOLS]

PKEY <- rbind(pkbam, pkbbs)
#rm(pkbam, pkbbs)
#gc()

## Map `METHOD` field from project summary table onto `PKEY$METHOD`
## so that duration and distance method can be carried forward to
## point count table
levels(PCODE$Method)[levels(PCODE$Method) == "QCAtlas:118"] <- "QCATLAS:118"
compare_sets(PCODE$Method, PKEY$METHOD)
setdiff(PKEY$METHOD, PCODE$Method)
setdiff(PCODE$Method, PKEY$METHOD)
PKEY$DURMETH <- PCODE$DURMETH[match(PKEY$METHOD, PCODE$Method)]
PKEY$DISMETH <- PCODE$DISTMETH[match(PKEY$METHOD, PCODE$Method)]

## Identifying roadside surveys
PKEY$ROAD <- 0L
treat.as.bbs <- c("HOBBS","CF","MNBBA", levels(pkbbs$PCODE))
PKEY$ROAD[PKEY$PCODE %in% treat.as.bbs] <- 1L

#### Offset specific variables

## Date/time components
PKEY$MIN[is.na(PKEY$MIN)] <- 0 # min that critical, is not -- said Yoda
MM <- ifelse(PKEY$MONTH < 10, paste0("0", PKEY$MONTH), as.character(PKEY$MONTH))
HH <- ifelse(PKEY$HOUR < 10, paste0("0", PKEY$HOUR), as.character(PKEY$HOUR))
mm <- ifelse(PKEY$MIN < 10, paste0("0", PKEY$MIN), as.character(PKEY$MIN))
#mm[is.na(mm) & !is.na(HH)] <- "00"
DD <- with(PKEY, paste0(YEAR, "-", MM, "-", DAY, " ", HH, ":", mm, ":00"))
DD <- strptime(DD, "%Y-%m-%e %H:%M:%S")
PKEY$DATE <- DD
## Julian day
PKEY$JULIAN <- DD$yday # this is kept as original
PKEY$JDAY <- DD$yday / 365
summary(PKEY$JDAY)
## prevent too far extrapolation
PKEY$JDAY[PKEY$JDAY < 0.35 | PKEY$JDAY > 0.55] <- NA
## TSSR = time since sunrise
Coor <- as.matrix(cbind(as.numeric(SS$X),as.numeric(SS$Y)))[match(PKEY$SS, rownames(SS)),]
JL <- as.POSIXct(DD)
subset <- rowSums(is.na(Coor))==0 & !is.na(JL)
sr <- sunriset(Coor[subset,], JL[subset], direction="sunrise", POSIXct.out=FALSE) * 24
PKEY$srise <- NA
PKEY$srise[subset] <- sr
PKEY$start_time <- PKEY$HOUR + PKEY$MIN/60
TZ <- SS$TZONE[match(PKEY$SS, rownames(SS))]
lttz <- read.csv("~/repos/bamanalytics//lookup/tzone.csv")
lttz <- nonDuplicated(lttz, Timezone, TRUE)
PKEY$MDT_offset <- lttz$MDT_offset[match(TZ, rownames(lttz))]
table(TZ, PKEY$MDT_offset)
PKEY$TSSR <- (PKEY$start_time - PKEY$srise + PKEY$MDT_offset) / 24
PKEY$TSSR_orig <- PKEY$TSSR # keep a full copy
PKEY$TSSR[PKEY$start_time > 12] <- NA ## after noon
summary(PKEY$TSSR)
summary(PKEY$start_time)

PKEY <- PKEY[PKEY$DURMETH != "J",] # unknown duration
PKEY <- PKEY[PKEY$DISMETH != "O",] # unknown distance
PKEY <- droplevels(PKEY)

## QC Atlas problem
#with(PKEY, table(PCODE, is.na(MAXDUR)))
#with(PKEY, table(PCODE, is.na(MAXDIS)))


#### Point count tables and methodology

## Some of these tables are also tweaked locally (see below)
BEH <- sqlFetch(con, "dbo_DD_DescripBEH")
#DISINT <- sqlFetch(con, "dbo_DD_DescripDistance")
#DISINT$SSMA_TimeStamp <- NULL
#DURINT <- sqlFetch(con, "dbo_DD_DescripPeriod")
DISMET <- sqlFetch(con, "dbo_DD_distance_codes_methodology")
DURMET <- sqlFetch(con, "dbo_DD_duration_codes_methodology")
## Trailing whitespace removed from factor levels
levels(DISMET$DISTANCECODE) <- gsub(" *$", "", levels(DISMET$DISTANCECODE))
levels(DURMET$DURATIONCODE) <- gsub(" *$", "", levels(DURMET$DURATIONCODE))

#### Point count tables
pcbam <- sqlFetch(con, "dbo_National_PtCount_V4_2015")
pcbam$SSMA_TimeStamp <- NULL

pcbbs <- sqlFetch(con2, "POINTCOUNT_BBS_V4_2016")
colnames(pcbbs) <- toupper(colnames(pcbbs))
colnames(pcbbs)[colnames(pcbbs) == "SPECIES_ID"] <- "SPECIES"
colnames(pcbbs)[colnames(pcbbs) == "PERIOD"] <- "DURATION"
pcbbs$PCODE <- "BBS"

## Columns to keep
pccols <- c("PCODE","SS","PKEY","DURATION","DISTANCE",
    "SPECIES","ABUND","BEH")
## Close the database connections
close(con)
close(con2)

## Duration and distance intervals (locally tweaked)
DURINT <- read.csv("~/repos/bamanalytics/lookup/durint.csv")
DISINT <- read.csv("~/repos/bamanalytics/lookup/disint.csv")
DURINT$dur <- paste0(DURINT$Dur_Start, "-", DURINT$DUR_end)
DISINT$dis <- paste0(DISINT$DIST_START, "-", DISINT$DIST_END)
rownames(DURINT) <- DURINT[,1]
rownames(DISINT) <- DISINT[,1]
## Combined point count table
PCTBL <- rbind(pcbam[,pccols], pcbbs[,pccols])
## Mapping duration and distance intervals
PCTBL$dur <- as.factor(DURINT$dur[match(PCTBL$DURATION, rownames(DURINT))])
PCTBL$dis <- as.factor(DISINT$dis[match(PCTBL$DISTANCE, rownames(DISINT))])
## Methodology
PCTBL$DISMETH <- droplevels(PKEY$DISMETH[match(PCTBL$PKEY, PKEY$PKEY)])
PCTBL$DURMETH <- droplevels(PKEY$DURMETH[match(PCTBL$PKEY, PKEY$PKEY)])

## Filtering surveys (need to exclude PKEY)
keeppkey <- rep(TRUE, nrow(PCTBL))
## 11=0-20
## 8=unk
keeppkey[PCTBL$DURATION %in% c(11,8)] <- FALSE
## Excluding unknown distance bands
keeppkey[PCTBL$DISTANCE %in% c(4,5,9)] <- FALSE
## Excluding unknown duration methodology
keeppkey[PCTBL$DURMETH == "J"] <- FALSE
## Excluding unknown distance methodology
keeppkey[PCTBL$DISMETH == "O"] <- FALSE
## Actual filtering -- but dropping PKEYs
PCTBL <- droplevels(PCTBL[keeppkey,])

## Filtering within survey (do not exclude PKEY)
## Filtering behaviour
#sort(100 * table(PCTBL$BEH) / sum(table(PCTBL$BEH)))
## 1=Heard
## 11=no birds observed at station - added 2011
## 6=seen and heard
## Excluding non-aerial detections
table(PCTBL$BEH, PCTBL$PCODE=="BBS")
keep <- rep(TRUE, nrow(PCTBL))
keep[!(PCTBL$BEH %in% c("1","6","11"))] <- FALSE
## this is fake, but there is no other option until a fix
#keep[is.na(PCTBL$BEH)] <- TRUE # dont know what this in -- FIXED in BBS_V4

## Excluding >10 min intervals
## 10=10-20
## 3=before or after
## 9=10-15
keep[PCTBL$DURATION %in% c(10,3,9)] <- FALSE
## Excluding NA values
keep[is.na(PCTBL$dur)] <- FALSE
keep[is.na(PCTBL$dis)] <- FALSE
keep[is.na(PCTBL$ABUND)] <- FALSE
## Actual filtering -- but keeping PKEYs (do not drop levels)
#PCTBL$keep <- keep
PCTBL <- PCTBL[keep,]

## Excluding/dropping species

PCTBL$SPECIES <- droplevels(PCTBL$SPECIES)
levels(PCTBL$SPECIES) <- toupper(levels(PCTBL$SPECIES))
compare_sets(PCTBL$SPECIES, TAX$Species_ID)
setdiff(PCTBL$SPECIES, TAX$Species_ID)
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "YWAR"] <- "YEWA"
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "SCJU"] <- "DEJU" # change SCJU to DEJU
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "MYWA"] <- "YRWA" # change MYWA to YRWA
levels(TAX$Species_ID)[levels(TAX$Species_ID) == "COSN"] <- "WISN" # change COSN to WISN

levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "YWAR"] <- "YEWA"
levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "SCJU"] <- "DEJU" # change SCJU to DEJU
levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "MYWA"] <- "YRWA" # change MYWA to YRWA
levels(PCTBL$SPECIES)[levels(PCTBL$SPECIES) == "COSN"] <- "WISN" # change COSN to WISN
PCTBL$SPECIES <- droplevels(PCTBL$SPECIES)
setdiff(PCTBL$SPECIES, TAX$Species_ID)


PCTBL$SPECIES_ALL <- PCTBL$SPECIES
sspp <- read.csv("~/repos/bamanalytics/lookup/singing-species.csv")
levels(PCTBL$SPECIES)[!(levels(PCTBL$SPECIES) %in% sspp$Species_ID[sspp$Singing_birds])] <- "NONE"

## Excluding columns
PCTBL$DURATION <- NULL
PCTBL$DISTANCE <- NULL
PCTBL$BEH <- NULL
PCTBL$dur <- droplevels(PCTBL$dur)
PCTBL$dis <- droplevels(PCTBL$dis)

compare_sets(SS$SS, PKEY$SS)
compare_sets(SS$SS, PCTBL$SS)
compare_sets(PKEY$PKEY, PCTBL$PKEY)

save(SS, PKEY, PCTBL, TAX,
    file=file.path(ROOT2, "out",
    #paste0("data_package_2016-04-18.Rdata")))
#    paste0("data_package_2016-07-05.Rdata")))
    paste0("data_package_2016-12-01.Rdata")))
