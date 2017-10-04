
## simplest way :
setwd(homeS)
require(rgrass7)

initGRASS(gisBase = "/usr/local/grass-7.1.svn/", 
	   home = tempdir(), 
           gisDbase = "/home/manuel/dev/grass7/",
           location = "GSM90L93", mapset = "1km",
           override = TRUE)

SOC <- readRAST("meanCContent")
SOCL <- readRAST("meanCContent30cmLower")@data
SOCU <- readRAST("meanCContent30cmUpper")@data

varSOC <- readRAST("var1km")

SOCDf <- SOC@data
SOCDf$var <- varSOC@data

SOCDf$SOCstock30cm <- SOCDf$meanCContent*0.98*0.3*(1.66-0.318*(SOCDf$meanCContent/10)^.5)*1000/1000

SOCDf$SOCLowerStock30cm <- SOCL$meanCContent30cmLower*0.98*0.3*(1.66-0.318*(SOCL$meanCContent30cmLower/10)^.5)*1000/1000
SOCDf$SOCUpperStock30cm <- SOCU$meanCContent30cmUpper*0.98*0.3*(1.66-0.318*(SOCU$meanCContent30cmUpper/10)^.5)*1000/1000

SOC <- readRAST("meanCContent")
SOCL <- readRAST("meanCContent30cmLower")@data
SOCU <- readRAST("meanCContent30cmUpper")@data

varSOC <- readRAST("var1km")

SOC@data <- SOCDf

summary((1.66-0.318*(SOCDf$meanCContent/10)^.5)*1000)
sum(SOCDf$stock, na.rm = T)*1e6

SOC@data <- SOCDf
writeRAST(x=SOC, vname="SOCstock30cm", overwrite = T)
writeRAST(x=SOC, vname="SOCLowerStock30cm", overwrite = T)
writeRAST(SOC, vname="SOCUpperStock30cm", overwrite = T)

a <- readRAST("SOCstock30cm")@data
b <- readRAST("SOCstock30cmLower")@data
c <- readRAST("SOCstock30cmUpper")@data

SOC <- readRAST("SOCstock30cm")
SOC$SOCstock30cmLower <- b$SOCstock30cmLower
SOC$SOCstock30cmUpper <- c$SOCstock30cmUpper
sum(a$SOCstock30cm, na.rm = T)*1e6
sum(b$SOCstock30cmLower, na.rm = T)*1e6
sum(c$SOCstock30cmUpper, na.rm = T)*1e6
toto <- as.data.frame(SOC)
require(lattice)
require(tydir)
totoT <- gather(toto, var, stockC, -s1, -s2)

pdf("GSPStocks.pdf")
trellis.par.set(regions=list(col=topo.colors(100)))
print(levelplot(stockC~s1+s2|var, data = totoT))
dev.off()

layout(matrix(1:4, 1, 4), widths = c(4, 4 ,4,1))
plot(SOC["SOCstock30cm"],what = "image", zlim = c(0,1))


SOC@stock <- NA
plot(SOC)


require(raster)
SOC <- raster()
plot(rVarAll)
