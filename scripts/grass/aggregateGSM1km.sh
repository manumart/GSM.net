


for wd in 0 5 15
do
for var in lower mean upper 
do
	r.import --overwrite \
	  input=/home/manuel/Y/communs_infosol/Projets/GlobalSoilMap.net/GSM_90m/SOC/soc_cor_trans/$var\
	  '/soc_cor_l'$wd'_'$var'.tif' output=soc_cor_l$wd'_'$var

done
done

r.mapcalc --overwrite  expression='logsoc_cor_l0_mean=log(soc_cor_l0_mean)'

for wd in 0 5 15
do
for var in lower mean upper 
do
  exp=logsoc_cor_l$wd'_'$var'=log(soc_cor_l'$wd'_'$var')'
  echo $exp
  r.mapcalc --overwrite  expression="$exp"

done
done

## variances
for wd in 0 5 15
do
  exp=varsoc_cor_l$wd' = ((logsoc_cor_l'$wd'_mean - logsoc_cor_l'$wd'_lower) / 1.6448)^2'
  echo $exp
  r.mapcalc --overwrite  expression="$exp"

done

#averaged C content in the top 3 layers
r.mapcalc --overwrite  expression='meanCContent=(soc_cor_l0_mean*5/30)+(soc_cor_l5_mean*10/30)+(soc_cor_l15_mean*15/30)'

## this seems to work better.
## we now can see that distributions are left vs. right skewed depending on the region 
r.mapcalc --overwrite  expression="diff2=(soc_cor_l15_upper-soc_cor_l15_lower)/2 - soc_cor_l15_mean"
## Thus have a conservative approach : the biggest difference with the mean
r.mapcalc --overwrite  expression="maxDiff15=max(soc_cor_l15_upper-soc_cor_l15_mean) - max(soc_cor_l15_mean-soc_cor_l15_lower)"

## compute the standard deviation as q95% / 1.6448
## max(mean - 0.05%pc, 0.95%pc - mean) / 1.6448
for wd in 0 5 15
do
  exp=sd$wd'=max(soc_cor_l'$wd'_upper-soc_cor_l'$wd'_mean,soc_cor_l'$wd'_mean-soc_cor_l'$wd'_lower)/1.6448'
  echo $exp
  r.mapcalc --overwrite  expression=$exp
done

## variance des profils : teneurs calculées au prorata du mélange 5/30, 10/30, 15/30
r.mapcalc --overwrite  expression='varAll=(sd0*5/30)^2+(sd5*10/30)^2+(sd15*15/30)^2'

## creates a region with 1km resolution
g.mapset location=GSM90L93 mapset=1km

r.resamp.stats --o input=varAll@PERMANENT output=sumVarAll method=sum
r.resamp.stats --o input=varAll@PERMANENT output=meanVarAll method=average
r.mapcalc --overwrite  expression="m=sumVarAll/meanVarAll"
r.mapcalc --overwrite  expression="var1km=sumVarAll/m^2"
r.colors -g map=var1km@1km color=bgyr
r.resamp.stats --o input=meanCContent@PERMANENT output=meanCContent30cm method=average
r.mapcalc --overwrite  expression="meanCContent30cmLower=max(0,meanCContent30cm-1.6448*sqrt(var1km))"
r.mapcalc --overwrite  expression="meanCContent30cmUpper=meanCContent30cm+1.6448*sqrt(var1km)"
r.mapcalc --overwrite  expression="SOCstock30cm=meanCContent30cm*0.98*0.3*(1.66-0.318*(meanCContent30cm/10)^.5)*1000/1000"
r.mapcalc --overwrite  expression="SOCstock30cmUpper=meanCContent30cmUpper*0.98*0.3*(1.66-0.318*(meanCContent30cmUpper/10)^.5)*1000/1000"
r.mapcalc --overwrite  expression="SOCstock30cmLower=meanCContent30cmLower*0.98*0.3*(1.66-0.318*(meanCContent30cmLower/10)^.5)*1000/1000"



## following Titia's method for computing stocks :
## 20% rock fragments
## 
r.mapcalc --overwrite  expression="stock1km=meanCContent30cm*0.3*0.8*1.2"

g.mapset location=GSP mapset=carbon30cm

### then some R code (aggregateSOCS.r) useless : the writeRAST func does not work well
#RScript "../../R/aggregateSOCS.r"
## Should continue within R with the raster package...
r.proj --o location=GSM90L93 mapset=1km input=SOCstock30cm                          
r.proj --o location=GSM90L93 mapset=1km input=SOCstock30cmUpper
r.proj --o location=GSM90L93 mapset=1km input=SOCstock30cmLower
r.proj --o location=GSM90L93 mapset=1km input=meanCContent30cm
r.proj --o location=GSM90L93 mapset=1km input=meanCContent30cmLower
r.proj --o location=GSM90L93 mapset=1km input=meanCContent30cmUpper
r.mapcalc --overwrite  expression="test1=SOCstock30cmUpper - SOCstock30cm"
r.mapcalc --overwrite  expression="test2=SOCstock30cmLower - SOCstock30cm"


d.mon wx2
d.rast map=meanCContent30cm
d.rast map=meanCContent30cmLower
d.rast meanCContent30cmUpper
d.rast SOCstock30cm
d.rast SOCstock30cmUpper
d.rast SOCstock30cmLower
d.rast test1
d.rast test1
d.rast test2

## export for GSP


for r in SOCstock30cm SOCstock30cmLower SOCstock30cmUpper meanCContent30cm meanCContent30cmLower meanCContent30cmUpper
do
  echo 
#carbon 30 cm
#  r.out.gdal --o input=$r@carbon30cm output=/home/manuel/dev/GSM.net/data/GSPExport/$r.tiff format=GTiff
#  r.out.gdal --o input=$r@carbon30cm output=/home/manuel/dev/GSM.net/data/GSPExport/$r.png format=PNG
	r.out.gdal --o input=$r@carbon30cm output=/home/manuel/dev/GSM.net/data/GSPExport/$r.png format=PNG
done


### other maps, and among others covariates
## LandUse
for var in bdforet ecoclim 
do
    r.import --overwrite \
        input=/home/manuel/Y/communs_infosol/Projets/GlobalSoilMap.net/GSM500m/Covariates_500m/$var'_500m.tif'\
	   output=$var
done


g.mapset location=GSM90L93 mapset=500m
for r in bdforet ecoclim
do
  echo 
#carbon 30 cm
  r.out.gdal --o input=$r@500m output=/home/manuel/dev/GSM.net/data/GSPExport/$r.tiff format=GTiff
#  r.out.gdal --o input=$r@carbon30cm output=/home/manuel/dev/GSM.net/data/GSPExport/$r.png format=PNG
done


