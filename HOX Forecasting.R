###############################################################
##  Pseudo out-of-sample forecasting: HOX Housing index      ##
###############################################################


library(stargazer)
library(AER)
library(lmtest)
library(tseries)
library(urca)
library(dynlm)
library(sandwich)
library(readxl)
library(forecast)
library(xts)


HOXSWE <- read_excel("Valueguardhela.xls")
str(HOXSWE)
class(HOXSWE) 


HOXSWE_dataframe <- HOXSWE[,"HOXSWE"]
View(HOXSWE_dataframe)
class(HOXSWE_dataframe)
HOXSWE_ts<- ts(HOXSWE_dataframe, start=c(2005,01), end = c(2016,01), frequency = 12)
HOXSWE <- ts(HOXSWE_ts, start=c(2005,01), end = c(2016,01), frequency = 12 )

ts.plot(HOXSWE)


HOXSWE_hela<- ts(HOXSWE_ts, start=c(2005,01), end = c(2018,01), frequency = 12)
plot(HOXSWE_hela)
is.ts(HOXSWE_hela)


class(HOXSWE) 
is.ts(HOXSWE) 

## Grafisk illustration av serien - Valueguard Boprisindexet
## från 2005-2016, d.v.s den mer begränsade tidsperioden.
plot(HOXSWE, main = "BOPRISINDEX HOXSWE", xlab = "Tid")
grid(col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
legend("topleft", "Boprisindex HOXSWE", lty=1, bg="lightgray",  cex = 0.75)

HOXSWE 
## Månadsdatan sträcker sig från januari månad 2005 (2005:M01) 
## ända fram till Januari år 2016 (2016:M01)


## Deskriptiv statistik i form av en tabell.
summary(HOXSWE_dataframe)
sum_HOXSWE <- summary(HOXSWE_dataframe)
stargazer(sum_HOXSWE,  type = 'latex', title = "Deskriptiv Statistik")

###############################################################
##                     ARIMA-modellering                     ##
##                                                           ##
###############################################################




################################################################
##
##  1) IDENTIFIERING 
##
################################################################


## En förutsättning för att genomföra en ARIMA-modellering med
## huvudsaklig grund i Box-Jenkin metodik inom ekonometri för att
## därmed genomföra prognostisering, är det först och främst av vikt
## att säkerställa att vår något begränsade serie är stationärt. Vi börjar 
## med att genomföra ett standardtest avseende enhetsrötter vid namn
## Augmented Dickey-Fuller test.


## Augmented Dickey-Fuller Test (ADF-test) - Test för enhetsrot.

ADF.test <- function(data){
  ADF <- function(type, data) {
    require(urca)
    result1 <- ur.df(data,
                     type = type,
                     lags = 3*frequency(data),
                     selectlags = "AIC")
    cbind(t(result1@teststat),result1@cval)
  }
  types <- c("trend", "drift", "none")
  result2 <- apply(t(types), 2, ADF, data)
  cat(rep("#", 17),'\n')
  cat("Augmented Dickey--Fuller test\n")
  cat(rep("#", 17),'\n')
  round(rbind(result2[[1]][c(1,3),],
              result2[[2]],
              result2[[3]]), 2)
}


## H0: Finns en enhetsrot
## HA: Finns ej en enhetsrot

ADF.test(HOXSWE) # URSPRUNGSDATA
## Då teststatistikan för tau3 i absoluta tal understiger
## de kritiska värdena på valfri signifikansnivå kan vi konstatera
## att vi har en enhetsrot. Vi kan i detta fall dra slutsatsen att vi
## kan förkasta determinstisk trend och drift samt att då nollhypotesen
## avseende om att en enhetsrot föreligger, inte kan förkastas!


## Efter ett kvantiferat konstaterande om att vår data inte är stationärt
## kan vi även genom autokorrelationsplotten (ACF) utröna ifall det visuellt
## är stationärt eller inte. Förslagsvis, ett mycket svagt avtagande ACF skulle 
## indikera på icke-stationäritet.


par(mfrow=c(2,1))
acf(HOXSWE, main = "ACF Boprisindex | Ursprunglig serie") 
pacf(HOXSWE, main = "PACF Boprisindex | Ursprunglig serie ")
## ACF uppvisar successivt avtagande s.k. spikar
## vid laggarna. Genom visuell inspektion kan vi konstatera
## att serien inte är stationär. 


## Simuleringen av autokorrelations(ACF)-samt partiella autokorrelationsplottar (PACF)
## dock primärt ACF indikerar att ursprungsserien "HOXSWE" är icke-stationär,
## som därtill karaktäriseras av autokorrelation. På grund av detta
## behöver serien genomgå en transformation. Då vi arbetar med ARIMA modeller
## i detta projekt, ska vi göra serien differensstationär genom att
## ta första differensen. Värt att notera, då vi har månadsdata (S=12),
## kan det finnas indikation på att en säsongs ARIMA (SARIMA) kan bli 
## aktuellt att modellera d.v.s kan kräva en säsongsdifferensning i detta fall
## utöver första differensen. Vi ska undersöka det senare.


## Genomför första differensen av urprungsserien för att uppnå stationäritet.
## (p,d,q) <- icke-säsongskompontenter.

par(mfrow=c(2,1))
d.HOXSWE <- diff(HOXSWE)
plot(d.HOXSWE, main = "FÖRSTA DIFFERENSEN | Boprisindex HOXSWE", xlab = "Tid",
     col ="black", las=1)
## Efter genomförd differensning, har väntevärdet och variansen av serien blivit 
## oberoende av tiden, dvs mer konstant än tidigare. Då vi arbetar med månadsdata
## går det därtill även att utröna inslag av säsongsvariation med skarpa och
## höga toppar i serien efter genomförda första differensen.

ADF.test(d.HOXSWE)
## Genom att genomföra Augmented Dickey Fuller testet återigen i syfte
## att få det bekräftat ifall vi fortfarande har en enhetsrot eller ej.
## Enligt utskriftet från ADF-testet kan vi utröna att teststatistikan
## för tau3 i absolutbelopp på 2.04 fortfarande understiger valfri signifikansnivå. 
## Vi kan konstatera att vi fortfarande har en enhetsrot. 
## Vi skall nu härnäst genomföra en säsongsdifferensning.


## Till följd av det vi kom fram till efter första differensen ska vi nu
## genomföra en säsongsdifferensning där S = 12 på grund av månadsdata.
## (P,D,Q) <- säsongskomponenter 

diff12.HOXSWE <- diff(d.HOXSWE,12)
plot(diff12.HOXSWE, main = " Säsongs Differensning | Boprisindex HOXSWE", xlab = "Tid",
     col ="red", las=1)
## Efter den genomförda säsongsdifferensning kan vi observera att
## de tidigare nämnda höga toppar har eliminerats mer med undantag
## för skarpa rusningen mot början av 2010. Säsongsdifferensning
## appliceras som en skillnad mellan ett värde och ett värde med lag som 
## är en multipel av "S" i detta fall S=12 därav X_t-12.


## Vi fortsätter vidare med att testa ifall det föreligger en enhetsrot
## på den säsongsdifferensade serien denna gång.
ADF.test(diff12.HOXSWE)
## Teststatistikan tau3 på 3.95 i absolutbelopp överstiger 3.43
## i absolutbelopp på 5% nivån. Vi kan förkasta nollhypotesen
## avseende enhetsrot. Serien är stationär samt har inte längre
## någon enhetsrot.

## Exporterar tablell från ADF-testet på säsongsdifferensade serien
## till latex.
exportADF.diff12.HOXSWE <- ADF.test(diff12.HOXSWE)
exportADF.diff12.HOXSWE 
stargazer(exportADF.diff12.HOXSWE, type ="latex")

####################################################################
## 
##      KPSS-testet (Kwiatkowski-Phillips-Schmidt-Shin)
##
####################################################################

## Nollhypotesen vid ett KPSS-test testar också för stationäritet.

## H0: Stationär
## HA: Icke-stationär 

## Teststatistikan är ett mått på variansen i testekvationen och om
## den är s.k. 'liten', är detta då en indikation på att vår serie
## är stationär. För att serien ska klassificeras som stationär krävs
## det att nollhypotesen i det genomförda KPSS-testet inte förkastas.

ur.kpss(HOXSWE, type="tau")@teststat
## Teststatiska: 0.3179527
ur.kpss(HOXSWE, type ="tau")@cval
## Får fram de kritiska värdena på 10%, 5%, 2.5% och 1% nivån.
## Vi observerar att den erhållna teststatistikan övertiger 
## samtliga kritiska värdena på valfri signifikans nivå.
## Med det sagt, ursprungsdatan är ej stationärt


## Vi undersökar vidare nu ifall den säsongsdifferensade
## serien har blivit stationärt.
ur.kpss(diff12.HOXSWE, type ="tau")@teststat
## Teststatistika: 0.04388372
ur.kpss(diff12.HOXSWE, type = "tau")@cval
## Teststatistikan på 0.04388372 understiger samtliga kritiska värden
## på valfri nivå. 0.04388372 < 0.146 på 5% nivån. Vi kan därav
## konstatera att vi ej kan förksata nollhypotesen avseende stationäritet.
## Den säsongsdifferensade serien är nu stationär på 5% nivån samt
## även valfri vald signifikans nivå.


## På grund genomförandet av första differensen som standardmetodik
## i syfte att uppnå stationäritet och även tillämpning av 
## säsongsdifferensning till följd av misstanke om behovet av säsongskomponenter,
## har vi till följd av första- och säsongsdifferensning erhållit följande
## integrationsordningar: "d = 1" & "D = 1" i modellen ARIMA(p,d,q)(P,D,Q)[s]
## där s=12, som skall estimeras.


#############################################################################
##
## SAMKÖRNING AV FIGURER. 
##
#############################################################################


## Samkör både den första differensade och urprungsdatan (Icke-stationära) ihop.

par(fin=c(100,100))
layout(matrix(c(1,1,2,3), 2, 2, byrow =TRUE))
plot.ts(d.HOXSWE, main = "FÖRSTA DIFFERENSEN OCH URSPRUNGSDATA | HOXSWE", axes = FALSE, col = "blue",
        ylab = "", type = "l", lwd = 1)
axis(side=4)
mtext(side=4, "BOPRISINDEX", at=c(155), line = 1)
axis(side=1)
par(new=TRUE)
plot.ts(HOXSWE, axes = FALSE, col = "red", ylab = "%", lty = 1)
axis(side = 2)
box()
plot.ts(d.HOXSWE, main = "FÖRSTA DIFFERENSEN (%)", xlab="Tid", col = "red")
plot.ts(HOXSWE, main = "BOPRISINDEX HOXSWE", xlab="Tid", col = "blue")


## Samkör både den säsongsdifferensade och urprungsdatan (Icke-stationära) ihop.

par(fin=c(100,100))
layout(matrix(c(1,1,2,3), 2, 2, byrow =TRUE))
plot.ts(diff12.HOXSWE, main = "SÄSONGSDIFFERENSAD OCH URSPRUNGSDATA | HOXSWE", axes = FALSE, col = "blue",
        ylab = "", type = "l", lwd = 1)
axis(side=4)
mtext(side=4, "BOPRISINDEX", at=c(155), line = 1)
axis(side=1)
par(new=TRUE)
plot.ts(HOXSWE, axes = FALSE, col = "red", ylab = "%", lty = 1)
axis(side = 2)
box()
plot.ts(diff12.HOXSWE, main = "SÄSONGSDIFFERENSAD (%)", xlab="Tid",ylab = "%", col = "red")
plot.ts(HOXSWE, main = "BOPRISINDEX HOXSWE", xlab="Tid", ylab = "%", col = "blue")



## Vi ska nu studera ACF och PACF i syfte mekaniskt
## estimera modellordningar för ARIMA-modellerna.
## Till följd av säsongsdifferensningen som vi genomförde
## kommer vi behöv att estimera en SARIMA(p,d,q)(P,D,Q)[12] där S=12
## p.g.a månadsdata.
## Där ordningen för de autoregressiva AR(p) och movinge average 
## komponenterna MA(q) kommer att estimeras genom visuell inspektion
## av ACF och PACF på första differensade correlogram och 
## på så vis estimera (p,d,q) medan samma metodik kommer att
## tillämpas för att estimera SAR(p) och SMA(q) fast
## på den säsongsdifferensade serien.

par(mfrow=c(2,1))
acf(d.HOXSWE, main ="ACF Första Differens")
pacf(d.HOXSWE, main = "PACF Första Differens") 
## Går att observera en signifikant s.k. spik vid första laggen i ACF
## vilket kan indikera en MA(1) samtidigt som man observerar
## spikar utanför konfidensbandet. Vad gäller PACF, betraktar
## man den mest signikanta och därtill intressanta spik vid
## första laggen, kanske en AR(1) d.v.s SARIMA(1,1,1)(P,D,Q)[12].
## Men vi fortsätter att undersöka.


## En tumregel är enligt (Enders,2015, p.216) att sätta laggen till 3 år
## d.v.s för månadsdata till 36. 
par(mfrow=c(2,1))
acf(d.HOXSWE, lag.max= 36, main ="ACF Första Differnens | lagmax 36")
pacf(d.HOXSWE, lag.max= 36, main ="PACF Första Differens | lagmax 36")
## Vi utökar antalet laggar och genom grafisk inspektion kan vi
## observera 2 signifikanta s.k. spikar vid runt lag 1(M12) och 2(M24) I PACF 
## vilket kan indikera en AR(2) medan ACF uppvisar ett successivt dämpande
## mönster. ACF indikerar på behovet av tillämpning av säsongskomponenter.


## Vi utökar antal laggar till något extremt som 50 för fortsatt inspektion.
acf(d.HOXSWE, lag.max = 50, main ="ACF Första Differens | lagmax 50")
pacf(d.HOXSWE, lag.max = 50, main ="PACF Första Differens| lagmax 50")
## ACF fortsätter att bekräfta ett successivt avtagande mönster vid varje enskild lag.
## PACF bekräftar  återigen på en AR(2)


## Ännu en gång utökar vi antal laggar till något extremt som 100 för inspektion.
acf(d.HOXSWE, lag.max = 100, main ="ACF Första Differens | lagmax 100")
pacf(d.HOXSWE, lag.max = 100, main ="PACF Första Differens| lagmax 100")

## Utan funktionen "lag.max", ser det ut som en ARIMA(1,1,1)
## integrerad av ordning 1 där d=1 med varsin autoregressiv komponent
## samt moving average komponent.
## När vi väl nyttjar funktionen "lag.max" kan vi observera
## av avtagande mönster i ACF vid varje enskild lagg vilket
## indikerar på säsongsvariation samtidigt PACFB
## uppvisar två signifikanta spikar vid första och andra laggen.
## Till följd av ACF successiva avtagande karaktär, kan det 
## uppfattas som ARIMA(2,1,0) men vi fortsätter med undersökningen.


## Men, på grund av indikation av säsongsvariation, ska vi nu undersöka samt
## mekaniskt estimera säsongskomponenterna (P,D,Q) ur SARIMA (p,d,q)(P,D,Q)[12].
## Vid estimation av säsongskomponenterna ska vi studera mönster vid de enskilda
## laggarna. För vårt månadsdata, ska vi studera vid lag 1 (M12), 2 (M24) 
## och 3 (M36) 

par(mfrow=c(2,1))
acf(diff12.HOXSWE, main ="ACF Säsongs Differensad")
pacf(diff12.HOXSWE, main = "PACF Säsongs Differensad") 
## ACF uppvisar ett kort segment med avtagande spikar
## som sedan faller in i bandet, dock uppvisas en signfikant
## spik borta vid lag 1 (M12).


## En tumregel är enligt (Enders,2015, p.216) att sätta laggen till 3 år
## d.v.s för månadsdata till 36. 
par(mfrow=c(2,1))
acf(diff12.HOXSWE, lag.max = 36, main ="ACF Säsongs Differensad | lagmax 36")
pacf(diff12.HOXSWE, lag.max = 36, main ="PACF Säsongs Differensad | lagmax 36")



## Utökar antal laggar till något extremt som 50 för inspektion.
par(mfrow=c(2,1))
acf(diff12.HOXSWE, lag.max = 50, main ="ACF Säsongs Differensad | lagmax 50")
pacf(diff12.HOXSWE, lag.max = 50, main ="PACF Säsongs Differensad | lagmax 50")



## Utökar antal laggar till något extremt som 100 för inspektion och tydlighetens
## skull.
par(mfrow=c(2,1))
acf(diff12.HOXSWE, lag.max = 100, main ="ACF Säsongs Differensad | lagmax 100")
pacf(diff12.HOXSWE, lag.max = 100, main ="PACF Säsongs Differensad | lagmax 100")


## Genom visuell inspektion kan vi observera ett 
## ACF uppvisar ett avtagande mönster under ett kort segment, och faller relativt 
## hastigt under konfidensbandet. 
## Autokorrelationsplotten uppvisar en signifikant spik vid lag 1 (M12).
## vilket kan indikera på en SMA(1) där Q=1 d.v.s
## en säsongs moving average komponent medan PACF på den säsongsdifferensade
## serien uppvisar ett något mer tvetydigt mönster.
## PACF uppvisar s.k. spikar vid lag 1(M12) och lag 3(M36) och likaså
## uppviser PACF intressant mönster vid lag 2 (M12).
## Man skulle kunna utröna ett avtagande mönster från lag 1 till lag 3.


## EXPORTERA till figur.
par(mfrow=c(2,2))
acf(d.HOXSWE, lag.max = 100, main ="ACF Första Differens")
pacf(d.HOXSWE, lag.max = 100, main ="PACF Första Differens")
acf(diff12.HOXSWE, lag.max = 100, main ="ACF Säsongsdifferensad")
pacf(diff12.HOXSWE, lag.max = 100, main ="PACF Säsongsdifferensad")


## Avslutningsvis tillämpar vi AUTO.ARIMA() funktionen som i försig
## väljer den ultimata modellen utifrån informationskriterier.

auto.arima(HOXSWE) # # ARIMA(1,1,1)(0,1,1)[12] # d=1 och D=1 på Ursprungsdata

auto.arima(HOXSWE, trace=TRUE) 
## Vill undersöka grundad på vilka modeller som ARIMA(1,1,1)(0,1,1)[12] valts ut.  


#############################################################################
##
##  2) ESTIMERING (AV OLIKA SARIMA MODELLORDNINGAR)
##
#############################################################################



## ARIMA(2,1,0)(0,1,1)[12]
HOXSWE.arima210011 <- arima(HOXSWE, order = c(2,1,0), seasonal = c(0,1,1))
coeftest(HOXSWE.arima210011) 
## Icke-signfikant AR(2) komponent, är på AR(1) och SMA(1)
## höggradigt signifikanta.


## ARIMA(1,1,0)(0,1,1)
HOXSWE.arima110011 <- arima(HOXSWE, order = c(1,1,0), seasonal = c(0,1,1))
coeftest(HOXSWE.arima110011)  
## Höggradirgt sigfnikanta parameterestimat.
## Den modellordning som ordningen ändras till för 2018:M01
## 2017:M12 samt 2017:M11
## Modellen som förelår modellordningar enligt auto.arima beroende på tidsperiod)


## ARIMA(1,1,1)(0,1,1)
HOXSWE.arima111011 <- arima(HOXSWE, order = c(1,1,1), seasonal = c(0,1,1))
coeftest(HOXSWE.arima111011) 
## Samtliga parameterestimaten är höggradigt signfifikanta.


###########################################################################
##
##  3) Diagnostik
##
###########################################################################


## Jarque-Bera Test (Test av normalitet)
## Detta test tillämpas för att undersöka ifall residualerna
## är normalfördelade.

## H0: Residualerna är normalfördelade
## HA: Residualerna är EJ normalfördelade

jarque.bera.test(residuals(HOXSWE.arima210011)) # P-värde: 0.01718
jarque.bera.test(residuals(HOXSWE.arima110011)) # P-värde: 0.0328
jarque.bera.test(residuals(HOXSWE.arima111011)) # P-värde: 0.02322

## Ingen av de mekaniskt estimerade modellerna uppfyllde antagandet om 
## normalitet bland residualerna.
## Dock så fortsätter vi med diagnostiken trots medvetenheten om det.

## LJUNG-BOX TEST ( Test av autokorrelation i residualer)

## Genom att utföra Ljung-Box Testet för autokorrelation,
## kan vi undersöka ifall det finns kvarvarande autokorrelation
## i residualerna.

## H0: Finns EJ kvarvarande autokorrelation
## HA: Finns kvarvarande autokorrelation.

Box.test(resid(HOXSWE.arima210011), lag=4, type="L", fitdf = 3)  # P-värde = 0.1288
Box.test(resid(HOXSWE.arima111011), lag=4, type="L", fitdf = 3)  # P-värde = 0.2028
Box.test(resid(HOXSWE.arima110011), lag=4, type="L", fitdf = 3)  # P-värde = 0.02802

## Enbart modellen SARIMA(1,1,0)(0,1,1)[12] uppvisar att den har kvarvarande
## autokorrelation då P-värdet på 0.02802 understiger 5% nivån vilket möjliggör
## att nollhypotesen avseende ingen kvarvarande autokorrelation
## i residualerna förkastas.Genom funktionen "coeftest" har vi satt "fitdf" 
## till antal koefficienter justerat för interceptet d.v.s utan interceptet och
## interceptet försvinner till följd av differensning. Så "fitdf" är satt lika med
## antal parametrar.

## Vi väljer att arbeta vidare med den mekaniskt estimerade modellen
## och den modell som auto.arima väljer.

## ARIMA(2,1,0)(0,1,1)[12]
par(mfrow=c(3,2))
plot(residuals(HOXSWE.arima210011), main ="Residualsplot | ARIMA(2,1,0)(0,1,1)[12]")
hist(residuals(HOXSWE.arima210011), main ="ARIMA(2,1,0)(0,1,1)[12]", col = "gray", xlab = "Residual")
acf(residuals(HOXSWE.arima210011), main = "ACF residuals | ARIMA(2,1,0)(0,1,1)[12]")
pacf(residuals(HOXSWE.arima210011), main = "PACF residuals | ARIMA(2,1,0)(0,1,1)[12]")
qqnorm(residuals(HOXSWE.arima210011), main = "Normal Q-Q Plot | ARIMA(2,1,0)(0,1,1)[12]")
abline(a = 0, b=1, h=0, v=0, col=c("gray", "red", "blue"),
       lty = c(1,2,2))


## ARIMA(1,1,1)(0,1,1)[12]
par(mfrow=c(3,2))
plot(residuals(HOXSWE.arima111011), main ="Residualsplot | SARIMA(1,1,1)(0,1,1)[12]")
hist(residuals(HOXSWE.arima111011), main= "Histogram of residuals | SARIMA(1,1,1)(0,1,1)[12]",col = "gray", xlab = "Residual")
acf(residuals(HOXSWE.arima111011), main ="ACF residals | SARIMA(1,1,1)(0,1,1)[12]")
pacf(residuals(HOXSWE.arima111011), main ="PACF residals | SARIMA(1,1,1)(0,1,1)[12]")
qqnorm(residuals(HOXSWE.arima111011), main = "Normal Q-Q Plot | SARIMA(1,1,1)(0,1,1)[12]")
abline(a = 0, b=1, h=0, v=0, col=c("gray", "red", "blue"),
       lty = c(1,2,2))


## Genom grafisk inspektion av QQ ploten, kan vi se att det föreligger
## några extremvärden, värden som ligger betydligt långt ifrån den
## den räta linjän som går ingenom origo. Då vi även tidigare
## förkastade nollhypotesen avseende normalitet vid JB-testet, för samtliga estimerade
## ARIMA-modeller, finns det skäl till att undersöka ifall vi har möjligheten
## att eliminera problemen med residualernas egenskaper. Vi väljer
## att göra det för de två modeller vi har slussat fram grundad på utfallen under
## olika diagnostiktesten.


## ARIMA(2,1,0)(0,1,1)[12] : p-value = 0.1541
RES.210011 <- residuals(HOXSWE.arima210011)
MAX.210011 <- RES.210011 == max(RES.210011)
MIN.210011 <- RES.210011 == min(RES.210011)
MAXMIN.210011 <- cbind(MAX.210011, MIN.210011)
HOXSWE.arima210011.X <- arima(HOXSWE, order = c(2,1,0), seasonal = c(0,1,1), xreg= MAXMIN.210011) 

## Genomför JB-testet återigen.
jarque.bera.test(residuals(HOXSWE.arima210011.X))
## NORMALFÖRDELAT då det observerade p-värdet på 0.1541 överstiger
## 5% nivån. 


## ARIMA(1,1,1)(0,1,1)[12] : p-value = 0.1161
RES.111011 <- residuals(HOXSWE.arima111011)
MAX.111011 <- RES.111011 == max(RES.111011)
MIN.111011 <- RES.111011 == min(RES.111011)
MAXMIN.111011 <- cbind(MAX.111011, MIN.111011)
HOXSWE.arima111011.X <- arima(HOXSWE, order = c(1,1,1), seasonal = c(0,1,1), xreg= MAXMIN.111011) 

## Genomför JB-testet återigen.
jarque.bera.test(residuals(HOXSWE.arima111011.X))
# NORMALFÖRDELAT också.

## Vi kan ej förkasta nollhypotesen avseende att residualerna är normala
## i båda fallen.


#############################################################################
##
## MODELL UTVÄRDERING MED INFORMATIONSKRITERIER
##
#############################################################################


## AIC (Akaike Information Criterion)

SARIMAaic = auto.arima(HOXSWE, trace=TRUE, ic="aic")
# Best model: SARIMA(1,1,1)(0,1,1)[12]  : 460.4947

## BIC:

SARIMAbic = auto.arima(HOXSWE, trace=TRUE, ic="bic")
## Best model: ARIMA(1,1,0)(0,1,1)[12] : 470.2072

## Genom BIC som IC, inser vi däremot att vår mekaniskt estimerade modell
## SARIMA(2,1,0)(0,1,1)[12] kommer på 3:e plats, dock så tar
## SARIMA(1,1,1)(0,1,1)[12] 2:a plats i detta fall istället för 1:a.
## Vi ser även att modellen ARIMA(1,1,0)(0,1,1)[12] är det ultimata enligt
## BIC. Värt att ha i beaktelse, är det facto att då bestraffningen är större,
## så föredras modell med färre parametrar.


## Med hjälp av den rangordning av modeller som informationskriterier genomför
## i syfte att utvärdera de estimerade modeller, går det att observera att
## det föreligger skillnader mellan ARIMA(1,1,1)(0,1,1)[12]
## och ARIMA(2,1,0)(0,1,1)[12] i termer av den rangordning man får.


## ARIMA(1,1,1)(0,1,1)[12]
## ARIMA(2,1,0)(0,1,1)[12]


########################################################################
##
## REKURSIVA PROGNOSER
##
########################################################################

## Vi har använt ett dataset fram till 2016:M01
## trots att hela serien HOX sträcker sig fram till 2018:M01.
## Det vi ska göra är att successivt använda denna hittils utnyttjade informationen.
## Vi ska ponera att vi ständigt får ny information.
## Detta innebär att för varje prognos vi gör fram till och med 2018:M01
## har vi ett utfall att benchmarka med fr.o.m år 2016. Detta arbetssätt
## kallas för pseudo-out-sample-forecasts på engelska, där vi ponerar att
## vi inte hade det faktiska utfallsdatan för de kommande två åren
## efter 2016.

## Sedan tidigare vet vi om att auto.arima(HOXSWE) identifierade
## en SARIMA (1,1,1)(0,1,1)[12] dessutom utifrån informationskriterier.


###################################################################
# REKURSIV PROGNOS - SARIMA (1,1,1)(0,1,1)[12]
#
# 
# auto.arima() -> SARIMA (1,1,1)(0,1,1)[12]
###################################################################


## Vi skapar en s.k. loop som ska i vårt fall innehålla våra prognoser.
## Vi uppmanar loopen att börja skapa prognoser från 2016:M02
## kommande 24 månader.

Prognos.HOXSWE.arima111011  <- ts(matrix(NA, nrow = 24, ncol = 1),
                                  start = c(2016, 2), frequency = 12)
for (i in 1:24){
  EndDate <- 2016.0083 + (i - 1) / 12
  Data <- window(HOXSWE_hela, end = EndDate)
  Result <- arima(Data, order = c(1, 1, 1),
                  seasonal= list(order = c(0, 1, 1)))
  Prognos.HOXSWE.arima111011[i] <- forecast(Result, h = 1)$mean
}
## 1/12 = 0.083 (Varje månad är 0.083-del av ett år)
## 2016.0083 = Januari 2016 dvs(2016:M01). Det var fram till
## 2016:M01 som vi estimerade ARIMA-modellerna och från och med 
## kommer vi låtsas som om vi ej har faktiska utfallsdatan.


Prognos.HOXSWE.arima111011 
## Får fram pseudo-out-sample-forecasts för de kommande 24 månader
## med start ifrån 2016:M02 fram till 2018:M01.


## Studera de kommande prognoserna med utfallsdata
Utfall.från201602<- window(HOXSWE_hela, start = c(2016,02))
window(cbind(Utfall.från201602, Prognos.HOXSWE.arima111011), end = c(2018,02))
ts.plot(Utfall.från201602, main ="Utfallsserie från 2016:M02")
ts.plot(Prognos.HOXSWE.arima111011, main="Psuedo-out-sample-forecasts från 2016:M02")


## Samkör Utfall från 2016 och Prognostiserade SARIMA(1,1,1)(0,1,1)[12]
ts.plot(Utfall.från201602, Prognos.HOXSWE.arima111011, 
        main ="Utfallsdatan & Pseudo SARIMA(1,1,1)(0,1,1)[12]",
        gpars=list(xlab = "Time",col =c("black","red")))
legend("topleft", legend = c("Utfall från 2016", "Pseudo SARIMA111011"), col = c("black", "red"), lty = 1, cex=0.6)


## Samkör hela HOX serien från 2005:M01 till 2018:M01 samt prognostiserade SARIMA(1,1,1)(0,1,1)[12]

## Utfall.hela = hela serien från 2005 till 2018 i fallet för SARIMA(1,1,0)(0,1,1)[12]
Utfall.hela05till18 <- window(HOXSWE_hela, start = c(2005,1))
plot(Utfall.hela05till18, main ="HOX")
window(cbind(Utfall.hela05till18, Prognos.HOXSWE.arima111011))
ts.plot(Utfall.hela05till18 , Prognos.HOXSWE.arima111011, 
        main ="HOX Utfallsdata & Pseudo SARIMA(1,1,1)(0,1,1)[12]", 
        gpars=list(xlab = "Time",col =c("black","red")))

legend("topleft", legend = c("HOXSWE", "Pseudo SARIMA111011"), 
       col = c("black", "red"), lwd =1, lty = 1,  cex=0.6)



#######################################################################
# REKURSIV PROGNOS - SARIMA (2,1,0)(0,1,1)[12]
# 
# Den mekaniskt estimerade modellordningen -> SARIMA (2,1,0)(0,1,1)[12]
#######################################################################

## Vi skapar återigen en loop här igen som ska i vårt fall innehålla 
## våra prognoser.
## Vi uppmanar loopen att börja skapa prognoser från 2016:M02 som ska köra
## kommande 24 månader.

Prognos.HOXSWE.arima210011  <- ts(matrix(NA, nrow = 24, ncol = 1),
                                  start = c(2016, 2), frequency = 12)

for (i in 1:24){
  EndDate <- 2016.0083 + (i - 1) / 12
  Data <- window(HOXSWE_hela, end = EndDate)
  Result <- arima(Data, order = c(2, 1, 0),
                  seasonal= list(order = c(0, 1, 1)))
  Prognos.HOXSWE.arima210011[i] <- forecast(Result, h = 1)$mean
}
## 1/12 = 0.083 (Varje månad är 0.083-del av ett år)
## 2016.0083 = Januari 2016 dvs(2016:M01). Det var fram till
## 2016:M01 som vi estimerade ARIMA-modellerna och från och med nu 
## kommer vi låtsas som om vi ej har faktiska utfallsdatan.


Prognos.HOXSWE.arima210011
## Får fram pseudo-out-sample-forecasts för de kommande 24 månader
## med start ifrån 2016:M02 fram till 2018:M01.


Utfall <- window(HOXSWE_hela, start = c(2016,2))
window(cbind(Utfall.från201602, Prognos.HOXSWE.arima210011), end = c(2018,2))


## Samkör Utfall och Prognostiserade SARIMA(2,1,0)(0,1,1)[12]

ts.plot(Utfall.från201602, Prognos.HOXSWE.arima210011,
        main ="Utfall & Prognostiserade SARIMA(2,1,0)(0,1,1)[12]", 
        gpars=list(xlab = "Time",col =c("black","purple")))
legend("topleft", legend = c("Utfall från 2016", "Pseudo SARIMA210011"), col = c("black", "purple"), lty = 1, cex=0.60)


#######################################################################
## 
## Samkör HELA serien från 2005:M01 till 2018:M01 samt prognostiserade 
## SARIMA(2,1,0)(0,1,1)[12]
##
#######################################################################


plot(Utfall.hela05till18, main="HOX",ylab="%")
window(cbind(Utfall.hela05till18, Prognos.HOXSWE.arima210011), end = c(2018,02))
ts.plot(Utfall.hela05till18, Prognos.HOXSWE.arima210011, 
        main ="Hela Utfall & Prognostiserade SARIMA(2,1,0)(0,1,1)[12]",
        gpars=list(xlab = "Time", col =c("black","purple")))
legend("topleft", legend = c("HOXSWE", "Pseudo SARIMA210011"), 
       col = c("black", "purple"), lty = 1, cex=0.6)


#######################################################################
## 
##  Samkör utfallsdata, SARIMA(1,1,1)(0,1,1) och SARIMA(2,1,0)(0,1,1)
##
#######################################################################

ts.plot(Utfall.från201602, Prognos.HOXSWE.arima210011, Prognos.HOXSWE.arima111011,
        main ="Utfallsdata & SARIMA111011 & SARIMA210011", 
        gpars=list(xlab = "Time",col =c("black","purple", "red")))
legend("topleft", legend = c("Utfallsdata från 2016", "Pseudo SARIMA210011",
                             "Pseudo SARIMA111011"), col = c("black", "purple" ,"red"), lty = 1, cex=0.55)

## Båda av våra estimerade modeller SARIMA(1,1,1)(0,1,1) och SARIMA(2,1,0)(0,1,1)
## följer ursprungsserien väl. Mycket marginella skillnader
## mellan de genomföra pseudo-out-of-forecasts av modellerna
## SARIMA(1,1,1)(0,1,1) och SARIMA(2,1,0)(0,1,1).



#######################################################################
## Optimalitets Test.
##
## Nu ska vi testa även för vilken modellordning som är optimal 
## d.v.s unbiased!
## Att prognosen är icke-tendentiös är ett viktigt villkor som är
## nödvändigt för en prognos att vara optimal.
#######################################################################

Prognos.arima.1B <- ts(matrix(NA, nrow = 24, ncol = 1),
                       start = c(2016, 2), frequency = 12)


M.Order2 <- ts(matrix(NA, nrow = 24, ncol = 7),
               start = c(2016, 2), frequency = 12,
               names=c("p","q","P","Q","F","d","D"))

for (i in 1:24){
  EndDate <- 2016.0083 + (i - 1) / 12
  Data <- window(HOXSWE_hela, end = EndDate)
  Result <- auto.arima(Data)
  Prognos.arima.1B[i] <- forecast(Result, h = 24)$mean
  M.Order2[i,] <- Result$arma
}

M.Order2
## Visar sig vara samma modellordning förutom sista 3 månader.
## Föreslår SARIMA(1,1,1)(0,1,1)[12] från 2016:M02 till 2017:M10.
## Från 2017:M11 ändras modellordningen till SARIMA(1,1,0)(0,1,1)[12].

## Värt att notera, modellen SARIMA(1,1,0)(0,1,1)[12] estimerades
## tidigare i projektet dock fallerade modellen både på
## normalitet samtidigt som den led av kvarvarande
## autokorrelation i residualerna. Därav valde vi inte att
## arbeta med den modell p.g.a den inte uppfyllde
## antagandena vid genomförandet av standardtesterna.


window(cbind(Utfall.från201602, Prognos.HOXSWE.arima111011 , Prognos.arima.1B),
       start = c(2016, 2), end = c(2018, 1))
## För månaden 2017:M11, fångar den rekursiva upp marginellt bättre.
## För månaden 2017:M12 fångar den rekursiva SARIMA(1,1,1)(0,1,1)
## något bättre upp än SARIMA(1,1,0)(0,1,1).
## Däremot för månaden 2018:M01 framstås SARIMA(1,1,0)(0,1,1)
## som att det fångar upp bättre då punktprognosen på 225.4712
## är närmre utfallsdatan på 228.71 än den rekursiva SARIMA(1,1,1)(0,1,1)
## vars punktprognos är 224.7824



#################################################################
# RULLANDE PROGNOSER
#
#################################################################

## Vi prövar även med att tillämpa rullande prognoser i syfte
## med att jämföra med rekursiva prognoser. Vid tillämpning av
## rullande prognoser, kommer nya observationer att adderas
## samtidigt som de äldsta observationerna exkluderas. På så vis
## kommer datasetet vara ideligen oförändrad.



Prognos.arima.rullande <- ts(matrix(NA, nrow = 24, ncol = 1),
                             start = c(2016, 2), frequency = 12)

M.Order4 <- ts(matrix(NA, nrow = 24, ncol = 7),
               start = c(2016, 2), frequency = 12,
               names=c("p","q","P","Q","F","d","D"))

for (i in 1:24){
  StartDate <- 2005.0083 + (i-1)/12
  EndDate <- 2016.0083 + (i-1)/12
  Data <- window(HOXSWE_hela, start = StartDate, end = EndDate)
  Result <- auto.arima(Data)
  Prognos.arima.rullande[i] <- forecast(Result, h = 24)$mean
  M.Order4[i,] <- Result$arma
}

M.Order4
window(cbind(Utfall.från201602,Prognos.arima.1B,
             Prognos.arima.rullande, Prognos.HOXSWE.arima111011, Prognos.HOXSWE.arima210011), end = c(2018, 1))

Samtligaprognoser <- window(cbind(Utfall.från201602,Prognos.arima.1B,
                                  Prognos.arima.rullande, Prognos.HOXSWE.arima111011,
                                  Prognos.HOXSWE.arima210011), end = c(2018, 1))


samtligaprognoser <- as.data.frame(window(cbind(Utfall.från201602,Prognos.arima.1B,
                                                Prognos.arima.rullande,
                                                Prognos.HOXSWE.arima111011, Prognos.HOXSWE.arima210011), end = c(2018, 1)))

## Exportera utfallsdata fr. 2016:M02 till 2018:M01,
## rekursiv prognos SARIMA(1,1,1)(0,1,1)[12], Pseudo 1B (modellen som förelår
## modellordningar enligt auto.arima beroende på tidsperiod) samt
## rullande prognos.
stargazer(samtligaprognoser, type="latex", summary=FALSE, font.size = "scriptsize")



#########################################################################
## Samkör serier på:
##
## Utfallsdata från 2016 till 2018,
## utfallsdata från 2005 till 2018, rullande prognos
## Pseudo 1B(rekursiv prognos där modellordningen ändras enligt auto.arima
## samt rekursiva SARIMA(1,1,1)(0,1,1)[12] och SARIMA(2,1,0)(0,1,1)[12]
##########################################################################

window(cbind(Utfall.hela05till18, Prognos.HOXSWE.arima111011, Prognos.arima.1B, 
             Prognos.arima.rullande, Prognos.HOXSWE.arima210011), end = c(2018, 1))

Allihopa <- window(cbind(Utfall.hela05till18, Prognos.HOXSWE.arima111011, 
                         Prognos.arima.1B, Prognos.arima.rullande, Prognos.HOXSWE.arima210011), end = c(2018, 1))

Allihopa
## Hela utfallsserien från 2005-2018 samt de rekursiva
## pseudo-out-sample-forecasts, rullande prognos.



## SAMKÖR ALLA 4 serier ihop med HOX från 2016 till 2018.
ts.plot(Utfall.från201602, Prognos.HOXSWE.arima111011,Prognos.arima.1B, 
        Prognos.arima.rullande, Prognos.HOXSWE.arima210011,  main ="Samkör samtliga serier fr. 2016",
        gpars=list(xlab = "Time", col =c("black","purple", "red", "green", "pink")))

legend("topleft", legend = c("Utfall från 2016", "SARIMA111011", "SARIMA210011", "Pseudo 1B","Rullande Prognos"), 
       col = c("black", "purple", "red", "green", "pink"), lty = 1, cex=0.6)



## SAMKÖR ALLA 4 serier ihop med HOX från 2005 till 2018.
ts.plot(Utfall.hela05till18, Prognos.HOXSWE.arima111011,Prognos.arima.1B, 
        Prognos.arima.rullande,Prognos.HOXSWE.arima210011,  main ="Samkör samtliga serier fr.2005",
        gpars=list(xlab = "Time", col =c("black","purple", "red", "green", "pink")))

legend("topleft", legend = c("Utfall HOX", "SARIMA111011", "SARIMA210011", "Pseudo 1B", " Rullande Prognos"), 
       col = c("black", "purple", "red", "green", "pink"), lty = 1, cex=0.6)



###################################################################
# PROGNOSUTVÄRDERING
#
#
################################################################### 

## Rekursiv -  SARIMA(2,1,0)(0,1,1)[12]
FE.Prognos.HOXSWE.arima111011 <- Utfall.från201602 - Prognos.HOXSWE.arima111011 

## Rekursiv -  SARIMA(2,1,0)(0,1,1)[12]
FE.Prognos.HOXSWE.arima210011 <- Utfall.från201602 - Prognos.HOXSWE.arima210011

## Rekursiv - Pseudo 1B: där modellordningen ändras utifrån auto.arima vid olika tidsperioder
FE.Prognos.arima.1B <- Utfall.från201602 - Prognos.arima.1B

## Rullande
FE.Prognos.arima.rullande <- Utfall.från201602 - Prognos.arima.rullande

###############################################################################

coeftest(dynlm(FE.Prognos.HOXSWE.arima111011 ~1), vcov. = NeweyWest) 
## SARIMA(1,1,1)(0,1,1)[12] med ett P-värde:  0.6624 
## (Prognosen är unbiased på 5% nivån.)


coeftest(dynlm(FE.Prognos.arima.1B  ~1), vcov. = NeweyWest) 
## P-värde: 0.5979 (Prognosen är unbiased)


coeftest(dynlm(FE.Prognos.HOXSWE.arima210011 ~ 1 ), vcov. = NeweyWest) 
## SARIMA(2,1,0)(0,1,1)[12] med ett P-värde: 0.6523
## (Prognosen är unbiased på 5% nivån.)
## trots att Pseudo 1B modell aldrig föreslog SARIMA(2,1,0)(0,1,1)[12] 
## någongång.

coeftest(dynlm(FE.Prognos.arima.rullande  ~ 1 ), vcov. = NeweyWest) 
## P-värde: 0.5046 (Prognosen är unbiased)


## I följande sekvens avseende prognosprecision, skall vi testa
## hypotesen om lika prognosprecision genom Diebold-Mariano-testet.

## H0: förlustdifferensen är noll 
## HA: förlustdifferensen är skild från noll 

prognoser <- c("FE.Prognos.HOXSWE.arima111011", "FE.Prognos.arima.1B",
               "FE.Prognos.arima.rullande", "FE.Prognos.HOXSWE.arima210011")

(testlista <- combn(prognoser, 2))
## Vi får en lista på samtliga prognoser
## samt därtill en list på parvisa testerna vi vill göra.


## Diebold och Mariano

DM.TEST <- function(x){
  require(dynlm)
  require(lmtest)
  require(sandwich)
  fe1 <- get(x[1])
  fe2 <- get(x[2])
  res <- dynlm(I(fe1^2 - fe2^2) ~ 1)
  round(coeftest(res, vcov. = NeweyWest)[4], 3)
}

apply(testlista, 2, DM.TEST)

## Ser att ingen av prognoserna förkastar nollhypotesen
## avseende lika förlustdifferens.
## Då samtliga p-värden överstiger den valda
## signifikansnivå på 5%, kan vi inte förkasta nollhypotesen
## om lika förlustdifferens avsende prognosfelen.


## För undersökning av prognosprecision, kan även standardiserade
## statistika mått såsom RMSE & MAE tillämpas.

accuracy(HOXSWE.arima111011)
## RMSE: 1.4512, MAE: 1.065833

accuracy(HOXSWE.arima210011)
## RMSE: 1.459459, MAE:1.070971

## Det marginella lägre root mean squared error och
## mean absolue error statistikor, ger indikation på
## SARIMA(1,1,1)(0,1,1)[12] är den datagenereringsprocess
## med lägst prognosfel och därmed högsta prognosprecision
## sinsemellan de två modellerna.

#######################################################################
# END
#######################################################################

