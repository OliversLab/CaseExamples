
install.packages("ggplot2", "gridExtra")

## Alles markieren und Knopf "Run" drücken 
## ODER
## Jede Zeile einzeln nacheinander ablaufen lassen mittels: STRG+ENTER

library(ggplot2)
library(gridExtra)

## Aus: Population pharmcokinetics and pharmacodynamics of piperacillin/tazobactam 
## in patients with complicated intra-abdominal infection
## Li et al. Journal of Antimicrobial Chemotherapy (2005) 56, 388-395

piperacillin_single_dose <- function(par = c(ETA1=0,           # par ist ein Vektor, enthält
                                             ETA2=0),          # Benutzerdefinierte Abweichung vom typischen Patienten
                                     omg = c(OMEGA1=0.277,     # Die im popPK Modell hinterlegten Standardabweichungen
                                             OMEGA2=0.252),    # der Modellparameter CL und Vd werden hier abgelegt
                                     rem = c(PROP=0.185,       # Das residuale Fehlermodell mit proportionalem
                                             ADD=1.77),        # und additivem Fehler stammt aus dem pop PK Modell, hier nicht benötigt
                                     amt,                      # amt: Dosis in mg
                                     dur,                      # dur: Infusionsdauer (Tinf) in Stunden
                                     time,                     # time:zu simulierende Zeitpunkt(e) oder Zeitspanne in Stunden
                                     CLcr,                     # CLcr: Creatinin-Clearance des Individuums (mL/min)
                                     WT,                       # WT: Körpergewicht des Individuums (kg)
                                     ID,                       # ID: Fortlaufende Nummerierung des Individuums
                                     noIIV = F) {              # noIIV: Flag schaltet inter-individuelle Variabilität (FALSE) an und aus (TRUE)
  
  ## Kovariatenmodell
  ## Typische Werte für Clearance und Verteilungsvolumen
  ## Werden anhand der popPK Modellparameter aus Kovariaten berechnet
  popCL <- 5.05 + 9.60 * CLcr/89
  popVd <- 22.3 * WT/81.8
  
  ## Statistisches Modell
  ## Wenn IIV ausgeschaltet ist
  if(noIIV){
    ecl = 0
    evd = 0
  } else { 
    ## sonst werden zufällige Werte für ecl und evd aus einer Normalverteilung mit
    ## Mittelwert = 0 und Standardabweichung = OMEGA1 für ecl und OMEGA2 für evd
    ecl = rnorm(n=1, mean=0, sd=(omg[["OMEGA1"]]))
    evd = rnorm(n=1, mean=0, sd=(omg[["OMEGA2"]]))
  }
  
  ## Da aus dem popPK Modell hervorgeht, dass die Modellparameter log-normal verteilt sind
  ## wird die IIV als exponent der e-Funktion angegeben
  ## der Exponent stammt nun entweder aus der Benutzereingabe (par)
  ## oder aus der Normalverteilung oben (ecl und evd)
  indCL <- popCL * exp(ecl+par[["ETA1"]])
  indVd <- popVd * exp(evd+par[["ETA2"]])
  
  
  ## Strukturmodell für ein 1-Kompartiment Modell mit Infusion (0. Ordnung)
  ## und Elimination 1. Ordnung
  ke = indCL/indVd
  
  conc <-  ifelse(time <= dur,    
                  amt/(dur*ke*indVd)*(1-exp(-ke*time)),
                  amt/(dur*ke*indVd)*(1-exp(-ke*dur))*exp(-ke*(time-dur))
  )
  
  ## Alle Ergebnisse und Daten in einer "Tabelle" zusammenfassen...
  results <- data.frame(TIME=time,
                        AMT=amt,
                        DUR=dur,
                        CLCR=CLcr,
                        WT=WT,
                        IPRED=conc,
                        ID=ID,
                        indCL=indCL,
                        indVd=indVd,
                        popCL=popCL,
                        popVd=popVd)
  
  ##... und als Gesamtergebnis der Funktion zurückgeben
  return(results)  
}

## Zeit unseren Patienten zu simulieren

TIME <- seq(0,8,by=0.01)   ## Zeitpunkte zu simulieren
AMT <- 4000                ## Applizierte Dosis in mg
DUR <- 1                   ## Infusionsdauer (Tinf) in Stunden
CLcr <- 131                ## Kreatininclearance in mL/h
WT <- 75                   ## Patientengewicht in kg
MIC <- 2                   ## minimale Hemmkonzentration in mg/L
PK_TARGET <- 4 * MIC       ## PK/PD-Ziel 4x MHK

## Simulation eines typischen Patienten mit den Kovariaten unseres Individuums
## interindividuelle Variabilität wird hier ausgeschaltet, da wir
## nur am wahrscheinlichsten Verlauf der Kurve interessiert sind
## => Zeile "noIIV = T" 

pop_prediction <- piperacillin_single_dose(amt=AMT, 
                                           dur=DUR,  
                                           time=TIME, 
                                           CLcr = CLcr, 
                                           WT=WT, ID=0,
                                           noIIV = T)

## Ergebnisse werden in einer Abbildung dargestellt

pl_1 <- ggplot(data = pop_prediction) +                        ## Definiere Daten für die Abbildung
  geom_line(mapping = aes(x=TIME, y=IPRED)) +          ## Art der Abbildung und was soll auf welche Achse
  theme_bw() +                                         ## Äußerlichkeiten der Abbildung
  geom_hline(mapping = aes(yintercept=PK_TARGET), linetype=2) +            # Eintragen des PK/PD-Zieles
  annotate(geom = "text", x=7, y=PK_TARGET*2, label="4-fache MHK", size=6 ) +   # und beschriften
  xlab("Zeit seit Dosis [h]") +                        ## Beschriftung von x- ...
  ylab("Konzentration [mg/L]") +                         ## und y-Achse
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), plot.title = element_text(size=18), plot.subtitle = element_text(size=16)) + 
  annotate(geom= "text", x=4, y= 170, label= paste("typische Piperacillin-Clearance:", round(pop_prediction[1,]$indCL,2), "L/h"), size=6) +
  annotate(geom= "text", x=4, y= 150, label= paste("typisches Verteilungsvolumen:", round(pop_prediction[1,]$indVd,2), "L"), size=6)


## Abbildung anzeigen

pl_1 

## MC Simulation
## Durch mehrfache Simulieren des gleichen Patienten
## Werden innerhalb der Funktion jedesmal andere Werte für ecl und evd 
## aus der Normalverteilung "gezogen" -> rnorm(...) zieht zufällige Stichproben aus der Normalverteilung
## Dies kann einige Minuten dauern

mc_data <- NULL ## Für die R-Profis: Das ist nicht besonders effizient, 
                ## der Speicher sollte vor der Schleife zugewiesen werden

for(i in 1:1000){
  tmp_data <- piperacillin_single_dose(amt=AMT, dur=DUR, time=TIME, CLcr = CLcr, WT=WT, ID=i)
  mc_data <- rbind(mc_data, tmp_data)
}


## Modellparameter (CL und Vd) der 1000 virtuellen Patienten aus den Daten herausziehen

mc_pop_cl <- (mc_data[mc_data$TIME==0,]$indCL)
mc_pop_vd <- (mc_data[mc_data$TIME==0,]$indVd)


## Histogramme vorbereiten

pl_dens_cl <- ggplot() + geom_histogram(aes(x=mc_pop_cl, y=..density..), fill="white", colour="black", bins=50) + theme_bw() +
  xlab("CL [L/h]") + ylab("Häufigkeit") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), plot.title = element_text(size=18), plot.subtitle = element_text(size=16)) +
  theme(axis.text.y = element_blank())
pl_dens_vd <-ggplot() + geom_histogram(aes(x=mc_pop_vd, y=..density..), fill="white", colour="black", bins=50) + theme_bw() +
  xlab("Vd [L]") + ylab("Häufigkeit") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), plot.title = element_text(size=18), plot.subtitle = element_text(size=16)) +
  theme(axis.text.y = element_blank())

## Aus dem generierten Datensatz von 1000 Individuuen werden nur die 
## Konzentrationswerte herausgezogen

all_IPRED <- NULL

for(i in 1:max(mc_data$ID)){
  all_IPRED <- rbind(all_IPRED,mc_data[mc_data$ID==i,]$IPRED)
}

## Und die Quantile zur Bildung der 95, 90, 85 und 80 % Intervalle
## Werden berechnet

s <- apply(all_IPRED,2,function(x) quantile(x,probs=c(0.025, 0.05, 0.075, 0.10, 0.9, 0.925, 0.95, 0.975, 0.5)))

## Alle Daten in eienr "Tabelle" zusammenfassen
mc_data_2 <- data.frame(TIME=TIME,
                        s1=s[1,],s2=s[8,], # 95% 
                        s3=s[2,],s4=s[7,], # 90%
                        s5=s[3,],s6=s[6,], # 85%
                        s7=s[4,],s8=s[5,], # 80%
                        median=s[9,])      # median => Sollte sich mit der pop_prediction decken

## Abbildung mit Vorhersageintervall erstellen

pl_2 <- ggplot(data=mc_data_2) + geom_line(mapping = aes(x=TIME, y=median)) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s1, ymax=s2), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s3, ymax=s4), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s5, ymax=s6), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s7, ymax=s8), alpha=0.15) + 
  geom_hline(mapping = aes(yintercept=PK_TARGET), linetype=2) +           
  annotate(geom = "text", x=7, y=PK_TARGET*2, label="4-fache MHK", size=6 ) +   
  theme_bw() + xlab("Zeit seit Dosis [h]") + ylab("Konzentration [mg/L]") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), plot.title = element_text(size=18), plot.subtitle = element_text(size=16))

## Abbildung anzeigen
grid.arrange(pl_2, grid.arrange(pl_dens_vd, pl_dens_cl, nrow=2), ncol=2, widths=c(2,1))

## Jetzt kommt TDM hinzu

## Konzentrationsbestimmung nach 2 und 5 Stunden
tdm_data <- data.frame(conc=c(50, 12), # gemessene Konzentration (mg/L)
                       times=c(2, 5))   # zum Zeitpunkt x in Stunden

pl_3 <- ggplot(data=mc_data_2) + geom_line(mapping = aes(x=TIME, y=median)) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s1, ymax=s2), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s3, ymax=s4), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s5, ymax=s6), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s7, ymax=s8), alpha=0.15) + 
  geom_point(data=tdm_data, mapping = aes(x=times, y=conc), size=3, colour="red") +   ## <- Messwerte "Neues Wissen"
  geom_hline(mapping = aes(yintercept=PK_TARGET), linetype=2) +           
  annotate(geom = "text", x=7, y=PK_TARGET*2, label="4-fache MHK", size=6 ) + 
  theme_bw() + xlab("Zeit seit Dosis [h]") + ylab("Konzentration [mg/L]") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), plot.title = element_text(size=18), plot.subtitle = element_text(size=16))

## Abbildung anzeigen
grid.arrange(pl_3, grid.arrange(pl_dens_vd, pl_dens_cl, nrow=2), ncol=2, widths=c(2,1))

## Anpassung nach dem maximum a posteriori Bayesischen Prinzip
## => Ermittlung der Parameter die am "besten" zu den neuen Daten (TDM-Daten)
## Als auch zum popPK Modell passen

## Aufstellen einer Objektiven Funktion, die es zu minimieren gilt

obj_func <- function(par,               # Hier werden verschiedene Werte für ETA1 und ETA2 "ausprobiert"
                     tdm_data,          # TDM Datensatz
                     dosage_data,       # Dosierungsdaten
                     patient_data,      # Patientencharakteristika
                     omg = c(OMEGA1=0.277,     # Die im popPK Modell hinterlegten Standardabweichungen
                             OMEGA2=0.252),    # der Modellparameter CL und Vd werden hier abgelegt
                     rem = c(PROP=0.185,       # Das residuale Fehlermodell mit proportionalem
                             ADD=1.77))        # und additivem Fehler stammt aus dem pop PK Modell, hier nicht benötigt)      
{
  
  ## Extrahieren der TDM-Daten
  measured_conc <- tdm_data$conc
  times <- tdm_data$times
  
  simulated_conc = piperacillin_single_dose(par=par, # Die aktuell probierte Kombination aus ETA1 und ETA2
                                            amt=dosage_data[["AMT"]], 
                                            dur=dosage_data[["DUR"]], 
                                            time=times,                     ## Nur die TDM Zeitpunkte müssen Simuliert werden
                                            CLcr = patient_data[["CLcr"]], 
                                            WT=patient_data[["WT"]], 
                                            ID=0,
                                            noIIV=T)$IPRED
  
  ## Berechnung des Funktionswertes:
  ## Um den Residuellen Fehler gewichtete Abweichungsquadrate für jeden Messwert werden summiert
  ## Dazu wird die Summe der mit der Varianz OMEGA1² und OMEGA2² gewichteten Abweichung der 
  ## Werte in ETA1 und ETA2 (=Abweichung vom typischen Patienten) addiert
  
  sig2 <- simulated_conc^2*rem[["PROP"]]+rem[["ADD"]]
  OFV <- sum( (measured_conc-simulated_conc)^2/(sig2)+log(sig2)) + sum( ((par-0)^2)/(omg^2)) ## <- Hier wird der OFV berechnet
  
  ## Bedeutet also, ein ETA darf umso weiter abweichen, je größer die im Modell hinterlegte Standardabweichung OMEGA ist
  ## => Zu starke Abweichung vom Modell wird dadurch "bestraft", dass ein hoher Funktionswert errechnet wird
  
  return(OFV)
}



## par beinhaltet jetzt die Startwerte für ETAs (0 ist das wahrscheinlichste ETA => modus der Verteilung)
## Die R-Funktion "optim" miniert die in "fn" übergebene objektive Funktion
optim_res <- optim(par=c(ETA1=0, ETA2=0), 
                   fn=obj_func, 
                   tdm_data=tdm_data,
                   dosage_data=c(AMT=AMT, DUR=DUR),
                   patient_data=c(CLcr=CLcr, WT=WT))


## Einsetzen der Ermittelten Werte für ETA1 und ETA2, welche
## die Objektive Funktion minimieren
## Es handelt sich jeweils um den Modus der a posteriori Verteilung von ETA1 und ETA2

ind_prediction <- piperacillin_single_dose(amt=AMT, 
                                           par = optim_res$par, ## <- Ermittelte Optimale Parameter
                                           dur=DUR,  
                                           time=TIME, 
                                           CLcr = CLcr, 
                                           WT=WT, ID=i,
                                           noIIV = T)


## Die Individuelle Vorhersage mit in die Abbildung einbinden
## Die Messwerte liegen NICHT auf der individuellen Vorhersage!
## => Das Modell selbst (a priori Wissen) begrenzt das Ausmaß der Anpassung
## Zu starke Abweichung vom Modell ist unwahrscheinlich

pl_4 <- ggplot(data=mc_data_2) + geom_line(mapping = aes(x=TIME, y=median)) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s1, ymax=s2), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s3, ymax=s4), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s5, ymax=s6), alpha=0.15) +
  geom_ribbon(mapping = aes(x=TIME, ymin = s7, ymax=s8), alpha=0.15) + 
  geom_point(data = tdm_data, mapping = aes(x=times, y=conc), size= 3, colour="red") +
  geom_line(data = ind_prediction, mapping = aes(x=TIME, y=IPRED), colour="red") + ## <- Individuelle Anpassung
  geom_hline(mapping = aes(yintercept=PK_TARGET), linetype=2) +           
  annotate(geom = "text", x=7, y=PK_TARGET*2, label="4-fache MHK", size=6 ) + 
  annotate(geom= "text", x=4, y= 170, label= paste("ind. Piperacillin-Clearance:", round(ind_prediction[1,]$indCL,2), "L/h") , size=6, colour="red") +
  annotate(geom= "text", x=4, y= 150, label= paste("ind. Verteilungsvolumen:", round(ind_prediction[1,]$indVd,2), "L") , size=6, colour="red" ) +
  theme_bw() + xlab("Zeit seit Dosis [h]") + ylab("Konzentration [mg/L]") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), plot.title = element_text(size=18), plot.subtitle = element_text(size=16))

## Vorbereiten der Histogramme
pl_dens_vd_extended <- pl_dens_vd + geom_vline(aes(xintercept=ind_prediction[1,]$indVd), size=3, alpha=0.5, colour="red")
pl_dens_cl_extended <- pl_dens_cl + geom_vline(aes(xintercept=ind_prediction[1,]$indCL), size=3, alpha=0.5, colour="red") 


## Abbildung anzeigen
grid.arrange(pl_4, grid.arrange(pl_dens_vd_extended, pl_dens_cl_extended, nrow=2), ncol=2, widths=c(2,1))


#######################################
## Abbildungen aus dem Artikel erzeugen

## Abbildung A
pl_1 

## Abbildung B
grid.arrange(pl_2, grid.arrange(pl_dens_vd, pl_dens_cl, nrow=2), ncol=2, widths=c(2,1))

## Abbildung C
grid.arrange(pl_3, grid.arrange(pl_dens_vd, pl_dens_cl, nrow=2), ncol=2, widths=c(2,1))

## Abbildung D
grid.arrange(pl_4, grid.arrange(pl_dens_vd_extended, pl_dens_cl_extended, nrow=2), ncol=2, widths=c(2,1))
