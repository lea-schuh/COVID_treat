# COVID_treat
contains code accompanying the paper: Schuh et al., (2024) medRxiv, https://doi.org/10.1101/2024.05.31.24308284  

## Within-host dynamics of antiviral treatment for SARS-CoV-2 infection
Lea Schuh<sup>1,\*</sup>, Peter V. Markov<sup>1,2</sup>, Ioanna Voulgaridi<sup>3</sup>, Zacharoula Bogogiannidou<sup>3</sup>, Varvara A. Mouchtouri<sup>3</sup>, Christos Hadjichristodoulou<sup>3</sup>, Nikolaos I. Stilianakis<sup>1,4,\*</sup> \
<sup>1</sup> Joint Research Centre (JRC), European Commission, Ispra, Italy \
<sup>2</sup> London School of Hygiene & Tropical Medicine, University of London, London, United Kingdom \
<sup>3</sup> Laboratory of Hygiene and Epidemiology, University of Thessaly, Larissa, Greece 
<sup>4</sup> Department of Biometry and Epidemiology, University of Erlangen–Nuremberg, Erlangen, Germany \
<sup>\*</sup> corresponding authors \
\
e-mails: lea.schuh@ec.europa.eu, nikolaos.stilianakis@ec.europa.eu

### Summary
The effectiveness of antiviral treatment with remdesivir against COVID-19 has been investigated in clinical trials suggesting earlier recovery. However, this effect seems to be rather modest. In this study, we assessed the clinical course of SARS-CoV-2 infections in 369 COVID-19 individuals across a spectrum of illness severities, including both untreated individuals and individuals who received antiviral treatment with remdesivir. Moreover, using a process-based mathematical model, we quantified and analyzed the within-host infection dynamics of 69 untreated and 19 antiviral-treated individuals. For untreated individuals, we found that those hospitalized exhibit significantly lower levels of early immune
response and higher cumulative viral loads than those who were not. For treated individuals, we found that those who died were on average hospitalized later after symptom onset than those who survived, underscoring the importance of early medical intervention for severe COVID-19. Our model estimates a rather limited antiviral activity of remdesivir and, consequently, comparable viral load dynamics between individuals responding and not responding to antiviral treatment. Our results provide valuable insights into the clinical course of COVID-19 during antiviral treatment with remdesivir and suggest the need for alternative treatment regimens. \
\
Required software: MATLAB (R2023a) 

### Folder structure
<ul>
  <li>Data</li>
  <li>Figures</li>
  <li>Results</li>
  <li>Scripts</li>
</ul> 

### How to run the code
The important script is main_script.m for running the workflow (parameter estimation and downstream analysis) and for plotting the resulting figures. 

#### Data
Some legal restrictions apply to the clinical and epidemiological data. They are available to interested researchers upon request (contact: xhatzi@med.uth.gr, mouchtourib@uth.gr).

#### Figures and ExtendedDataFigures
All raw and final subfigures and figures in png, pdf, and svg formats. 

#### Results
The individual-specific estimated parameters with and without treatemnt using our model. 

#### Scripts
All functions used by the main_script.m can be found here. More information on the functions can be found in the scripts themselves eg. input, output.

