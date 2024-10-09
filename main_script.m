%%  main script for manuscript: Within-host dynamics of antiviral treatment for SARS-CoV-2 infection
%Lea Schuh 1,* , Peter V. Markov 1,2 , Ioanna Voulgaridi 3 , Zacharoula Bogogiannidou 3 , 
%Varvara A. Mouchtouri 3 , Christos Hadjichristodoulou3, Nikolaos I. Stilianakis1,4,*
%1 Joint Research Centre (JRC), European Commission, Ispra, Italy
%2 London School of Hygiene & Tropical Medicine, University of London, London, UK
%3 Laboratory of Hygiene and Epidemiology, University of Thessaly, Larissa, Greece
%4 Department of Biometry and Epidemiology, University of Erlangen-Nuremberg, Erlangen, Germany
%* corresponding authors

%% prep
addpath(genpath(pwd))
clearvars;
clc;

%% model fitting and parameter estimation for untreated and treated COVID-19 individuals 

%get fits for untreated individuals (model without treatment)
get_fits(0,'wo_treat',0)

%get fits for treated individuals (model with treatment)
get_fits(1,'w_treat',0)

%get fits for treated individuals (model without treatment - for comparison)
get_fits(0,'w_treat',0)

%% comparison analysis bteween non-hosp and hosp untreated, survived and deceased treated, and resp and non-resp individuals



