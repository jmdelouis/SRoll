#!/bin/sh

RELEASE=RD12_RC2

mkdir /pscratch1/delouis/${RELEASE}
mkdir /pscratch1/delouis/${RELEASE}/100
mkdir /pscratch1/delouis/${RELEASE}/143_SWB
mkdir /pscratch1/delouis/${RELEASE}/217_SWB
mkdir /pscratch1/delouis/${RELEASE}/353_SWB
mkdir /pscratch1/delouis/${RELEASE}/545
mkdir /pscratch1/delouis/${RELEASE}/857

mkdir /pscratch1/delouis/${RELEASE}/100/MAP
mkdir /pscratch1/delouis/${RELEASE}/143_SWB/MAP
mkdir /pscratch1/delouis/${RELEASE}/217_SWB/MAP
mkdir /pscratch1/delouis/${RELEASE}/353_SWB/MAP
mkdir /pscratch1/delouis/${RELEASE}/545/MAP
mkdir /pscratch1/delouis/${RELEASE}/857/MAP

mkdir /pscratch1/delouis/${RELEASE}/100/VECT
mkdir /pscratch1/delouis/${RELEASE}/143_SWB/VECT
mkdir /pscratch1/delouis/${RELEASE}/217_SWB/VECT
mkdir /pscratch1/delouis/${RELEASE}/353_SWB/VECT
mkdir /pscratch1/delouis/${RELEASE}/545/VECT
mkdir /pscratch1/delouis/${RELEASE}/857/VECT

mkdir ${RELEASE}
sed "s;fitang;"${RELEASE}";g" param100.txt > ${RELEASE}/param100.txt
sed "s;fitang;"${RELEASE}";g" param143.txt > ${RELEASE}/param143.txt
sed "s;fitang;"${RELEASE}";g" param143_noswb.txt > ${RELEASE}/param143_noswb.txt
sed "s;fitang;"${RELEASE}";g" param217.txt > ${RELEASE}/param217.txt
sed "s;fitang;"${RELEASE}";g" param217_noswb.txt > ${RELEASE}/param217_noswb.txt
sed "s;fitang;"${RELEASE}";g" param217_notf.txt > ${RELEASE}/param217_notf.txt
sed "s;fitang;"${RELEASE}";g" param353.txt > ${RELEASE}/param353.txt
sed "s;fitang;"${RELEASE}";g" param353_noswb.txt > ${RELEASE}/param353_noswb.txt
sed "s;fitang;"${RELEASE}";g" param545.txt > ${RELEASE}/param545.txt
sed "s;fitang;"${RELEASE}";g" param545_kcmb.txt > ${RELEASE}/param545_kcmb.txt
sed "s;fitang;"${RELEASE}";g" param857.txt > ${RELEASE}/param857.txt
sed "s;fitang;"${RELEASE}";g" do100.qsub > ${RELEASE}/do100.qsub
sed "s;fitang;"${RELEASE}";g" do143.qsub > ${RELEASE}/do143.qsub
sed  "s;_all;_noswb;g;s;fitang;"${RELEASE}";g;s;.txt;_noswb.txt;g;s;.log;_noswb.log;g" do143.qsub > ${RELEASE}/do143_noswb.qsub
sed "s;fitang;"${RELEASE}";g" do217.qsub > ${RELEASE}/do217.qsub
sed  "s;_all;_noswb;g;s;fitang;"${RELEASE}";g;s;.txt;_noswb.txt;g;s;.log;_noswb.log;g" do217.qsub > ${RELEASE}/do217_noswb.qsub
sed  "s;_all;_notf;g;s;fitang;"${RELEASE}";g;s;.txt;_notf.txt;g;s;.log;_notf.log;g" do217.qsub > ${RELEASE}/do217_notf.qsub
sed "s;fitang;"${RELEASE}";g" do353.qsub > ${RELEASE}/do353.qsub
sed  "s;_all;_noswb;g;s;fitang;"${RELEASE}";g;s;.txt;_noswb.txt;g;s;.log;_noswb.log;g" do353.qsub > ${RELEASE}/do353_noswb.qsub
sed "s;fitang;"${RELEASE}";g" do545.qsub > ${RELEASE}/do545.qsub
sed "s;fitang;"${RELEASE}";g" do857.qsub > ${RELEASE}/do857.qsub

#================================== ODD
sed "s;fitang;"${RELEASE}";g;s;v64_extended;v64_extended_Odd;g;s;RD12_REP6;RD12_REP6_Odd;g;s;REP6_RD12;REP6_RD12_Odd;g" param100_XX.txt > ${RELEASE}/param100_Odd.txt
sed "s;fitang;"${RELEASE}";g;s;v64_extended;v64_extended_Odd;g;s;RD12_REP6;RD12_REP6_Odd;g;s;REP6_RD12;REP6_RD12_Odd;g" param143_XX.txt > ${RELEASE}/param143_Odd.txt
sed "s;fitang;"${RELEASE}";g;s;v64_extended;v64_extended_Odd;g;s;RD12_REP6;RD12_REP6_Odd;g;s;REP6_RD12;REP6_RD12_Odd;g" param217_XX.txt > ${RELEASE}/param217_Odd.txt
sed "s;fitang;"${RELEASE}";g;s;v61ter;v61ter_Odd;g;s;RD12_REP6;RD12_REP6_Odd;g;s;REP6_RD12;REP6_RD12_Odd;g" param353_XX.txt > ${RELEASE}/param353_Odd.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;_v61;_v61_Odd;g;s;RD12_REP6;RD12_REP6_Odd;g;s;REP6_RD12;REP6_RD12_Odd;g" param545.txt > ${RELEASE}/param545_Odd.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;_v61;_v61_Odd;g;s;RD12_REP6;RD12_REP6_Odd;g;s;REP6_RD12;REP6_RD12_Odd;g" param545_kcmb.txt > ${RELEASE}/param545_kcmb_Odd.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;_v61;_v61_Odd;g;s;RD12_REP6;RD12_REP6_Odd;g;s;REP6_RD12;REP6_RD12_Odd;g" param857.txt > ${RELEASE}/param857_Odd.txt
sed "s;_all;_odd;g;s;fitang;"${RELEASE}";g;s;.txt;_Odd.txt;g;s;.log;_Odd.log;g" do100.qsub > ${RELEASE}/do100_Odd.qsub
sed "s;_all;_odd;g;s;fitang;"${RELEASE}";g;s;.txt;_Odd.txt;g;s;.log;_Odd.log;g" do143.qsub > ${RELEASE}/do143_Odd.qsub
sed "s;_all;_odd;g;s;fitang;"${RELEASE}";g;s;.txt;_Odd.txt;g;s;.log;_Odd.log;g" do217.qsub > ${RELEASE}/do217_Odd.qsub
sed "s;_all;_odd;g;s;fitang;"${RELEASE}";g;s;.txt;_Odd.txt;g;s;.log;_Odd.log;g" do353.qsub > ${RELEASE}/do353_Odd.qsub
sed "s;_all;_odd;g;s;fitang;"${RELEASE}";g;s;.txt;_Odd.txt;g;s;.log;_Odd.log;g" do545.qsub > ${RELEASE}/do545_Odd.qsub
sed "s;_all;_odd;g;s;fitang;"${RELEASE}";g;s;.txt;_Odd.txt;g;s;.log;_Odd.log;g" do857.qsub > ${RELEASE}/do857_Odd.qsub


#================================== EVEN
sed "s;fitang;"${RELEASE}";g;s;v64_extended;v64_extended_Even;g;s;RD12_REP6;RD12_REP6_Even;g;s;REP6_RD12;REP6_RD12_Even;g" param100_XX.txt > ${RELEASE}/param100_Even.txt
sed "s;fitang;"${RELEASE}";g;s;v64_extended;v64_extended_Even;g;s;RD12_REP6;RD12_REP6_Even;g;s;REP6_RD12;REP6_RD12_Even;g" param143_XX.txt > ${RELEASE}/param143_Even.txt
sed "s;fitang;"${RELEASE}";g;s;v64_extended;v64_extended_Even;g;s;RD12_REP6;RD12_REP6_Even;g;s;REP6_RD12;REP6_RD12_Even;g" param217_XX.txt > ${RELEASE}/param217_Even.txt
sed "s;fitang;"${RELEASE}";g;s;v61ter;v61ter_Even;g;s;RD12_REP6;RD12_REP6_Even;g;s;REP6_RD12;REP6_RD12_Even;g" param353_XX.txt > ${RELEASE}/param353_Even.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;_v61;_v61_Even;g;s;RD12_REP6;RD12_REP6_Even;g;s;REP6_RD12;REP6_RD12_Even;g" param545.txt > ${RELEASE}/param545_Even.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;_v61;_v61_Even;g;s;RD12_REP6;RD12_REP6_Even;g;s;REP6_RD12;REP6_RD12_Even;g" param545_kcmb.txt > ${RELEASE}/param545_kcmb_Even.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;_v61;_v61_Even;g;s;RD12_REP6;RD12_REP6_Even;g;s;REP6_RD12;REP6_RD12_Even;g" param857.txt > ${RELEASE}/param857_Even.txt
sed "s;_all;_even;g;s;fitang;"${RELEASE}";g;s;.txt;_Even.txt;g;s;.log;_Even.log;g" do100.qsub > ${RELEASE}/do100_Even.qsub
sed "s;_all;_even;g;s;fitang;"${RELEASE}";g;s;.txt;_Even.txt;g;s;.log;_Even.log;g" do143.qsub > ${RELEASE}/do143_Even.qsub
sed "s;_all;_even;g;s;fitang;"${RELEASE}";g;s;.txt;_Even.txt;g;s;.log;_Even.log;g" do217.qsub > ${RELEASE}/do217_Even.qsub
sed "s;_all;_even;g;s;fitang;"${RELEASE}";g;s;.txt;_Even.txt;g;s;.log;_Even.log;g" do353.qsub > ${RELEASE}/do353_Even.qsub
sed "s;_all;_even;g;s;fitang;"${RELEASE}";g;s;.txt;_Even.txt;g;s;.log;_Even.log;g" do545.qsub > ${RELEASE}/do545_Even.qsub
sed "s;_all;_even;g;s;fitang;"${RELEASE}";g;s;.txt;_Even.txt;g;s;.log;_Even.log;g" do857.qsub > ${RELEASE}/do857_Even.qsub

#================================== HM1
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 13144;g;s;RD12_REP6;RD12_REP6_HM1;g;s;REP6_RD12;REP6_RD12_HM1;g" param100_XX.txt > ${RELEASE}/param100_HM1.txt
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 13144;g;s;RD12_REP6;RD12_REP6_HM1;g;s;REP6_RD12;REP6_RD12_HM1;g" param143_XX.txt > ${RELEASE}/param143_HM1.txt
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 13144;g;s;RD12_REP6;RD12_REP6_HM1;g;s;REP6_RD12;REP6_RD12_HM1;g" param217_XX.txt > ${RELEASE}/param217_HM1.txt
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 13144;g;s;RD12_REP6;RD12_REP6_HM1;g;s;REP6_RD12;REP6_RD12_HM1;g" param353_XX.txt > ${RELEASE}/param353_HM1.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 13144;g;s;RD12_REP6;RD12_REP6_HM1;g;s;REP6_RD12;REP6_RD12_HM1;g" param545.txt > ${RELEASE}/param545_HM1.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 13144;g;s;RD12_REP6;RD12_REP6_HM1;g;s;REP6_RD12;REP6_RD12_HM1;g" param545_kcmb.txt > ${RELEASE}/param545_kcmb_HM1.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 13144;g;s;RD12_REP6;RD12_REP6_HM1;g;s;REP6_RD12;REP6_RD12_HM1;g" param857.txt > ${RELEASE}/param857_HM1.txt
sed "s;_all;_hm1;g;s;fitang;"${RELEASE}";g;s;.txt;_HM1.txt;g;s;.log;_HM1.log;g" do100.qsub > ${RELEASE}/do100_HM1.qsub
sed "s;_all;_hm1;g;s;fitang;"${RELEASE}";g;s;.txt;_HM1.txt;g;s;.log;_HM1.log;g" do143.qsub > ${RELEASE}/do143_HM1.qsub
sed "s;_all;_hm1;g;s;fitang;"${RELEASE}";g;s;.txt;_HM1.txt;g;s;.log;_HM1.log;g" do217.qsub > ${RELEASE}/do217_HM1.qsub
sed "s;_all;_hm1;g;s;fitang;"${RELEASE}";g;s;.txt;_HM1.txt;g;s;.log;_HM1.log;g" do353.qsub > ${RELEASE}/do353_HM1.qsub
sed "s;_all;_hm1;g;s;fitang;"${RELEASE}";g;s;.txt;_HM1.txt;g;s;.log;_HM1.log;g" do545.qsub > ${RELEASE}/do545_HM1.qsub
sed "s;_all;_hm1;g;s;fitang;"${RELEASE}";g;s;.txt;_HM1.txt;g;s;.log;_HM1.log;g" do857.qsub > ${RELEASE}/do857_HM1.qsub

#================================== HM2
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 13145;g;s;RD12_REP6;RD12_REP6_HM2;g;s;REP6_RD12;REP6_RD12_HM2;g" param100_XX.txt > ${RELEASE}/param100_HM2.txt
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 13145;g;s;RD12_REP6;RD12_REP6_HM2;g;s;REP6_RD12;REP6_RD12_HM2;g" param143_XX.txt > ${RELEASE}/param143_HM2.txt
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 13145;g;s;RD12_REP6;RD12_REP6_HM2;g;s;REP6_RD12;REP6_RD12_HM2;g" param217_XX.txt > ${RELEASE}/param217_HM2.txt
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 13145;g;s;RD12_REP6;RD12_REP6_HM2;g;s;REP6_RD12;REP6_RD12_HM2;g" param353_XX.txt > ${RELEASE}/param353_HM2.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 13145;g;s;RD12_REP6;RD12_REP6_HM2;g;s;REP6_RD12;REP6_RD12_HM2;g" param545.txt > ${RELEASE}/param545_HM2.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 13145;g;s;RD12_REP6;RD12_REP6_HM2;g;s;REP6_RD12;REP6_RD12_HM2;g" param545_kcmb.txt > ${RELEASE}/param545_kcmb_HM2.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 13145;g;s;RD12_REP6;RD12_REP6_HM2;g;s;REP6_RD12;REP6_RD12_HM2;g" param857.txt > ${RELEASE}/param857_HM2.txt
sed "s;_all;_hm2;g;s;fitang;"${RELEASE}";g;s;.txt;_HM2.txt;g;s;.log;_HM2.log;g" do100.qsub > ${RELEASE}/do100_HM2.qsub
sed "s;_all;_hm2;g;s;fitang;"${RELEASE}";g;s;.txt;_HM2.txt;g;s;.log;_HM2.log;g" do143.qsub > ${RELEASE}/do143_HM2.qsub
sed "s;_all;_hm2;g;s;fitang;"${RELEASE}";g;s;.txt;_HM2.txt;g;s;.log;_HM2.log;g" do217.qsub > ${RELEASE}/do217_HM2.qsub
sed "s;_all;_hm2;g;s;fitang;"${RELEASE}";g;s;.txt;_HM2.txt;g;s;.log;_HM2.log;g" do353.qsub > ${RELEASE}/do353_HM2.qsub
sed "s;_all;_hm2;g;s;fitang;"${RELEASE}";g;s;.txt;_HM2.txt;g;s;.log;_HM2.log;g" do545.qsub > ${RELEASE}/do545_HM2.qsub
sed "s;_all;_hm2;g;s;fitang;"${RELEASE}";g;s;.txt;_HM2.txt;g;s;.log;_HM2.log;g" do857.qsub > ${RELEASE}/do857_HM2.qsub


#================================== YR1
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 11194;g;s;RD12_REP6;RD12_REP6_YR1;g;s;REP6_RD12;REP6_RD12_YR1;g" param100_XX.txt > ${RELEASE}/param100_YR1.txt
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 11194;g;s;RD12_REP6;RD12_REP6_YR1;g;s;REP6_RD12;REP6_RD12_YR1;g" param143_XX.txt > ${RELEASE}/param143_YR1.txt
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 11194;g;s;RD12_REP6;RD12_REP6_YR1;g;s;REP6_RD12;REP6_RD12_YR1;g" param217_XX.txt > ${RELEASE}/param217_YR1.txt
sed "s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 11194;g;s;RD12_REP6;RD12_REP6_YR1;g;s;REP6_RD12;REP6_RD12_YR1;g" param353_XX.txt > ${RELEASE}/param353_YR1.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 11194;g;s;RD12_REP6;RD12_REP6_YR1;g;s;REP6_RD12;REP6_RD12_YR1;g" param545.txt > ${RELEASE}/param545_YR1.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 11194;g;s;RD12_REP6;RD12_REP6_YR1;g;s;REP6_RD12;REP6_RD12_YR1;g" param545_kcmb.txt > ${RELEASE}/param545_kcmb_YR1.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;EndRing   = 26050;EndRing   = 11194;g;s;RD12_REP6;RD12_REP6_YR1;g;s;REP6_RD12;REP6_RD12_YR1;g" param857.txt > ${RELEASE}/param857_YR1.txt
sed "s;_all;_yr1;g;s;fitang;"${RELEASE}";g;s;.txt;_YR1.txt;g;s;.log;_YR1.log;g" do100.qsub > ${RELEASE}/do100_YR1.qsub
sed "s;_all;_yr1;g;s;fitang;"${RELEASE}";g;s;.txt;_YR1.txt;g;s;.log;_YR1.log;g" do143.qsub > ${RELEASE}/do143_YR1.qsub
sed "s;_all;_yr1;g;s;fitang;"${RELEASE}";g;s;.txt;_YR1.txt;g;s;.log;_YR1.log;g" do217.qsub > ${RELEASE}/do217_YR1.qsub
sed "s;_all;_yr1;g;s;fitang;"${RELEASE}";g;s;.txt;_YR1.txt;g;s;.log;_YR1.log;g" do353.qsub > ${RELEASE}/do353_YR1.qsub
sed "s;_all;_yr1;g;s;fitang;"${RELEASE}";g;s;.txt;_YR1.txt;g;s;.log;_YR1.log;g" do545.qsub > ${RELEASE}/do545_YR1.qsub
sed "s;_all;_yr1;g;s;fitang;"${RELEASE}";g;s;.txt;_YR1.txt;g;s;.log;_YR1.log;g" do857.qsub > ${RELEASE}/do857_YR1.qsub

#================================== YR2
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 11195;g;s;RD12_REP6;RD12_REP6_YR2;g;s;REP6_RD12;REP6_RD12_YR2;g" param100_XX.txt > ${RELEASE}/param100_YR2.txt
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 11195;g;s;RD12_REP6;RD12_REP6_YR2;g;s;REP6_RD12;REP6_RD12_YR2;g" param143_XX.txt > ${RELEASE}/param143_YR2.txt
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 11195;g;s;RD12_REP6;RD12_REP6_YR2;g;s;REP6_RD12;REP6_RD12_YR2;g" param217_XX.txt > ${RELEASE}/param217_YR2.txt
sed "s;fitang;"${RELEASE}";g;s;BeginRing = 240;BeginRing = 11195;g;s;RD12_REP6;RD12_REP6_YR2;g;s;REP6_RD12;REP6_RD12_YR2;g" param353_XX.txt > ${RELEASE}/param353_YR2.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;BeginRing   = 240;BeginRing   = 11195;g;s;RD12_REP6;RD12_REP6_YR2;g;s;REP6_RD12;REP6_RD12_YR2;g" param545.txt > ${RELEASE}/param545_YR2.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;BeginRing   = 240;BeginRing   = 11195;g;s;RD12_REP6;RD12_REP6_YR2;g;s;REP6_RD12;REP6_RD12_YR2;g" param545_kcmb.txt > ${RELEASE}/param545_kcmb_YR2.txt
sed "s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;BeginRing   = 240;BeginRing   = 11195;g;s;RD12_REP6;RD12_REP6_YR2;g;s;REP6_RD12;REP6_RD12_YR2;g" param857.txt > ${RELEASE}/param857_YR2.txt
sed "s;_all;_yr2;g;s;fitang;"${RELEASE}";g;s;.txt;_YR2.txt;g;s;.log;_YR2.log;g" do100.qsub > ${RELEASE}/do100_YR2.qsub
sed "s;_all;_yr2;g;s;fitang;"${RELEASE}";g;s;.txt;_YR2.txt;g;s;.log;_YR2.log;g" do143.qsub > ${RELEASE}/do143_YR2.qsub
sed "s;_all;_yr2;g;s;fitang;"${RELEASE}";g;s;.txt;_YR2.txt;g;s;.log;_YR2.log;g" do217.qsub > ${RELEASE}/do217_YR2.qsub
sed "s;_all;_yr2;g;s;fitang;"${RELEASE}";g;s;.txt;_YR2.txt;g;s;.log;_YR2.log;g" do353.qsub > ${RELEASE}/do353_YR2.qsub
sed "s;_all;_yr2;g;s;fitang;"${RELEASE}";g;s;.txt;_YR2.txt;g;s;.log;_YR2.log;g" do545.qsub > ${RELEASE}/do545_YR2.qsub
sed "s;_all;_yr2;g;s;fitang;"${RELEASE}";g;s;.txt;_YR2.txt;g;s;.log;_YR2.log;g" do857.qsub > ${RELEASE}/do857_YR2.qsub

#================================== HR1
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1;g;s;REP6_RD12;REP6_RD12_HR1;g" param100_XX.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param100_HR1.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1;g;s;REP6_RD12;REP6_RD12_HR1;g" param143_XX.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param143_HR1.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1;g;s;REP6_RD12;REP6_RD12_HR1;g" param217_XX.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param217_HR1.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1;g;s;REP6_RD12;REP6_RD12_HR1;g" param353_XX.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param353_HR1.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1_NOSWD;g;s;REP6_RD12;REP6_RD12_HR1_NOSWD;g" param143_XX_noswb.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param143_HR1_noswb.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1_NOSWD;g;s;REP6_RD12;REP6_RD12_HR1_NOSWD;g" param217_XX_noswb.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param217_HR1_noswb.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1_NOSWD;g;s;REP6_RD12;REP6_RD12_HR1_NOSWD;g" param353_XX_noswb.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param353_HR1_noswb.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1;g;s;REP6_RD12;REP6_RD12_HR1;g" param545.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param545_HR1.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1;g;s;REP6_RD12;REP6_RD12_HR1;g" param545_kcmb.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param545_kcmb_HR1.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR1;g;s;REP6_RD12;REP6_RD12_HR1;g" param857.txt | sed "s;REP6;REP6_1st;g" > ${RELEASE}/param857_HR1.txt
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1.txt;g;s;.log;_HR1.log;g" do100.qsub > ${RELEASE}/do100_HR1.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1.txt;g;s;.log;_HR1.log;g" do143.qsub > ${RELEASE}/do143_HR1.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1.txt;g;s;.log;_HR1.log;g" do217.qsub > ${RELEASE}/do217_HR1.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1.txt;g;s;.log;_HR1.log;g" do353.qsub > ${RELEASE}/do353_HR1.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1_noswb.txt;g;s;.log;_HR1_noswb.log;g" do143.qsub > ${RELEASE}/do143_HR1_noswb.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1_noswb.txt;g;s;.log;_HR1_noswb.log;g" do217.qsub > ${RELEASE}/do217_HR1_noswb.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1_noswb.txt;g;s;.log;_HR1_noswb.log;g" do353.qsub > ${RELEASE}/do353_HR1_noswb.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1.txt;g;s;.log;_HR1.log;g" do545.qsub > ${RELEASE}/do545_HR1.qsub
sed "s;_all;_hr1;g;s;fitang;"${RELEASE}";g;s;.txt;_HR1.txt;g;s;.log;_HR1.log;g" do857.qsub > ${RELEASE}/do857_HR1.qsub

#================================== HR2
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2;g;s;REP6_RD12;REP6_RD12_HR2;g" param100_XX.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param100_HR2.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2;g;s;REP6_RD12;REP6_RD12_HR2;g" param143_XX.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param143_HR2.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2;g;s;REP6_RD12;REP6_RD12_HR2;g" param217_XX.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param217_HR2.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2;g;s;REP6_RD12;REP6_RD12_HR2;g" param353_XX.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param353_HR2.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2_NOSWD;g;s;REP6_RD12;REP6_RD12_HR2_NOSWD;g" param143_XX_noswb.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param143_HR2_noswb.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2_NOSWD;g;s;REP6_RD12;REP6_RD12_HR2_NOSWD;g" param217_XX_noswb.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param217_HR2_noswb.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2_NOSWD;g;s;REP6_RD12;REP6_RD12_HR2_NOSWD;g" param353_XX_noswb.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param353_HR2_noswb.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2;g;s;REP6_RD12;REP6_RD12_HR2;g" param545.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param545_HR2.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2;g;s;REP6_RD12;REP6_RD12_HR2;g" param545_kcmb.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param545_kcmb_HR2.txt
sed "s;adutot;adu;g;s;DOMAXVRAIE = 0;DOMAXVRAIE = 1;g;s;fitang;"${RELEASE}";g;s;RD12_REP6;RD12_REP6_HR2;g;s;REP6_RD12;REP6_RD12_HR2;g" param857.txt | sed "s;REP6;REP6_2nd;g" > ${RELEASE}/param857_HR2.txt
sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2.txt;g;s;.log;_HR2.log;g" do100.qsub > ${RELEASE}/do100_HR2.qsub
sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2.txt;g;s;.log;_HR2.log;g" do143.qsub > ${RELEASE}/do143_HR2.qsub
sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2.txt;g;s;.log;_HR2.log;g" do217.qsub > ${RELEASE}/do217_HR2.qsub
sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2.txt;g;s;.log;_HR2.log;g" do353.qsub > ${RELEASE}/do353_HR2.qsub

sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2_noswb.txt;g;s;.log;_HR2_noswb.log;g" do143.qsub > ${RELEASE}/do143_HR2_noswb.qsub
sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2_noswb.txt;g;s;.log;_HR2_noswb.log;g" do217.qsub > ${RELEASE}/do217_HR2_noswb.qsub
sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2_noswb.txt;g;s;.log;_HR2_noswb.log;g" do353.qsub > ${RELEASE}/do353_HR2_noswb.qsub

sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2.txt;g;s;.log;_HR2.log;g" do545.qsub > ${RELEASE}/do545_HR2.qsub
sed "s;_all;_hr2;g;s;fitang;"${RELEASE}";g;s;.txt;_HR2.txt;g;s;.log;_HR2.log;g" do857.qsub > ${RELEASE}/do857_HR2.qsub


#================================== DS1
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS1;g" param100_ds1.txt > ${RELEASE}/param100_ds1.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS1;g" param143_ds1.txt > ${RELEASE}/param143_ds1.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS1;g" param217_ds1.txt > ${RELEASE}/param217_ds1.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS1;g" param353_ds1.txt > ${RELEASE}/param353_ds1.txt
sed "s;_all;_ds1;g;s;fitang;"${RELEASE}";g;s;.txt;_ds1.txt;g;s;.log;_ds1.log;g" do100.qsub > ${RELEASE}/do100_ds1.qsub
sed "s;_all;_ds1;g;s;fitang;"${RELEASE}";g;s;.txt;_ds1.txt;g;s;.log;_ds1.log;g" do143.qsub > ${RELEASE}/do143_ds1.qsub
sed "s;_all;_ds1;g;s;fitang;"${RELEASE}";g;s;.txt;_ds1.txt;g;s;.log;_ds1.log;g" do217.qsub > ${RELEASE}/do217_ds1.qsub
sed "s;_all;_ds1;g;s;fitang;"${RELEASE}";g;s;.txt;_ds1.txt;g;s;.log;_ds1.log;g" do353.qsub > ${RELEASE}/do353_ds1.qsub


#================================== DS2
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS2;g" param100_ds2.txt > ${RELEASE}/param100_ds2.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS2;g" param143_ds2.txt > ${RELEASE}/param143_ds2.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS2;g" param217_ds2.txt > ${RELEASE}/param217_ds2.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS2;g" param353_ds2.txt > ${RELEASE}/param353_ds2.txt
sed "s;_all;_ds2;g;s;fitang;"${RELEASE}";g;s;.txt;_ds2.txt;g;s;.log;_ds2.log;g" do100.qsub > ${RELEASE}/do100_ds2.qsub
sed "s;_all;_ds2;g;s;fitang;"${RELEASE}";g;s;.txt;_ds2.txt;g;s;.log;_ds2.log;g" do143.qsub > ${RELEASE}/do143_ds2.qsub
sed "s;_all;_ds2;g;s;fitang;"${RELEASE}";g;s;.txt;_ds2.txt;g;s;.log;_ds2.log;g" do217.qsub > ${RELEASE}/do217_ds2.qsub
sed "s;_all;_ds2;g;s;fitang;"${RELEASE}";g;s;.txt;_ds2.txt;g;s;.log;_ds2.log;g" do353.qsub > ${RELEASE}/do353_ds2.qsub

#================================== DS3
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS3;g" param143_ds3.txt > ${RELEASE}/param143_ds3.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS3;g" param217_ds3.txt > ${RELEASE}/param217_ds3.txt
sed "s;fitang;"${RELEASE}";g;s;REP6_RD12;REP6_RD12_DS3;g" param353_ds3.txt > ${RELEASE}/param353_ds3.txt
sed "s;_all;_ds3;g;s;fitang;"${RELEASE}";g;s;.txt;_ds3.txt;g;s;.log;_ds3.log;g" do143.qsub > ${RELEASE}/do143_ds3.qsub
sed "s;_all;_ds3;g;s;fitang;"${RELEASE}";g;s;.txt;_ds3.txt;g;s;.log;_ds3.log;g" do217.qsub > ${RELEASE}/do217_ds3.qsub
sed "s;_all;_ds3;g;s;fitang;"${RELEASE}";g;s;.txt;_ds3.txt;g;s;.log;_ds3.log;g" do353.qsub > ${RELEASE}/do353_ds3.qsub
