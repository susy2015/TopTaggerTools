#ifndef WEIGHT_H
#define WEIGHT_H

class Weight{
 public:

  static float MuTriWgt_2016(int ptbin);
  static float MuTriWgt_2017(int ptbin);
  static float MuTriWgt_2018(int ptbin);
  static float MuTriWgt(int etabin);
  static float HtTriWgt(int metbin, int htbin);
  static float TopEffSF_2016(int topptbin);
  static float TopMisSF_2016(int topptbin);
  static float TopSF_2017(int topptbin);
  static float TopSF_2018(int topptbin);
  static float TopEffSF_2016(int topptbin, int jetbin);
  static float TopMisSF_2016(int topptbin, int jetbin);

  static int ptbin(float pt);
  static int etabin(float eta);
  static int metbin(float met);
  static int htbin(float ht);
  static int Topptbin(float toppt);
  static int jetbin(int jet);

  static float TopEffSF_16;
  static float TopMisSF_16;

  static float TopEffSF_17;
  static float TopMisSF_17;

  static float TopEffSF_18;
  static float TopMisSF_18;

};

float Weight::TopEffSF_16 = 1.01964;
float Weight::TopMisSF_16 = 0.907106;

float Weight::TopEffSF_17 = 0.957049;
float Weight::TopMisSF_17 = 0.992094;

float Weight::TopEffSF_18 = 0.935941;
float Weight::TopMisSF_18 = 0.903555;

int Weight::ptbin(float muPt){
  int bin =0;
  if(10 <= muPt && muPt < 20) bin =1;
  if(20 <= muPt && muPt < 30) bin =2;
  if(30 <= muPt && muPt < 40) bin =3;
  if(40 <= muPt && muPt < 50) bin =4;
  if(50 <= muPt && muPt < 60) bin =5;
  if(60 <= muPt && muPt < 70) bin =6;
  if(70 <= muPt && muPt < 80) bin =7;
  if(80 <= muPt && muPt < 90) bin =8;
  if(90 <= muPt && muPt < 100) bin =9;
  if(100 <= muPt && muPt < 110) bin =10;
  if(110 <= muPt && muPt < 120) bin =11;
  if(120 <= muPt && muPt < 130) bin =12;
  if(130 <= muPt && muPt < 140) bin =13;
  if(140 <= muPt && muPt < 150) bin =14;
  if(150 <= muPt && muPt < 200) bin =15;
  if(200 <= muPt && muPt < 250) bin =16;
  if(250 <= muPt && muPt < 300) bin =17;
  if(300 <= muPt && muPt < 350) bin =18;
  if(350 <= muPt && muPt < 400) bin =19;
  if(400 <= muPt)               bin =20;
  return bin;
}

int Weight::etabin(float muEta){
  int bin =0;
  if(-2.2 <= muEta && muEta < -1.8) bin =1;
  if(-1.8 <= muEta && muEta < -1.4) bin =2;
  if(-1.4 <= muEta && muEta < -1.0) bin =3;
  if(-1.0 <= muEta && muEta < -0.6) bin =4;
  if(-0.6 <= muEta && muEta < -0.2) bin=5;
  if(-0.2 <= muEta && muEta < 0.2) bin=6;
  if( 0.2 <= muEta && muEta < 0.6) bin=7;
  if( 0.6 <= muEta && muEta < 1.0) bin=8;
  if( 1.0 <= muEta && muEta < 1.4) bin=9;
  if( 1.4 <= muEta && muEta < 1.8) bin=10;
  if( 1.8 <= muEta && muEta < 2.2) bin=11;
  if( 2.2 <= muEta && muEta < 2.6) bin=12;
  return bin;
}

int Weight::htbin(float ht){

  int bin =0;
  if(ht>=1000) bin =1;
  return bin;
}

int Weight::metbin(float met){
  int bin =0;
  if(met>=25 && met<50) bin=1;
  if(met>=50 && met<75) bin=2;
  if(met>=75 && met<100) bin=3;
  if(met>=100 && met<125) bin=4;
  if(met>=125 && met<150) bin=5;
  if(met>=150 && met<175) bin=6;
  if(met>=175 && met<200) bin=7;
  if(met>=200 && met<275) bin=8;
  if(met>=275 && met<400) bin=9;
  if(met>=400 && met<600) bin=10;
  if(met>=600 && met<1000) bin=11;
  if(met>=1000) bin=12;
  return bin;
}

int Weight::Topptbin(float topPt){
  int bin =0;
  if(250 <= topPt && topPt < 300) bin =1;
  if(300 <= topPt && topPt < 350) bin =2;
  if(350 <= topPt && topPt < 450) bin =3;
  if(450 <= topPt)                bin =4;

  return bin;
}

int Weight::jetbin(int jet){
  int bin =0;
  if(5 <= jet && jet < 6) bin =1;
  if(6 <= jet && jet < 7) bin =2;
  if(7 <= jet && jet < 8) bin =3;
  if(8 <= jet)            bin =4;

  return bin;
}

float Weight::MuTriWgt_2016(int MuPtBin){
  float TriWgt[21] = {0.0516629 ,0.0219888 ,0.630774 ,0.801294 ,0.824612 ,0.906815 ,0.911391 ,0.913463 ,0.916105 ,0.916083 ,0.917895 ,0.918171 ,0.917234 ,0.915588 ,0.914291 ,0.912883 ,0.904298 ,0.899571 ,0.902224 ,0.897059 ,0.884043};
  float Wgt = TriWgt[MuPtBin];
  return Wgt;
}

float Weight::MuTriWgt_2017(int MuPtBin){
  float TriWgt[21] = {0.0169829 ,0.0144944 ,0.393525 ,0.755547 ,0.791661 ,0.889515 ,0.898893 ,0.90139 ,0.90033 ,0.90281 ,0.903197 ,0.902215 ,0.902212 ,0.898964 ,0.898811 ,0.895749 ,0.88714 ,0.885658 ,0.879651 ,0.878254 ,0.873968};
  float Wgt = TriWgt[MuPtBin];
  return Wgt;
}

float Weight::MuTriWgt_2018(int MuPtBin){
  float TriWgt[21] = {0.0331628 ,0.0212502 ,0.425248 ,0.775295 ,0.807893 ,0.916184 ,0.918452 ,0.921058 ,0.92099 ,0.92114 ,0.921068 ,0.919682 ,0.919659 ,0.918916 ,0.917512 ,0.914624 ,0.90798 ,0.906 ,0.900144 ,0.900776 ,0.886938};
  float Wgt = TriWgt[MuPtBin];
  return Wgt;
}

float Weight::MuTriWgt(int MuEtaBin){
  float TriWgt[13] = {0.8020833, 0.8113949, 0.8111837, 0.8824405, 0.9024091, 0.8737864, 0.9186085, 0.8759649, 0.8940410, 0.8848286, 0.8293217, 0.8263979, 0.7605634};
  float Wgt = TriWgt[MuEtaBin];
  return Wgt;
}

float Weight::HtTriWgt(int metbin, int htbin){

  float TriWgt[2][13]= {{0.001542561, 0.003222389, 0.00987073, 0.03865682, 0.1387231, 0.3564816, 0.6276442, 0.8154821, 0.9340538, 0.9858562, 0.9931507, 1.0, 1.0},{0.02067183, 0.02504944, 0.04486466, 0.07434402, 0.1518288, 0.2802669, 0.4642409, 0.6596434, 0.8510453, 0.9563492, 0.9874214, 0.9736842, 0.9736842}};
  float Wgt = TriWgt[htbin][metbin];
  return Wgt;
}

float Weight::TopEffSF_2016(int topptbin){
  float SF[9] = {1.01964, 1.01964, 1.01964, 1.01964, 1.01964, 1.01964, 1.01964, 1.01964, 1.01964};
  float sf = SF[topptbin];
  return sf;
}

float Weight::TopMisSF_2016(int topptbin){
  float SF[9] = {0.645982 ,0.827829 ,0.898575 ,0.90248 ,0.933803 ,0.965716 ,1.0716 ,1.12506 ,1.21877};
  float sf = SF[topptbin];
  return sf;
}

float Weight::TopEffSF_2016(int topptbin, int jetbin){
  float SF[5][5] = {{1.39735, 1.01588, 1.1131, 1.02869, 1.13937},{1.14548, 1.00906, 0.866989, 1.28356, 0.889261},{1.04679, 0.999998, 1.01572, 1.0111, 1.18435},{1.11069, 1.03944, 0.914347, 0.941984, 0.939497},{1.05261, 0.850017, 1.05786, 1.01094, 1.22246}};
  float sf = SF[topptbin][jetbin];
  return sf;
}

float Weight::TopMisSF_2016(int topptbin, int jetbin){
  float SF[5][5] = {{0.675671, 0.709362, 0.71803, 0.787825, 0.836342},{0.911494, 0.774873, 0.863443, 0.987031, 0.999224},{0.805839, 0.802888, 0.863861, 1.01889, 1.0865},{0.853396, 0.815337, 0.917457, 1.06121, 1.23127},{1.01642, 1.02978, 1.13235, 1.31968, 1.75918}};
  float sf = SF[topptbin][jetbin];
  return sf;
}

#endif
