#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"

using namespace clas12;

void eppi0_sim()
{  
auto db=TDatabasePDG::Instance();

TLorentzVector beam(0,0,6.535,6.535); //Объявление четырехвекторов, описывающих кинематику
TLorentzVector el;
TLorentzVector pr;
TLorentzVector g1;
TLorentzVector g2;
TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());  
TLorentzVector q;
TLorentzVector W;
TLorentzVector pi0;
TLorentzVector missPi0; 
double angle_g1_g2;
//для генерированных событий
TLorentzVector elmc;
TLorentzVector pi0mc;
TLorentzVector g1mc;
TLorentzVector g2mc;
TLorentzVector gbigmc;
TLorentzVector prmc;
TLorentzVector qmc;
TLorentzVector Wmc;

gBenchmark->Start("timer");
int counter=0;

std::ofstream out;
out.open("/home/golodovka/Programm/Clas12Tool/RunRoot/1pack/result.txt");
vector<int> result; 

auto start = std::chrono::high_resolution_clock::now();

gBenchmark->Start("timer");
  
std::cout << " reading file example program (HIPO) " << std::endl;

vector<string> data;	

data.push_back("out_out1.hipo");
// data.push_back("out_out2.hipo");
// data.push_back("out_out3.hipo");
// data.push_back("out_out4.hipo");
// data.push_back("out_out5.hipo");
// data.push_back("out_out6.hipo");
// data.push_back("out_out7.hipo");
// data.push_back("out_out8.hipo");
// data.push_back("out_out9.hipo");
// data.push_back("out_out10.hipo");

for(int r=0;r<data.size();r++)
{ 
  hipo::reader  reader;
  reader.open(data[r].c_str());

  hipo::dictionary factory;
  reader.readDictionary(factory);

  hipo::event event; 
  hipo::bank PART(factory.getSchema("REC::Particle"));
  hipo::bank PARTMC(factory.getSchema("MC::Lund"));
  hipo::bank ECAL(factory.getSchema("REC::Calorimeter"));

  while(reader.next()==true)
  {
    counter++;
    if (counter==1) continue;
    if (counter % 100000 == 0) cout<<counter<<"\n"; 

    reader.read(event);
    event.getStructure(PART);
    event.getStructure(ECAL);
    event.getStructure(PARTMC);

    int nrows = PART.getRows();  // количество частиц в рассматриваемом событии
    int nrowsmc = PARTMC.getRows();

    vector<int> proton;  //вектора с номерами, соответствующие каждому из типов частиц
    vector<int> gamma;
    vector<int> electron;

    int secte=-10;  //переменные для сравнивания секторов, в которых регистрировались частицы
    int sectp=-10;
    int sectg=-10;

    //generated events
    elmc.SetXYZM(PARTMC.getFloat("px",0),PARTMC.getFloat("py",0),PARTMC.getFloat("pz",0),0);
    qmc=beam-elmc;
    Wmc=beam-elmc+target;

    int ngammc=0;
    for (int i = 1; i < nrowsmc; i++)
    {
      int pidmc = PARTMC.getInt("pid",i);
      if(pidmc==2212) prmc.SetXYZM(PARTMC.getFloat("px",i),PARTMC.getFloat("py",i),PARTMC.getFloat("pz",i),0.9383); 
      if(pidmc==22)
      {
        if(ngammc==0)
        {
          g1mc.SetXYZM(PARTMC.getFloat("px",i),PARTMC.getFloat("py",i),PARTMC.getFloat("pz",i),0.);
          ngammc++;
        }else  g2mc.SetXYZM(PARTMC.getFloat("px",i),PARTMC.getFloat("py",i),PARTMC.getFloat("pz",i),0.);
      }
    }

    if(g2mc.E()<g1mc.E())
    { gbigmc=g1mc; 
      g1mc=g2mc; 
      g2mc=gbigmc;
    }
    pi0mc=g1mc+g2mc;

    //reconstructed events
    el.SetXYZM(PART.getFloat("px",0),PART.getFloat("py",0),PART.getFloat("pz",0),0);
    q=beam-el;
    W=beam-el+target;

    //заполнение векторов с номерами частиц

    for (int i = 0; i < nrows; i++)
    {
      int pid = PART.getInt("pid",i);
      if(pid==11) electron.push_back(i);
      if(pid==2212) proton.push_back(i);
      if(pid==22) gamma.push_back(i);
    } 

    TLorentzVector Egamtest1(0,0,0,0);  //некие вспомогательные вектора, с помощью которых я сравниваю
    vector<float> g1mc_simi; 
    vector<float> g2mc_simi;
    float perc = 0.1;
    
    if(electron.size()==1 && proton.size()==1 && gamma.size()>=1)
    {
      for(int k=0; k<gamma.size(); k++)
      { 
        Egamtest1.SetXYZM(PART.getFloat("px",gamma[k]),PART.getFloat("py",gamma[k]),PART.getFloat("pz",gamma[k]),0);
        if (((abs(Egamtest1.E()-g1mc.E())/g1mc.E()) < perc) && (Egamtest1.Theta()<0.1))
        {
          g1mc_simi.push_back(k);
        } else {}
        if ((abs(Egamtest1.E()-g2mc.E())/g2mc.E()) < perc)
        {
          g2mc_simi.push_back(k);
        } else {}
      }
    } else {}

    if(g1mc_simi.size()!=0)
    {
      for(int m=0; m<g1mc_simi.size(); m++)
      {
        out<<g1mc.Theta()<<','<<g2mc.Phi()<<std::endl;
        out<< PART.getFloat("px",gamma[m]) << ',' << PART.getFloat("py",gamma[m])<< ',' << PART.getFloat("pz",gamma[m])<< ',' << std::endl;
      }
    } else {}

    if(g2mc_simi.size()!=0)
    { 
      for(int l=0; l<g2mc_simi.size(); l++)
      {
        out<< PART.getFloat("px",gamma[l]) << ',' << PART.getFloat("py",gamma[l])<< ',' << PART.getFloat("pz",gamma[l])<< ',' << std::endl;
      }
    } else {}

    
    // if g1mc_simi.size() != 0
    // {
    // }
    // for (int i = 0; i < ; i++)
    // {
    //     if (words[i] == contra)
    //     {
    //         result = i;
    //         break;
    //     }
    // }

    // [0.1,0.2,0.3] g1mc g2mc
    // [0.2,0.4,0.1] g1 g2 g3 g4 ..... 

    // if(abs(g1mc.E() - g1.E())<0.1 && abs(g2mc.E() - g2.E())<0.1)
    // {
    //   result.push_back(1);
    // } else{
    //         result.push_back(0);  
    //       }
    // if (out.is_open())
    //   {   
    //     out<< g1mc.E()<<','<< g2mc.E()<< ',' << result[result.size() - 1] << std::endl;
    //   } 
  } 
}

  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

  out.close();
  /// s gitom vse okey i prom
}
