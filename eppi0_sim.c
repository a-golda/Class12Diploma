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
#include <algorithm>
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

out<< "px" << ',' << "py" << ',' << "pz" << ',' << "E" << ',' << "Theta" << ',' << "Phi" << ',' << "E_e" 
<<','<< "e_theta" <<',' << "e_phi"<< ','<< "E_p" <<','<< "p_theta" <<',' << "p_phi"<< ','<< "E_cal_max" << ',' << "E_cal_min" << ','<< "E_cal_sum" <<','<< "Cal_n" << ',' << "is_good" <<std::endl;

auto start = std::chrono::high_resolution_clock::now();

gBenchmark->Start("timer");
  
std::cout << " reading file example program (HIPO) " << std::endl;

vector<string> data;	

data.push_back("out_out1.hipo");
data.push_back("out_out2.hipo");
data.push_back("out_out3.hipo");
data.push_back("out_out4.hipo");
data.push_back("out_out5.hipo");
data.push_back("out_out6.hipo");
data.push_back("out_out7.hipo");
data.push_back("out_out8.hipo");
data.push_back("out_out9.hipo");
data.push_back("out_out10.hipo");

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
    int nrowscal = ECAL.getRows();

    vector<int> proton;  //вектора с номерами, соответствующие каждому из типов частиц
    vector<int> gamma;
    vector<int> electron;

    vector<int> proton_id; 
    vector<int> gamma_id;
    vector<int> electron_id;


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

    TLorentzVector Egamtest1(0,0,0,0);
    TLorentzVector Temp1(0,0,0,0); 
    TLorentzVector Temp2(0,0,0,0); //некие вспомогательные вектора, с помощью которых я сравниваю
    vector<float> g1mc_simi; 
    vector<float> g2mc_simi;
    vector<float> Cal_E;

    float perc = 0.1;
    
    if(electron.size()==1 && proton.size()==1 && gamma.size()>=1)
    {
      pr.SetXYZM(PART.getFloat("px",proton[0]),PART.getFloat("py",proton[0]),PART.getFloat("pz",proton[0]),0.938);
      
      for(int k=0; k<gamma.size(); k++)
      { 
        Egamtest1.SetXYZM(PART.getFloat("px",gamma[k]),PART.getFloat("py",gamma[k]),PART.getFloat("pz",gamma[k]),0);
        if (((abs(Egamtest1.E()-g1mc.E())/g1mc.E()) < perc) && (abs(Egamtest1.Theta()-g1mc.Theta())/g1mc.Theta()<perc) && (abs(Egamtest1.Phi()-g1mc.Phi())/g1mc.Phi()<perc))
        {
          g1mc_simi.push_back(k);
        } else {}
        if ((abs(Egamtest1.E()-g2mc.E())/g2mc.E()) < perc && (abs(Egamtest1.Theta()-g2mc.Theta())/g2mc.Theta()<perc) && (abs(Egamtest1.Phi()-g2mc.Phi())/g2mc.Phi()<perc))
        {
          g1mc_simi.push_back(k);
        } else {}

        for(int r_cal=0; r_cal<nrowscal; r_cal++)
        {
          if (k==ECAL.getInt("pindex", r_cal))
          {
            Cal_E.push_back(ECAL.getFloat("energy", r_cal));
          } else {}
        }

      }
      

      for(int m=0; m<gamma.size(); m++)
      {
        int Cal_n = Cal_E.size();
        float E_cal_max = *std::max_element(Cal_E.begin(), Cal_E.end());
        float E_cal_min = *std::min_element(Cal_E.begin(), Cal_E.end());
        float E_cal_sum = 0;
        for (int s=0; s<Cal_E.size(); s++)
        {
          E_cal_sum+=Cal_E[s];
        }

        Temp1.SetXYZM(PART.getFloat("px",gamma[m]),PART.getFloat("py",gamma[m]),PART.getFloat("pz",gamma[m]),0);
        if(std::find(g1mc_simi.begin(), g1mc_simi.end(), m) != g1mc_simi.end()) // like 'in' func
        {
          out<< PART.getFloat("px",gamma[m]) << ',' << PART.getFloat("py",gamma[m])<< ',' << PART.getFloat("pz",gamma[m])  << ',' << Temp1.E() << ',' << Temp1.Theta()<< ','<< Temp1.Phi()<<','<< el.E() <<','<< el.Theta()<<',' << el.Phi()<<',' << pr.E() <<','<< pr.Theta()<<',' << pr.Phi() <<','<< E_cal_max << ',' << E_cal_min << ','<< E_cal_sum <<','<< Cal_n << ',' <<'1'<<std::endl;
        }
        else
        {
          out<< PART.getFloat("px",gamma[m]) << ',' << PART.getFloat("py",gamma[m])<< ',' << PART.getFloat("pz",gamma[m])  << ',' << Temp1.E() << ',' << Temp1.Theta()<< ','<< Temp1.Phi()<<','<< el.E() <<','<< el.Theta()<<',' << el.Phi()<<',' << pr.E() <<','<< pr.Theta()<<',' << pr.Phi() <<','<< E_cal_max << ',' << E_cal_min << ','<< E_cal_sum <<','<< Cal_n << ',' << '0'<<std::endl;
        }
      }

    } else {}
  } 
}

  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

  out.close();
}
