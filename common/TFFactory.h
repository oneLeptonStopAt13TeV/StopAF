#ifndef TFFactory_h
#define TFFactory_h
// C/C++ headerrs
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <exception>

//SS headers
#include "../sonicScrewdriver/interface/Table.h"
#include "../sonicScrewdriver/interface/Figure.h"

using namespace std;
using namespace theDoctor;

struct TFnames
{
    string CR;
    string SR;
    string data;
    string MC;
};

class TFFactory
{
public:

TFnames TFmem;
string TFRegion;

private:

Table tab;
Figure dataCR;
Figure MCSR;
Figure MCCR;
//TODO
double additionalError;

public:

TFFactory(string inputTable, string region = "standard", double error = 0)
{
    tab = Table(inputTable);
    dataCR = Figure(0,0);
    MCSR = Figure(0,0);
    MCCR = Figure(0,0);
    additionalError = error;
    TFmem.CR = "CR2l";
    TFmem.SR = "SR1l";
    TFmem.data = "data";
    TFmem.MC = "totalSM";
    TFRegion = region;
}

void defineFactors()
{
    if(TFRegion == "standard")
        return;
    else if(TFRegion == "2jMET250to350")
    {
        TFmem.CR = "CR2l2jMET250to350";
        TFmem.SR = "SR1l2jMET250to350";
    }
    else if(TFRegion == "2jMET350to450")
    {
        TFmem.CR = "CR2l2jMET350to450";
        TFmem.SR = "SR1l2jMET350to450";
    }
    else if(TFRegion == "2jMET450toInf")
    {
        TFmem.CR = "CR2l2jMET450toInf";
        TFmem.SR = "SR1l2jMET450toInf";
    }
    else if(TFRegion == "3jMET250to350")
    {
        TFmem.CR = "CR2l3jMET250to350";
        TFmem.SR = "SR1l3jMET250to350";
    }
    else if(TFRegion == "3jMET350to450")
    {
        TFmem.CR = "CR2l3jMET350to450";
        TFmem.SR = "SR1l3jMET350to450";
    }
    else if(TFRegion == "3jMET450to550")
    {
        TFmem.CR = "CR2l3jMET450to550";
        TFmem.SR = "SR1l3jMET450to550";
    }
    else if(TFRegion == "3jMET550toInf")
    {
        TFmem.CR = "CR2l3jMET550toInf";
        TFmem.SR = "SR1l3jMET550toInf";
    }
    else if(TFRegion == "4jMET250to350lowMT2W")
    {
        TFmem.CR = "CR2l4jMET250to350lowMT2W";
        TFmem.SR = "SR1l4jMET250to350lowMT2W";
    }
    else if(TFRegion == "4jMET350to450lowMT2W")
    {
        TFmem.CR = "CR2l4jMET350to450lowMT2W";
        TFmem.SR = "SR1l4jMET350to450lowMT2W";
    }
    else if(TFRegion == "4jMET450toInflowMT2W")
    {
        TFmem.CR = "CR2l4jMET450toInflowMT2W";
        TFmem.SR = "SR1l4jMET450toInflowMT2W";
    }
    else if(TFRegion == "4jMET250to350highMT2W")
    {
        TFmem.CR = "CR2l4jMET250to350highMT2W";
        TFmem.SR = "SR1l4jMET250to350highMT2W";
    }
    else if(TFRegion == "4jMET350to450highMT2W")
    {
        TFmem.CR = "CR2l4jMET350to450highMT2W";
        TFmem.SR = "SR1l4jMET350to450highMT2W";
    }
    else if(TFRegion == "4jMET450to550highMT2W")
    {
        TFmem.CR = "CR2l4jMET450to550highMT2W";
        TFmem.SR = "SR1l4jMET450to550highMT2W";
    }
    else if(TFRegion == "4jMET550to650highMT2W")
    {
        TFmem.CR = "CR2l4jMET550to650highMT2W";
        TFmem.SR = "SR1l4jMET550to650highMT2W";
    }
    else if(TFRegion == "4jMET650toInfhighMT2W")
    {
        TFmem.CR = "CR2l4jMET650toInfhighMT2W";
        TFmem.SR = "SR1l4jMET650toInfhighMT2W";
    }
    else
        throw std::runtime_error("region not known");

}

void readFactors()
{

    defineFactors();

    dataCR = tab.Get(TFmem.CR,TFmem.data);
    MCSR = tab.Get(TFmem.SR,TFmem.MC);
    MCCR = tab.Get(TFmem.CR,TFmem.MC);

}

Figure returnTF()
{
    readFactors();
    if(MCCR.value() != 0)
        return dataCR*(MCSR/MCCR);
    else 
        return Figure(-1,0);
}
//return TF
};

class TFProducer
{
public:

vector<string> TFRegions;
string TFBkg;

public:

TFProducer(vector<string> regions, string bkgType)
{
    TFRegions.clear();
    TFRegions = regions;
    TFBkg = bkgType;
}

void produceTFTable(string inputName, string outputName)
{
    vector<string> colId;
    colId.push_back("TF");

    vector<string> rowId;
    vector<Figure> value;

    for(uint32_t i=0; i<TFRegions.size(); i++)
    {
        TFFactory fact(inputName, TFRegions.at(i));
        fact.TFmem.MC = TFBkg;
        Figure res = fact.returnTF();

        value.push_back(res);
        rowId.push_back(TFRegions.at(i));
    }

    Table outTab(colId, rowId);

    for(uint32_t j=0; j<rowId.size(); j++)
    {
        outTab.Set(colId.at(0), rowId.at(j), value.at(j));
    }

    outTab.Print(outputName + ".txt");
    outTab.PrintLatex(outputName + ".tex");
}

};

#endif //TFFactory_h
