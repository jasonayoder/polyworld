#pragma once

#include "SheetsModel.h"
#include "brain/Brain.h"
#include "genome/NeuronType.h"

namespace genome { class SheetsGenome; }

class SheetsBrain : public Brain
{
 public:
	static struct Configuration
	{
		sheets::Vector3f minBrainSize;
		sheets::Vector3f maxBrainSize;
		float minSynapseProbabilityX;
		float maxSynapseProbabilityX;
		float minLearningRate;
		float maxLearningRate;
		int minVisionNeuronsPerSheet;
		int maxVisionNeuronsPerSheet;
		int minInternalSheetsCount;
		int maxInternalSheetsCount;
		float minInternalSheetSize;
		float maxInternalSheetSize;
		int minInternalSheetNeuronCount;
		int maxInternalSheetNeuronCount;
	} config;

	static void processWorldfile( proplib::Document &doc );
	static void init();

	SheetsBrain( NervousSystem *cns, genome::SheetsGenome *genome, sheets::SheetsModel *model );
	virtual ~SheetsBrain();

	int getNumInternalSheets();
	int getNumInternalNeurons();
	int getNumSynapses( sheets::Sheet::Type from, sheets::Sheet::Type to );
	
	int numNeuronsOfType(genome::NeuronType type);				//ALIFE14
	long numSynapsesOfType(genome::NeuronType type);			//ALIFE14
	int numInternalNeuronsOfType(genome::NeuronType type);		//ALIFE14
	long numInternalSynapsesOfType(genome::NeuronType type);	//ALIFE14
	

 private:
    void grow( genome::SheetsGenome *genome, sheets::SheetsModel *model );

	int _numModulatoryNeuronsInternal = 0;	//ALIFE14
    int _numInhibitoryNeuronsInternal = 0;	//ALIFE14
    int _numExcitatoryNeuronsInternal = 0;	//ALIFE14
    
    long _numModulatorySynapsesInternal = 0;	//ALIFE14
    long _numInhibitorySynapsesInternal = 0;	//ALIFE14
    long _numExcitatorySynapsesInternal = 0;	//ALIFE14

	int _numInternalSheets;
	int _numInternalNeurons;
	int _numSynapses[sheets::Sheet::__NTYPES][sheets::Sheet::__NTYPES];
};

inline int SheetsBrain::getNumInternalSheets() { return _numInternalSheets; }
inline int SheetsBrain::getNumInternalNeurons() { return _numInternalNeurons; }
inline int SheetsBrain::getNumSynapses( sheets::Sheet::Type from, sheets::Sheet::Type to ) { return _numSynapses[from][to]; }
