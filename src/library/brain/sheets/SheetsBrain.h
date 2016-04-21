#pragma once

#include "brain/Brain.h"
#include "brain/sheets/SheetsModel.h"

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
	int getNumInternalNeuronsOfType(int type); //for stats
	int getNumGasChannelsOfType( int gasIndex ); //for stats
	
	int getNumSynapses( sheets::Sheet::Type from, sheets::Sheet::Type to );
	int getNumGasChannels( sheets::Sheet::Type from, sheets::Sheet::Type to );

 private:
    void grow( genome::SheetsGenome *genome, sheets::SheetsModel *model );

	int _numInternalSheets;
	int _numInternalNeurons;
	
	int _numInternalNeuronsOfType [10];// for stats
	int _numGasChannelsOfType [10]; // for stats
	
	int _numSynapses[sheets::Sheet::__NTYPES][sheets::Sheet::__NTYPES];
	int _numGasChannels[sheets::Sheet::__NTYPES][sheets::Sheet::__NTYPES];
	
};

inline int SheetsBrain::getNumInternalSheets() { return _numInternalSheets; }
inline int SheetsBrain::getNumInternalNeurons() { return _numInternalNeurons; }
inline int SheetsBrain::getNumInternalNeuronsOfType(int type) {   return _numInternalNeuronsOfType[type]; }

inline int SheetsBrain::getNumSynapses( sheets::Sheet::Type from, sheets::Sheet::Type to ) { return _numSynapses[from][to]; }
inline int SheetsBrain::getNumGasChannels( sheets::Sheet::Type from, sheets::Sheet::Type to ) { return _numGasChannels[from][to]; }
inline int SheetsBrain::getNumGasChannelsOfType( int gasIndex ) { return _numGasChannelsOfType[gasIndex]; }


