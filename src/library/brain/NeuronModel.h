#pragma once

#include <stdio.h>

#include <iostream>

// forward decls
class AbstractFile;

#define DebugDumpAnatomical false
#if DebugDumpAnatomical
	#define daPrint( x... ) printf( x );
#else
	#define daPrint( x... )
#endif


class NeuronModel
{
 public:
	
	struct Dimensions
	{
		Dimensions() { numGasChannels = numNeurons = numInputNeurons = numOutputNeurons = numSynapses = 0; }

		int numNeurons;
		int numInputNeurons;
		int numOutputNeurons;

		long numGasChannels; //gasnets
		long numSynapses;
		int numInternalNeuronsOfType[10]; //for stats
		int numGasChannelsOfType[10];  //for stats

		inline int getFirstInputNeuron() { return 0; }
		inline int getFirstOutputNeuron() { return numInputNeurons; }
		inline int getFirstInternalNeuron() { return numInputNeurons + numOutputNeurons; }

		inline int getNumNonInputNeurons() { return numNeurons - numInputNeurons; }
	};

	virtual ~NeuronModel() {}

	virtual void init( Dimensions *dims,
					   float initial_activation ) = 0;

	virtual void set_neuron( int index,
							 void *attributes,
							 int startsynapses = -1,
							 int endsynapses = -1,
							 int type = -1,							//gasnets
							 int activatedByGas = -1,				//gasnets
							 float receptorStrength = -1,			//gasnets
							 float emissionRate = -1,				//gasnets
							 int startTargetGaschannels = -1,		//gasnets
							 int endTargetGaschannels = -1 ,		//gasnets
							 int startSourceGaschannels = -1,		//gasnets
							 int endSourceGaschannels = -1 ) = 0; 	//gasnets
							 
	virtual void set_neuron_endsynapses( int index,
										 int endsynapses ) = 0;
										 
    //gaschannel
    virtual void set_neuron_end_target_gaschannels( int index,
										 int endgaschannels ) = 0;

    //gaschannel
    virtual void set_neuron_end_source_gaschannels( int index,
										 int endgaschannels ) = 0;

    //gaschannel										 
    virtual void set_target_gaschannel(int index,
                                int from,
                                int to,
                                int length) = 0;
                                

   //gaschannel										 
    virtual void set_source_gaschannel(int index,
                                int from,
                                int to,
                                int length) = 0;                                
                                
	
	virtual void set_synapse( int index,
							  int from,
							  int to,
							  float efficacy,
							  float lrate ) = 0;

	virtual void update( bool bprint ) = 0;

	virtual void dumpAnatomical( AbstractFile *file ) = 0;

	virtual void startFunctional( AbstractFile *file ) = 0;
	virtual void writeFunctional( AbstractFile *file ) = 0;
};
