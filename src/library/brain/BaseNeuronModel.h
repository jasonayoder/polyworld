#pragma once

#include <assert.h>
#include <gl.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "sim/globals.h"
#include "utils/AbstractFile.h"
#include "utils/misc.h"
#include "Brain.h"
#include "NervousSystem.h"
#include "NeuronModel.h"

template <typename T_neuron, typename T_neuronattrs, typename T_synapse, typename T_gaschannel>
class BaseNeuronModel : public NeuronModel
{
 public:
	BaseNeuronModel( NervousSystem *cns )
	{
		this->cns = cns;

		neuron = NULL;
		neuronactivation = NULL;
		newneuronactivation = NULL;
		synapse = NULL;
		
		neurongasesconcentrations = NULL;	//gas
		target_gaschannels = NULL;			//gas
		source_gaschannels = NULL;			//gas
		
		

#if PrintBrain
		bprinted = false;
#endif
	}

	virtual ~BaseNeuronModel()
	{
		free( neuron );
		free( neuronactivation );
		free( newneuronactivation );
		
		free (neurongasesconcentrations);
		
		free( synapse );
        free( target_gaschannels );
        free( source_gaschannels );
        
	}

	virtual void init_derived( float initial_activation ) = 0;

	virtual void init( Dimensions *dims,
					   float initial_activation )
	{
		this->dims = dims;
		
		if (Brain::config.gasnetsDebugMode > 1) { 
		    cout << "  GasnetsDebugMode[2]: " ;
    		cout << "Brain created with dims->numGasChannels: " << dims->numGasChannels << "\n";
        }
		
#define __ALLOC(NAME, TYPE, N) if(NAME) free(NAME); NAME = (TYPE *)calloc(N, sizeof(TYPE)); assert(NAME);

		__ALLOC( neuron, T_neuron, dims->numNeurons );
		__ALLOC( neuronactivation, float, dims->numNeurons );
		__ALLOC( newneuronactivation, float, dims->numNeurons );
		__ALLOC( synapse, T_synapse, dims->numSynapses );
		
		__ALLOC( target_gaschannels, T_gaschannel, dims->numGasChannels ); //gasnets  ZERO this results in a segfault down the line - only if youi don't add gaschannels
		__ALLOC( source_gaschannels, T_gaschannel, dims->numGasChannels ); //gasnets  ZERO this results in a segfault down the line - only if youi don't add gaschannels
		

#undef __ALLOC


        // Allocate space for 3d array of neuron gas concentrations
        // neurongasconcentrations[neuron][gas][timestep]  
        // Where gas is index between 0 and *gasnetsNumGases* for gas type and timestep is from 0 to *gasnetsGasChannelSize*
        // and reflects the amount of gas ON THE WAY from various gas channels
        //segfault/free pointer area of concern
        //Imitated other code and follow tutorials online
        
        if (Brain::config.gasnetsEnabled) {
        
            int numberGases = Brain::config.gasnetsNumGases; 
            int maxLengthGasChannel = Brain::config.gasnetsGasChannelSize; 
        
            if (neurongasesconcentrations) {
                free(neurongasesconcentrations);
            }
        
            neurongasesconcentrations = (float *** )calloc(dims->numNeurons, sizeof(float **));
            for (int n = 0; n < dims->numNeurons; n++) {
                neurongasesconcentrations[n] = (float **)calloc(numberGases, sizeof(float *));
                for(int i = 0; i < numberGases; i++) { 
                   neurongasesconcentrations[n][i] = (float *)calloc(maxLengthGasChannel, sizeof(float));
                }
            }
        }

		citfor( NervousSystem::NerveList, cns->getNerves(), it )
		{
			Nerve *nerve = *it;

			nerve->config( &(this->neuronactivation), &(this->newneuronactivation) );  //not needed for gaschannels
		}		

		init_derived( initial_activation );
		
	}

	virtual void set_neuron( int index,
							 void *attributes,
							 int startsynapses,
							 int endsynapses,
							 int type,					//gasnets
							 float receptorStrength,  	//gasnets
							 float emissionRate,		//gasnets
							 int activatedByGas,		//gasnets
							 int startTargetGaschannels,//gaschannel
							 int endTargetGaschannels,	//gasnets
							 int startSourceGaschannels,//gasnets
							 int endSourceGaschannels)	//gaschannel
	{
		T_neuron &n = neuron[index];
		T_neuronattrs *attrs = (T_neuronattrs *)attributes;

		assert( !isnan(attrs->bias) );

		n.bias = attrs->bias;
		
		n.type = type;								  	//gasnets
		n.receptorStrength = attrs->receptorStrength;	//gasnets
		n.emissionRate = attrs->emissionRate;		   	//gasnets
		n.emissionRadius = attrs->emissionRadius;		//gasnets
		n.activatedByGas = attrs->activatedByGas;		//gasnets
		
		n.startsynapses = startsynapses;			  
		n.endsynapses = endsynapses;
		
		n.startTargetGaschannels = startTargetGaschannels; 	//gaschannel
		n.endTargetGaschannels = endTargetGaschannels;		//gaschannel
		n.startSourceGaschannels = startSourceGaschannels; 	//gaschannel
		n.endSourceGaschannels = endSourceGaschannels;		//gaschannel
		

#if DebugBrainGrow
		if( DebugBrainGrowPrint )
		{
            cout << "    bias = " << n.bias nlf;
            cout << "    startsynapses = " << n.startsynapses nlf;
		}
#endif
	}

	virtual void set_neuron_endsynapses( int index,
										 int endsynapses )
	{
		T_neuron &n = neuron[index];

		n.endsynapses = endsynapses;
	}
	
	//gaschannel 
	virtual void set_neuron_end_target_gaschannels( int index,
										 int endTargetGaschannels )
	{
		T_neuron &n = neuron[index];

		n.endTargetGaschannels = endTargetGaschannels;
	}


	//gaschannel 
	virtual void set_neuron_end_source_gaschannels( int index,
										 int endSourceGaschannels )
	{
		T_neuron &n = neuron[index];

		n.endSourceGaschannels = endSourceGaschannels;
	}	
	
	
	
	//gaschannel
	virtual void set_target_gaschannel(int index,
	                            int from, 
	                            int to,
	                            int length ) { //TODO FIXME gasnet length parameter necessary? - possibly add distance here instead
	
      T_gaschannel &g = target_gaschannels[index];

      //Will segfault here if the size is not set properly at the top of this file
      g.fromneuron = from;
      g.toneuron = to;
      g.length = length;
      
	}
	
	
	//gaschannel
	virtual void set_source_gaschannel(int index,
	                            int from, 
	                            int to,
	                            int length ) { //TODO FIXME gasnet length parameter necessary? - possibly add distance here instea
	
      T_gaschannel &g = source_gaschannels[index];

      //Will segfault here if the size is not set properly at the top of this file
      g.fromneuron = from;
      g.toneuron = to;
      g.length = length;
      
	}
	
	
	
	
	
	
	
	

	virtual void set_synapse( int index,
							  int from,
							  int to,
							  float efficacy,
							  float lrate )
	{
		T_synapse &s = synapse[index];

		assert( !isnan(efficacy) );
		assert( !isnan(lrate) );

		s.fromneuron = from;
		s.toneuron = to;
		s.efficacy = efficacy;
		s.lrate = lrate;

#if DebugBrainGrow
		if( DebugBrainGrowPrint )
		{
			cout << "        synapse[" << index
				 << "].toneur, fromneur, efficacy = "
				 << s.toneuron cms s.fromneuron cms s.efficacy nlf;
		}
#endif

	}

	virtual void dumpAnatomical( AbstractFile *file )
	{
		size_t	dimCM;
		float*	connectionMatrix;
		short	i,j;
		long	s;
		float	maxWeight = max( Brain::config.maxWeight, Brain::config.maxbias );
		double	inverseMaxWeight = 1. / maxWeight;
		long imin = 10000;
		long imax = -10000;

		dimCM = (dims->numNeurons+1) * (dims->numNeurons+1);	// +1 for bias neuron
		connectionMatrix = (float*) calloc( sizeof( *connectionMatrix ), dimCM );
		if( !connectionMatrix )
		{
			fprintf( stderr, "%s: unable to alloca connectionMatrix\n", __FUNCTION__ );
			return;
		}

		daPrint( "%s: before filling connectionMatrix\n", __FUNCTION__ );

		// compute the connection matrix
		// columns correspond to presynaptic "from-neurons"
		// rows correspond to postsynaptic "to-neurons"
		for( s = 0; s < dims->numSynapses; s++ )
		{
			int cmIndex;
			cmIndex = abs(synapse[s].fromneuron)  +  abs(synapse[s].toneuron) * (dims->numNeurons + 1);	// +1 for bias neuron
			if( cmIndex < 0 )
			{
				fprintf( stderr, "cmIndex = %d, s = %ld, fromneuron = %d, toneuron = %d, numneurons = %d\n", cmIndex, s, synapse[s].fromneuron, synapse[s].toneuron, dims->numNeurons );
			}
			daPrint( "  s=%d, fromneuron=%d, toneuron=%d, cmIndex=%d\n", s, synapse[s].fromneuron, synapse[s].toneuron, cmIndex );
			connectionMatrix[cmIndex] += synapse[s].efficacy;	// the += is so parallel excitatory and inhibitory connections from input and output neurons just sum together
			if( cmIndex < imin )
				imin = cmIndex;
			if( cmIndex > imax )
				imax = cmIndex;
		}
	
		// fill in the biases
		for( i = 0; i < dims->numNeurons; i++ )
		{
			int cmIndex = dims->numNeurons  +  i * (dims->numNeurons + 1);
			connectionMatrix[cmIndex] = neuron[i].bias;
			if( cmIndex < imin )
				imin = cmIndex;
			if( cmIndex > imax )
				imax = cmIndex;
		}

		if( imin < 0 )
			fprintf( stderr, "%s: cmIndex < 0 (%ld)\n", __FUNCTION__, imin );
		if( imax > (dims->numNeurons+1)*(dims->numNeurons+1) )
			fprintf( stderr, "%s: cmIndex > (numneurons+1)^2 (%ld > %d)\n", __FUNCTION__, imax, (dims->numNeurons+1)*(dims->numNeurons+1) );

		daPrint( "%s: imin = %ld, imax = %ld, numneurons = %d (+1 for bias)\n", __FUNCTION__, imin, imax, dims->numNeurons );

		daPrint( "%s: file = %08lx, index = %ld, fitness = %g\n", __FUNCTION__, (char*)file, index, fitness );

		// print the network architecture
		for( i = 0; i <= dims->numNeurons; i++ )	// running over post-synaptic neurons + bias ('=' because of bias)
		{
			for( j = 0; j <= dims->numNeurons; j++ )	// running over pre-synaptic neurons + bias ('=' because of bias)
				file->printf( "%+06.4f ", connectionMatrix[j + i*(dims->numNeurons+1)] * inverseMaxWeight );
			file->printf( ";\n" );
		}
		
		free( connectionMatrix );
	}

	virtual void startFunctional( AbstractFile *file )
	{
		file->printf( " %d %d %d %ld",
					  dims->numNeurons, dims->numInputNeurons, dims->numOutputNeurons, dims->numSynapses );		
	}

	virtual void writeFunctional( AbstractFile *file )
	{
		for( int i = 0; i < dims->numNeurons; i++ )
		{
			file->printf( "%d %g\n", i, neuronactivation[i] );
		}
	}
	
	//protected:
	NervousSystem *cns;
	Dimensions *dims;

	T_neuron *neuron;
	float *neuronactivation;
	float *newneuronactivation;
	
	float ***neurongasesconcentrations;		//gas
	T_synapse *synapse;						//gas
	T_gaschannel *target_gaschannels;		//gas
	T_gaschannel *source_gaschannels;		//gas

#if PrintBrain
	bool bprinted;
#endif

};

