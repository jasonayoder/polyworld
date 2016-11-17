#include "FiringRateModel.h"

#include "sim/debug.h"
#include "genome/Genome.h"
#include "genome/GenomeSchema.h"
#include "utils/misc.h"

#include "Brain.h" // temporary

using namespace std;

using namespace genome;


#define GaussianOutputNeurons 0
#if GaussianOutputNeurons
	#define GaussianActivationMean 0.0
	#define GaussianActivationStandardDeviation 1.5
	#define GaussianActivationVariance (GaussianActivationStandardDeviation * GaussianActivationStandardDeviation)
#endif


FiringRateModel::FiringRateModel( NervousSystem *cns )
: BaseNeuronModel<Neuron, NeuronAttrs, Synapse>( cns )
{
}

FiringRateModel::~FiringRateModel()
{
}

void FiringRateModel::init_derived( float initial_activation )
{

	for( int i = 0; i < dims->numNeurons; i++ )
		neuronactivation[i] = initial_activation;
}

void FiringRateModel::set_neuron( int index,
								  void *attributes,
								  int startsynapses,
								  int endsynapses )
{
	BaseNeuronModel<Neuron, NeuronAttrs, Synapse>::set_neuron( index,
															   attributes,
															   startsynapses,
															   endsynapses );
	NeuronAttrs *attrs = (NeuronAttrs *)attributes;
	Neuron &n = neuron[index];

	assert( !isnan(attrs->tau) );
	n.tau = attrs->tau;
}

void FiringRateModel::update( bool bprint )
{
    debugcheck( "(firing-rate brain) on entry" );

    short i;
    long k;
    if ((neuron == NULL) || (synapse == NULL) || (neuronactivation == NULL))
        return;

	IF_BPRINTED
	(
        printf("neuron (toneuron)  fromneuron   synapse   efficacy\n");
        
        for( i = dims->getFirstOutputNeuron(); i < dims->numNeurons; i++ )
        {
            for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ )
            {
				printf("%3d   %3d    %3d    %5ld    %f\n",
					   i, synapse[k].toneuron, synapse[k].fromneuron,
					   k, synapse[k].efficacy); 
            }
        }
	)


	for( i = dims->getFirstOutputNeuron(); i < dims->getFirstInternalNeuron(); i++ )
	{
        newneuronactivation[i] = neuron[i].bias;
		num_m_synapses[i] = 0;		//ALIFE14 m-neuron
		modulatorysignal[i] = 0;	//ALIFE14 modulation
		modulation[i] = 0.0;		//ALIFE14 test

        for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ )
        {
        //ALIFE14 start
        
            //normal operation done if neuromodulation is not enabled or its not a MODULATORY synapse
            if (!Brain::config.neuromodulation || synapse[k].type != MODULATORY ) {
                newneuronactivation[i] += synapse[k].efficacy *
                    neuronactivation[synapse[k].fromneuron];
                        
            } else { 
            
                    if (Brain::config.staticModulatorySynapseWeight ){
                        modulation[i] += synapse[k].modulatoryweight *
                            neuronactivation[synapse[k].fromneuron];
                    } else {
                        modulation[i] += synapse[k].efficacy *
                            neuronactivation[synapse[k].fromneuron];
                    }
                    num_m_synapses[i]++;
                    
            }
        //ALIFE14 end
		}
		

		if (num_m_synapses[i] > 0) {													//ALIFE14 modulation pre-squashing modulation
			modulation[i] /= num_m_synapses[i];											//ALIFE14 modulation
			//modulating activity ONLY													//ALIFE14 modulation
			if (Brain::config.neuromodulation && !Brain::config.modulatePlasticity ) {	//ALIFE14 modulation
				newneuronactivation[i] *= modulation[i];								//ALIFE14 modulation
			}
		}
		

	#if GaussianOutputNeurons
        newneuronactivation[i] = gaussian( newneuronactivation[i], GaussianActivationMean, GaussianActivationVariance );
	#else
		if( Brain::config.neuronModel == Brain::Configuration::TAU )
		{
			float tau = neuron[i].tau;
			newneuronactivation[i] = (1.0 - tau) * neuronactivation[i]  +  tau * logistic( newneuronactivation[i], Brain::config.logisticSlope );
		}
		else
		{
			newneuronactivation[i] = logistic( newneuronactivation[i], Brain::config.logisticSlope );
		}
	#endif
	}

	long numneurons = dims->numNeurons;
	float logisticSlope = Brain::config.logisticSlope;

	//ALIFE14 - purely for debug purposes
    float highmod=0;	//ALIFE14
    float midmod=0;		//ALIFE14
    float lowmod=0;		//ALIFE14
    float maxmod = -1;	//ALIFE14
    float minmod = 2;	//ALIFE14

#if DebugNeuromodulation
    if (DebugNeuromodulationPrint) {
        cout << "\nNumber of Neurons: " << (numneurons - dims->getFirstInternalNeuron()  ) << "\n"; //ALIFE14
    }
#endif
	
    for( i = dims->getFirstInternalNeuron(); i < numneurons; i++ )
    {
		float modulation = 0.0;										//ALIFE14 modulation
		int num_m_synapses = 0;										//ALIFE14 modulation
		float newactivation = neuron[i].bias;
       	for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ )
		{
			if (synapse[k].type == MODULATORY) {					//ALIFE14 modulation
				if (Brain::config.staticModulatorySynapseWeight ){	//ALIFE14 modulation
					modulation += synapse[k].modulatoryweight *		//ALIFE14 modulation
						neuronactivation[synapse[k].fromneuron];	//ALIFE14 modulation
				} else {											//ALIFE14 modulation
					modulation += synapse[k].efficacy *				//ALIFE14 modulation
						neuronactivation[synapse[k].fromneuron];	//ALIFE14 modulation
				}													//ALIFE14 modulation
				num_m_synapses++;									//ALIFE14 modulation
			}														//ALIFE14 modulation
			else {													//ALIFE14 modulation
	            newactivation += synapse[k].efficacy *
    	           neuronactivation[synapse[k].fromneuron];
			}														//ALIFE14 modulation
		}
		
		if (num_m_synapses > 0) {	//ALIFE14 pre-squashing modulation
#if DebugNeuromodulation
            if (DebugNeuromodulationPrint) {		
	    	    cout << num_m_synapses << " ";
            }
#endif
			modulation /= num_m_synapses;												//ALIFE14 modulation
			//if neural activity modulation ONLY										//ALIFE14
			if (Brain::config.neuromodulation && !Brain::config.modulatePlasticity ) {  //ALIFE14 TODO
				newactivation *= modulation;    										//ALIFE14 modulation
			}

		} 
		
        //newneuronactivation[i] = logistic(newneuronactivation[i], Brain::config.logisticSlope);

		if( Brain::config.neuronModel == Brain::Configuration::TAU )
		{
			float tau = neuron[i].tau;
			newactivation = (1.0 - tau) * neuronactivation[i]  +  tau * logistic( newactivation, logisticSlope );
		}
		else
		{
			newactivation = logistic( newactivation, logisticSlope );
		}
        newneuronactivation[i] = newactivation;
		modulatorysignal[i] = modulation; //ALIFE14 modulate plasticity
		
#if DebugNeuromodulation
		if (DebugNeuromodulationPrint) {
		    
		    if (modulation > maxmod) {
		        maxmod = modulation;
		    }
		    if (modulation < minmod) {
		        minmod = modulation;
		    }
		    
		    if (modulation > 1.0) {
		        highmod += 1;//modulation[i];
            }
            if (modulation > 0.1 && modulation <= 1.0 ) {
		        midmod += 1; //modulation[i];
		    } 
		    if ( modulation <= 0.1 ){
		        lowmod += 1;
		    }
		}
#endif						

    }


#if DebugNeuromodulation
    if (DebugNeuromodulationPrint && midmod != -5) {
        cout << "\nhighmod: " << highmod << ", ";
        cout << "midmod: " << midmod << ", ";        
        cout << "lowmod: " << lowmod << ", ";
        cout << "maxmod: " << maxmod << ", ";
        cout << "minmod: " << minmod << "\n";
    }
#endif



    debugcheck( "after updating neurons" );

	IF_BPRINT
	(
        printf("  i neuron[i].bias neuronactivation[i] newneuronactivation[i]\n");
        for (i = 0; i < dims->numNeurons; i++)
            printf( "%3d  %1.4f  %1.4f  %1.4f\n", i, neuron[i].bias, neuronactivation[i], newneuronactivation[i] );
	)

//	printf( "yaw activation = %g\n", newneuronactivation[yawneuron] );


    
    float learningrate;
	int type;
	long numsynapses = dims->numSynapses;
    for (k = 0; k < numsynapses; k++)
    {
		FiringRateModel__Synapse &syn = synapse[k];

		learningrate = syn.lrate;
		type = syn.type;  //ALIFE14 m-neuron
		
		//ALIFE14 m-neuron - definitely add something here #TODO FIX ME
		
#if DebugNeuromodulation
		if (DebugNeuromodulationPrint) {
		    if (modulatorysignal[syn.fromneuron] > highmod) {
		        highmod = modulatorysignal[syn.fromneuron];
            } else if (modulatorysignal[syn.fromneuron] > 0 && modulatorysignal[syn.fromneuron] < lowmod) {
		        lowmod = modulatorysignal[syn.fromneuron];
		    } else {
		        midmod = modulatorysignal[syn.fromneuron];
		    }
//		    cout << "syn.efficacy: " << syn.efficacy << ", ";
//		    cout << "modulatorysignal[syn.fromneuron]: " << modulatorysignal[syn.fromneuron] << ", ";
		}
#endif
	
		float efficacy = 0.0;															//ALIFE14
		if (Brain::config.modulatePlasticity) {											//ALIFE14
			if (modulatorysignal[syn.fromneuron] == 0) {								//ALIFE14
				modulatorysignal[syn.fromneuron] = 1;     ///varies between 0 and 2.0	//ALIFE14
			}																			//ALIFE14
			efficacy = syn.efficacy +  modulatorysignal[syn.fromneuron] * learningrate	//ALIFE14 modulatorySignal pre-synaptic neuron IMPORTANT
				* (newneuronactivation[syn.toneuron]-0.5f)								//ALIFE14
				* (   neuronactivation[syn.fromneuron]-0.5f);							//ALIFE14
		} else {																		//ALIFE14
			efficacy = syn.efficacy + learningrate										//ALIFE14
				* (newneuronactivation[syn.toneuron]-0.5f)								//ALIFE14
				* (   neuronactivation[syn.fromneuron]-0.5f);							//ALIFE14
		}																				//ALIFE14

#if DebugNeuromodulation
		if (DebugNeuromodulationPrint) {
//		    cout << "efficacy: " << efficacy << "\n";
		}
#endif		

        if (fabs(efficacy) > (0.5f * Brain::config.maxWeight))
        {
            efficacy *= 1.0f - (1.0f - Brain::config.decayRate) *
                (fabs(efficacy) - 0.5f * Brain::config.maxWeight) / (0.5f * Brain::config.maxWeight);
            if (efficacy > Brain::config.maxWeight)
                efficacy = Brain::config.maxWeight;
            else if (efficacy < -Brain::config.maxWeight)
                efficacy = -Brain::config.maxWeight;
        }
        else
        {
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
            // not strictly correct for this to be in an else clause,
            // but if lrate is reasonable, efficacy should never change
            // sign with a new magnitude greater than 0.5 * Brain::config.maxWeight
		    // ALIFE14 m-neuron, don't think anything is necessary to change here, double check TODO
            if (learningrate >= 0.0f)  // excitatory
                efficacy = MAX(0.0f, efficacy);
            if (learningrate < 0.0f)  // inhibitory
                efficacy = MIN(-1.e-10f, efficacy);
        }

		syn.efficacy = efficacy;
    }

#if DebugNeuromodulation
    if (DebugNeuromodulationPrint) {
//        cout << "highmod: " << highmod << ",";
//        cout << "midmod: " << midmod << ",";        
//        cout << "lowmod: " << lowmod << "\n";
    }
#endif

    debugcheck( "after updating synapses" );

    float* saveneuronactivation = neuronactivation;
    neuronactivation = newneuronactivation;
    newneuronactivation = saveneuronactivation;
}
