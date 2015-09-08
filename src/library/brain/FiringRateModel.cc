#include "brain/FiringRateModel.h"

#include "sim/debug.h"
#include "genome/Genome.h"
#include "genome/GenomeSchema.h"
#include "utils/misc.h"

#include "brain/Brain.h" // temporary

using namespace std;

using namespace genome;


#define GaussianOutputNeurons 0
#if GaussianOutputNeurons
	#define GaussianActivationMean 0.0
	#define GaussianActivationStandardDeviation 1.5
	#define GaussianActivationVariance (GaussianActivationStandardDeviation * GaussianActivationStandardDeviation)
#endif


FiringRateModel::FiringRateModel( NervousSystem *cns )
: BaseNeuronModel<Neuron, NeuronAttrs, Synapse, GasChannel>( cns )
{
}

FiringRateModel::~FiringRateModel()
{
}

void FiringRateModel::init_derived( float initial_activation )
{
	for( int i = 0; i < dims->numNeurons; i++ ) {
		neuronactivation[i] = initial_activation;
		if (Brain::config.gasnetsEnabled && !Brain::config.gasnetsGlobalFlatMode) {
    		for (int j=0;  j < Brain::config.gasnetsNumGases; j++) {
	    	    for (int t=0; t < Brain::config.gasnetsGasChannelSize; t++) {
		            neurongasesconcentrations[i][j][t] = 0;
                }
                
    		}
        }
    }
}

void FiringRateModel::set_neuron( int index,
								  void *attributes,
								  int startsynapses,
								  int endsynapses, 
								  int type,						//gasnets
								  int activatedByGas,			//gasnets
								  float receptorStrength,		//gasnets
								  float emissionRate,			//gasnets
								  int startTargetGaschannels,	//gasnets
								  int endTargetGaschannels,		//gasnets
								  int startSourceGaschannels,	//gaschannels
								  int endSourceGaschannels)		//gaschannels
{
	BaseNeuronModel<Neuron, NeuronAttrs, Synapse, GasChannel>::set_neuron( index,
															   attributes,
															   startsynapses,
															   endsynapses,
															   type,				//gasnets
															   activatedByGas,		//gasnets
															   receptorStrength, 	//gasnets
															   emissionRate,		//gasnets
															   startTargetGaschannels,   //gasnets
															   endTargetGaschannels,     //gasnets
															   startSourceGaschannels,	//gaschannels
															   endSourceGaschannels);		//gaschannels
	NeuronAttrs *attrs = (NeuronAttrs *)attributes;
	Neuron &n = neuron[index];

	assert( !isnan(attrs->tau) );
	n.tau = attrs->tau;
}

// This method takes the index of a neuron (the standard order) and then either returns
// 1.0 if this is a non-gasnets simulation
// or return the difference of gas concentrations based upon local vs. global mode
// and the ordering of the gas types
// This then serves as a multiplication coefficient for the delta of the synaptic weights
// The impact of receptors is handled separately after a call to this method
float FiringRateModel::getGasPlasticityImpactAtNeuron(int neuronIndex) {
    float gasPlasticityImpact  = 1.0;

    if (Brain::config.gasnetsEnabled) {
        
        //GLOBAL MODE
        if (Brain::config.gasnetsGlobalFlatMode)
        {
            switch (Brain::config.gasnetsModulationTypeOrder)
            {
                case Brain::Configuration::PlaAct:
                    gasPlasticityImpact = 1 + getGasConcentration(0) - getGasConcentration(1);
                    break;
                case Brain::Configuration::ActPla:
                    if (Brain::config.gasnetsNumGases >= 4) {
                        gasPlasticityImpact = 1 + getGasConcentration(2) - getGasConcentration(3);
                    }
                    break;
                default:
                    assert(false);
            }
        } else { //LOCAL MODE
            switch (Brain::config.gasnetsModulationTypeOrder)
            {
                case Brain::Configuration::PlaAct:
                    //neurongasesconcentration[neuron][gas#][timestep]
                    gasPlasticityImpact = 1 + neurongasesconcentrations[neuronIndex][1][0] - neurongasesconcentrations[neuronIndex][0][0];
                    break;
                case Brain::Configuration::ActPla:
                    if (Brain::config.gasnetsNumGases >= 4) {
                        gasPlasticityImpact = 1 + neurongasesconcentrations[neuronIndex][2][0] - neurongasesconcentrations[neuronIndex][3][0];
                    }
                    break;
                //TODO FIXME - add more gases (> 4) when appropriate
                default:
                    assert(false);
            }
        }
    }
    return gasPlasticityImpact;
}


// This method takes the index of a neuron (the standard order) and then either returns
// 1.0 if this is a non-gasnets simulation
// or return the difference of gas concentrations based upon local vs. global mode
// and the ordering of the gas types
// This then serves as a multiplication coefficient for neural activations sent to this neuron
// The impact of receptors is handled separately after a call to this method
float FiringRateModel::getGasActivationImpactAtNeuron(int neuronIndex) {

    float gasActivationImpact = 1.0;
    
    if (Brain::config.gasnetsEnabled) {

        if (!Brain::config.gasnetsGlobalFlatMode) {
            switch (Brain::config.gasnetsModulationTypeOrder)
            {
                case Brain::Configuration::ActPla:
                    gasActivationImpact = 1 + neurongasesconcentrations[neuronIndex][0][0] - neurongasesconcentrations[neuronIndex][1][0];
                    break;
                case Brain::Configuration::PlaAct:
                    if (Brain::config.gasnetsNumGases >= 4) {
                        //neurongasesconcentration[neuron][gas#][timestep]
                        gasActivationImpact = 1 + neurongasesconcentrations[neuronIndex][2][0] - neurongasesconcentrations[neuronIndex][3][0];
                    }
                    break;
                //TODO FIXME - add more gases (> 4) when the time is right 
                default:
                    assert(false);
            }
        } else {
            //GLOBAL MODE
            switch (Brain::config.gasnetsModulationTypeOrder)
            {
                case Brain::Configuration::ActPla:
                    gasActivationImpact = 1 + getGasConcentration(0) - getGasConcentration(1);
                    break;
                case Brain::Configuration::PlaAct:
                    if (Brain::config.gasnetsNumGases >= 4) {
                        gasActivationImpact = 1 + getGasConcentration(2) - getGasConcentration(3);
                    }
                    break;
                default:
                    assert(false);
            }
        }
    }
    return gasActivationImpact;
}


void FiringRateModel::update( bool bprint )
{
    debugcheck( "(firing-rate brain) on entry" );
    
    short i;
    long k;
    
    if ((neuron == NULL) || (synapse == NULL) || (neuronactivation == NULL) ) 
        return;
    
    //Only when gasnets enabled
    if (Brain::config.gasnetsEnabled) {

        //LOCAL non-global mode
        if (!Brain::config.gasnetsGlobalFlatMode) {
            if ( neurongasesconcentrations == NULL) {
                cout << "neurongasesconcentrations == NULL\n";
                return;
        	}
        } else { 
        //GLOBAL FLAT MODE
            for (int i=0; i< Brain::config.gasnetsNumGases; i++) {
                float rate = Brain::config.gasnetsDecayRate;
                increaseGasConcentration(i,  -rate * getGasConcentration(i) );
            }
        }
    }
    
	IF_BPRINTED
	(
        printf("neuron (toneuron)  fromneuron   synapse   efficacy\n");
        
        for( i = dims->getFirstOutputNeuron(); i < dims->numNeurons; i++ )
        {
            for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ )
            {
				printf("%3d   %3d        %3d        %5ld         %f\n",
					   i, synapse[k].toneuron, synapse[k].fromneuron,
					   k, synapse[k].efficacy); 
            }
        }
        
        if (Brain::config.gasnetsEnabled && Brain::config.gasnetsDebugMode > 4) {
            for( i = dims->getFirstInternalNeuron(); i < dims->numNeurons; i++ )
            {
                for(int g = neuron[i].startSourceGaschannels; g < neuron[i].endSourceGaschannels; g++ )
                {
                    printf("%3d   %3d        %3d        %5ld         %f\n",
                           i, source_gaschannels[g].toneuron, source_gaschannels[g].fromneuron,
                           g, neuron[source_gaschannels[g].fromneuron].emissionRate); 
                }
            }
        }
	)
	
	//used for ultra verbosity to see the various indicies of the neurons and their categories
	if (Brain::config.gasnetsDebugMode > 5) {
        cout << "         GasnetsDebugMode[6]: \n";
        cout << "           dims->getFirstOutputNeuron()  : " << dims->getFirstOutputNeuron() << "\n";
        cout << "           dims->getFirstInternalNeuron(): " <<dims->getFirstInternalNeuron() << "\n"; 
        cout << "           dims->numNeurons              : " << dims->numNeurons << "\n";
    }

    //purely for debugging purposes
    int debugFirstOut = -1;

	for( i = dims->getFirstOutputNeuron(); i < dims->getFirstInternalNeuron(); i++ )
	{
        newneuronactivation[i] = neuron[i].bias;
        
        if (Brain::config.gasnetsEnabled) {
        
            for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ ) {
                if (i != synapse[k].toneuron ) {
                    assert(false);
                }
                //get the relative impact that the gases should have on this neuron
                float gasActivationImpact = getGasActivationImpactAtNeuron( synapse[k].toneuron );

                //setup receptor level - this is based on the gasnets design
                float receptorLevel = 1;
                if (Brain::config.gasnetsReceptors) {
                    receptorLevel = neuron[synapse[k].toneuron].receptorStrength;
                }
                
                //include receptorLevel in calculation, if its receptorLevel=1 it won't change anything
                gasActivationImpact = 1 + (receptorLevel * (gasActivationImpact - 1)); //this will adjust based on the amount of gas and receptors
                
                if (Brain::config.gasnetsDebugMode > 4 && gasActivationImpact != 1 && debugFirstOut == -1) {
                    debugFirstOut = i;
                    //print stuff about this 
                    cout << "       GasnetsDebugMode[5]: OUTPUT NEURON gasActivationImpact= " << gasActivationImpact << "\n" ;
                    cout << "       GasnetsDebugMode[5]:    W/  gas, newneuronactivation["<< i <<"] +=" << ( synapse[k].efficacy * neuronactivation[synapse[k].fromneuron]) << "\n";
                    cout << "       GasnetsDebugMode[5]:    W/O gas, newneuronactivation["<< i <<"] +=" << ( synapse[k].efficacy * gasActivationImpact * neuronactivation[synapse[k].fromneuron]) << "\n";
                }
                newneuronactivation[i] += synapse[k].efficacy * gasActivationImpact * neuronactivation[synapse[k].fromneuron];
            }

		} else { //standard mode, the above will multiply by 1.0, but this seems to be extra safe and less computational intensive
		
		    for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ ) {
    		    newneuronactivation[i] += synapse[k].efficacy * neuronactivation[synapse[k].fromneuron];
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

	//purely for debugging, but extremely useful and important to trace activity
	//neural activation tracing
	int firstGasActivationImpactedNeuron  = -1; 
	int firstGasActivationImpactingNeuron = -1;
	//synaptic plasticity tracing
	int firstGasPlasticityImpactedNeuron  = -1; 
	int firstGasPlasticityImpactingNeuron = -1;
	

    for( i = dims->getFirstInternalNeuron(); i < numneurons; i++ )
    {
     
		float newactivation = neuron[i].bias;
		
		//STANDARD MODE
		if (!Brain::config.gasnetsEnabled ) {
		
		    //add activation from each incoming synapse
		    for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ ) {
		        newactivation += synapse[k].efficacy * neuronactivation[synapse[k].fromneuron];
		    }
		
		//GASNETS MODE WITH GAS-ACTIVATED NEURON
		} else if ( neuron[i].activatedByGas > 0 ) {
		    
		    float maximumActivationByGas = 1.0; //TODO -  global setting for range of activation - probably not needed after all

		    if (Brain::config.gasnetsGlobalFlatMode) {
		        //global flat mode - activatedByGas starts at 0 for standard neuron
		        newactivation += getGasConcentration( neuron[i].activatedByGas-1) * maximumActivationByGas;
		    } else {
                //local mode
		        newactivation += neurongasesconcentrations[i][ neuron[i].activatedByGas-1 ][0] * maximumActivationByGas;
		    }
		    
		} 
		else //GASNETS WITH NON GAS-ACTIVATED NEURON
		{
            //for each incoming synapse make sure to modulate it according to gas concentrations
            for( k = neuron[i].startsynapses; k < neuron[i].endsynapses; k++ )
            {
                if ( neuron[synapse[k].fromneuron].type >= 3) {
                    //Modulatory synapse must have been created, should NEVER happen
                    assert(false);
                }
                
                //setup gasActivationImpact variable to be used to modulate activity, set based on the local/global mode below
                float gasActivationImpact = FiringRateModel::getGasActivationImpactAtNeuron( synapse[k].toneuron );

                //setup receptor level - this is based on the gasnets design
                float receptorLevel = 1;
                if (Brain::config.gasnetsReceptors) {
                    receptorLevel = neuron[synapse[k].toneuron].receptorStrength;
                }
                
                //include receptorLevel in calculation, if its receptorLevel=1 it won't change anything
                gasActivationImpact = 1 + (receptorLevel * (gasActivationImpact - 1)); //this will adjust based on the amount of gas and receptors
                
                // Demonstrates the difference in the final value according to the gases present, shows what the value would be without gas
                if (Brain::config.gasnetsDebugMode > 3 && gasActivationImpact != 1 && firstGasActivationImpactedNeuron == -1) {
                    firstGasActivationImpactedNeuron = i;
                    cout << "     GasnetsDebugMode[4]:  gasActivationImpact= " << gasActivationImpact << "\n" ;
                    cout << "     GasnetsDebugMode[4]:     W/  gas, newneuronactivation["<< i <<"] +=" << ( synapse[k].efficacy * neuronactivation[synapse[k].fromneuron]) << "\n";
                    cout << "     GasnetsDebugMode[4]:     W/O gas, newneuronactivation["<< i <<"] +=" << ( synapse[k].efficacy * gasActivationImpact * neuronactivation[synapse[k].fromneuron]) << "\n";
                }
                
                newactivation += synapse[k].efficacy * /*added here*/ gasActivationImpact * //gasnets
                   neuronactivation[synapse[k].fromneuron];
                }
                   
		}
		
        //newneuronactivation[i] = logistic(newneuronactivation[i], Brain::config.logisticSlope);


        float storePreSquashedNewActivation = newactivation;

        ////gives reasonable values for difference:  0.35 to -0.62 generally towards the extremes
		float activationThreshold = Brain::config.gasnetsFiringRateThreshold;
		

		if( Brain::config.neuronModel == Brain::Configuration::TAU )
		{
			float tau = neuron[i].tau;
			newactivation = (1.0 - tau) * neuronactivation[i]  +  tau * logistic( newactivation, logisticSlope );
		}
		else
		{
			newactivation = logistic( newactivation, logisticSlope );
		}
		
		//only need to setup emissionMultipliers for neurons that emit gases!
		if (neuron[i].type > 2) {
		
            //basically the greatest emission concentration is at zero distance,  we don't want let it build more than 20 time steps
            float worstCaseDistancePenalty = 1.0;
            
            //diffusion abstraction code here
            if (Brain::config.gasnetsConcentrationGradient) {
                
                float e = 2.718281828;
                //Gasnets uses this: e ^ ((-2*d)/r)
                //consider adding distance to the gaschannel, would make life easier
                //here trying to calculate the worst case scenario so that we can cap multiplier when
                //it would cause the concentration to hit the maximum value
                float r = neuron[i].emissionRadius;
                //the maximum distance that the target can be hit from is the maximum distance the neuron can emit it
                float d = r;// (maxDistance / (maxLength / maxDistance));
                worstCaseDistancePenalty = pow(e, (-2.0*d/r) );
            }
            
            float maximumEmissionMultiplier = (1.0 / neuron[i].emissionRate) / worstCaseDistancePenalty ;
            
            //TODO - actually the gasnets model would have allowed this to go up indefinitely, according to the equations,
            //its just the effects that are capped - must decide how to actually handle this...
            maximumEmissionMultiplier = 100; //min(maximumEmissionMultiplier, (float)20.0);
            
            //store the amount of buildup of gas for each neuron
            if (Brain::config.gasnetsDiscreteEmissionBuildup) {
                //determine if m-neurons are emitting or not and adjusting concentration multiplier
                if (newactivation > activationThreshold ) {
                    neuron[i].emissionMultiplier += 1;
                    if (neuron[i].emissionMultiplier > maximumEmissionMultiplier) { 
                        neuron[i].emissionMultiplier = maximumEmissionMultiplier;
                    }
                } else {
                     neuron[i].emissionMultiplier -= 1;
                     if (neuron[i].emissionMultiplier < 0) {
                         neuron[i].emissionMultiplier = 0;
                     }		
                }
            } else { //floating point buildup/decay
            
                //to keep values consistent with the discrete max/min
                float difference = 2.0 * (newactivation - activationThreshold); //negative when under activationThreshold, positive when OVER
                
                
                //TODO - eliminate Brain::config.gasnetsFiringRateThresholdBiasBased mode?
                if (Brain::config.gasnetsFiringRateThresholdBiasBased) {
                    difference = 2.0 * (logistic( storePreSquashedNewActivation - activationThreshold, logisticSlope) - Brain::config.gasnetsFiringRateThreshold);
                }
                
                if (Brain::config.gasnetsDebugMode > 3) {
                    cout << "     GasnetsDebugMode[4]: " ;
                    if (neuron[i].activatedByGas > 0) {
                        cout << "  neuron["<<i <<"].emissionMultiplier increased by: " << difference << " for a gas activated neuron\n";
                    } else {
                        cout << "  neuron["<< i<<"].emissionMultiplier increased by: " << difference << " for a standard activated neuron\n";
                    }
                }
                
                neuron[i].emissionMultiplier += difference;
                
                if (neuron[i].emissionMultiplier > maximumEmissionMultiplier) {
                    neuron[i].emissionMultiplier = maximumEmissionMultiplier;
                }
                
                if (neuron[i].emissionMultiplier < 0) {
                    neuron[i].emissionMultiplier = 0;
                }
            
            }
		}
		
        /*		LOCAL MODE
        *
        *   We are trying to look at the gas channels associated with this neuron (a gas neuron itself)
        *
        *   We need to send a gas signal based on the activation of this current neuron neuron[i]
        *   
        *	We need to use the currently calculated activation to determine gas output newactivation  -> neuron[i]
        *	We have a list of the gas channels in which the current neuron is a target  neuron[i].startTargetGaschannels -> neuron[i].endTargetGaschannels
        *	We ALSO NOW HAVE a list of the gas channels in which the current neuron is the origin    neuron[i].startSourceGaschannels -> neuron[i].endSourceGaschannels
        *   
        */
        //index i refers to the neuron we have activations for
        //index target refers to the neuron which is a target of the gas channel g
        

        if (Brain::config.gasnetsEnabled) 
        {
            if (!Brain::config.gasnetsGlobalFlatMode) {
            
                
                //*
                
                //Here we use the data structure to go directly to the GasChannels where i is a source
                //The startSourceGaschannels was created specifically here for this purpose
                //If memory space becomes a bigger bottleneck, this could be revert to the very inefficient
                //search through all gaschannels and find only the ones that match
                //This happens for EVERY single neuron though
                //for( int target = dims->getFirstInternalNeuron(); target < numneurons; target++ ) 
                //{
                //    for( int g = neuron[target].startTargetGaschannels; g < neuron[target].endTargetGaschannels; g++ ) 
                //    {
                //        //does the neuron we currently have activations for match the origin for the current synapse in question?
                //        if (i == target_gaschannels[g].fromneuron ) 
                
                
                for (int g = neuron[i].startSourceGaschannels; g < neuron[i].endSourceGaschannels; g++){
                
                    //E, I, EI -> 0, 1, 2 
                    int gasIndex = neuron[source_gaschannels[g].fromneuron].type - 3;
                    
                    if (gasIndex < 0 || gasIndex > 1) {
                        //only handle two gases for now
                        assert(false);
                    }
                    
                    float distancePenalty = 1;
                    
                    //diffusion abstraction code here
                    if (Brain::config.gasnetsConcentrationGradient) {
                        
                        float e = 2.718281828;
                        //Gasnets uses this: e ^ ((-2*d)/r)
                        float r = neuron[i].emissionRadius;
                        float d = (source_gaschannels[g].length-1) * Brain::config.gasnetsDistancePerTimestep;

                        distancePenalty = pow(e, (-2.0*d/r) );
                        
                        //something is wrong if a gaschannel tries to send gas outside of its radius
                        if (d > r) {
                            cout << "d> r\n";
                            cout << "d: " << d << "\n";
                            cout << "r: " << r << "\n";
                            assert(false);
                        }
                        
                    }
                    
                    //multiply by how long neuron has been firing!
                    //multiplier can be based on discrete or floating building
                    
                    //we want to update the target neuron of the emitting neuron
                    //the source_gaschannels[i] variable refers to gaschannels that have neuron[i] as the source neuron
                    //so we want the properties from neuron[i] to generate values, but the concentration should be at that gaschannel's TARGET
                    neurongasesconcentrations[source_gaschannels[g].toneuron][gasIndex][source_gaschannels[g].length] += neuron[i].emissionRate * neuron[i].emissionMultiplier * distancePenalty;
                    
                    //super useful debug print out to see what neuron is impacting neural activations
                    if (Brain::config.gasnetsDebugMode > 3 && Brain::config.gasnetsModulationTypeOrder == Brain::Configuration::ActPla && firstGasActivationImpactingNeuron == -1 ) {
                        firstGasActivationImpactingNeuron = i;
                        cout << "     GasnetsDebugMode[4]: FIRST GAS-NEURON FIRING " ;
                        cout << "neuron[" << i << "] :";
                        cout << " [emissionRate: " <<neuron[i].emissionRate << "] ";
                        cout << " [emissionMultiplier: " << neuron[i].emissionMultiplier << "] ";
                        cout << " [distancePenalty: " << distancePenalty << "] ";
                        cout << "\n";
                        
                        cout << "     GasnetsDebugMode[4]: FIRST GAS-NEURON FIRING " ;
                        cout << "neuron["<<i<<"].type: " << neuron[i].type << " " ;
                        cout << "[bias: " << neuron[i].bias << "] ";
                        cout << "[newactivation: " << newactivation << "] ";
                        cout << "[tau: " << neuron[i].tau << "] ";
                        cout << "[tau - newactivation: " << neuron[i].tau-newactivation << "]\n";
                        
                        //show the indices of the target sites of all the gas channel originating from this neuron
                        cout << "     GasnetsDebugMode[4]: gaschannels from: " << i  << " to [";
                        for (int g2 = neuron[i].startSourceGaschannels; g2 < neuron[i].endSourceGaschannels; g2++) {
                            cout << source_gaschannels[g2].toneuron << ", ";
                        }
                        cout << "]\n";
                    
                    }
                    
                    //maximum concentration of 1
                    if (neurongasesconcentrations[source_gaschannels[g].toneuron][gasIndex][source_gaschannels[g].length] > 1) {
                        neurongasesconcentrations[source_gaschannels[g].toneuron][gasIndex][source_gaschannels[g].length] = 1;
                    }
                }
                
                //only set new activation for neuron if it is NOT an m-neuron
                if (neuron[i].type < 3 ) {
                    newneuronactivation[i] = newactivation;
                } else {
                    //check to see if there is a value here  - when M-neuron originating synapses should have zero activation
                    if (newneuronactivation[i] != 0) {
                        assert(false);
                    }
                }
            } else { 

                //GLOBAL MODE - for gas emitting neurons
                if (neuron[i].type >= 3 ) {
                    // E=0, I=1, EI=2, G1=3, G2=4, G3=5, G4=6   Gas1 should refer to 0 (3-3), and so on
                    int gasIndex = neuron[i].type-3;
                    
                    if (Brain::config.gasnetsDebugMode > 4) {
                        cout << "     GasnetsDebugMode[5]: " ;
                        cout << "neuron[i].type: " << neuron[i].type << " " ;
                        cout << "getGasConcentration("<< (neuron[i].type - 3) << "): " << getGasConcentration(gasIndex) << " ";
                        cout << "newactivation: " << newactivation << " ";
                        cout << "neuron[i].tau: " << neuron[i].tau << " ";
                        cout << "neuron[i].tau-newactivation: " << neuron[i].tau-newactivation << "\n";
                    }
                    /*
                    *  	GAS EMISSION BASED ON NEURAL ACTIVATION
                    *
                    *	GLOBAL MODE: 
                    *
                    * 	one concentration of gas is globally distributed for each gas - we raise it or lower it accordingly
                    */
                    increaseGasConcentration(gasIndex, newactivation  );
                    //gas emitting neurons take the input signal and use it to regulate their gas emissions, 
                    //not their neural activations, which should not do anything anyway
                    newneuronactivation[i] = 0;
                    
                } else {
                    //standard neurons
                    newneuronactivation[i] = newactivation;
                }
            }
        } else {  
            //gasnets not enabled
            newneuronactivation[i] = newactivation;
        }
    }
    
    debugcheck( "after updating neurons" );

	IF_BPRINT
	(
        printf("  i neuron[i].bias neuronactivation[i] newneuronactivation[i]\n");
        for (i = 0; i < dims->numNeurons; i++)
            printf( "%3d  %1.4f  %1.4f  %1.4f\n", i, neuron[i].bias, neuronactivation[i], newneuronactivation[i] );
	)

//	printf( "yaw activation = %g\n", newneuronactivation[yawneuron] );

    float learningrate;
	long numsynapses = dims->numSynapses;
    for (k = 0; k < numsynapses; k++)
    {
		FiringRateModel__Synapse &syn = synapse[k];

		learningrate = syn.lrate;
		
		float efficacy;
		
		if (Brain::config.gasnetsEnabled) {
            
            //get the relative impact the gases should have on the synaptic plasticity
            //A value between 0 and 2
            float gasPlasticityImpact = FiringRateModel::getGasPlasticityImpactAtNeuron(syn.toneuron);
            
            //setup receptor level - this is based on the gasnets design
            float receptorLevel = 1;
            if (Brain::config.gasnetsReceptors) {
                receptorLevel = neuron[syn.toneuron].receptorStrength;
            }
            
            //include receptorLevel in calculation, if its receptorLevel=1 it won't change anything
            gasPlasticityImpact = 1 + (receptorLevel * (gasPlasticityImpact - 1)); //this will adjust based on the amount of gas and receptors

            //consider making this A LOT stronger - so it can swing a synaptic weight WAY UP (they are almost all incredibly small deltas)
            //float plasticityRange =  10.0;
            //gasPlasticityImpact = gasPlasticityImpact * plasticityRange;
            
            //Able to see impact of actual gases on synapses plasticity 
            if (Brain::config.gasnetsDebugMode > 3 && gasPlasticityImpact != 1 && firstGasPlasticityImpactedNeuron == -1 && learningrate != 0 ) {
                firstGasPlasticityImpactedNeuron = syn.toneuron;
                cout << "     GasnetsDebugMode[4]:  gasPlasticityImpact= " << gasPlasticityImpact << "\n" ;
                cout << "     GasnetsDebugMode[4]:     W/  gas, syn.efficiacy= " << (syn.efficacy + learningrate * gasPlasticityImpact * (newneuronactivation[syn.toneuron]-0.5f) * (   neuronactivation[syn.fromneuron]-0.5f) ) << "\n";
                cout << "     GasnetsDebugMode[4]:     W/O gas, syn.efficiacy= " << (syn.efficacy + learningrate * (newneuronactivation[syn.toneuron]-0.5f) * (   neuronactivation[syn.fromneuron]-0.5f) ) << "\n";
            }

            //this is set to be 1 and should not impact the results at all
            //TODO - consider making this additive instead of multiplicative on the delta, might be much more powerful and useful....
            efficacy = syn.efficacy + learningrate /*added here*/ * gasPlasticityImpact //jasonayoder
                * (newneuronactivation[syn.toneuron]-0.5f)
                * (   neuronactivation[syn.fromneuron]-0.5f);
                
        } else {
            //standard unchanged original code, the above code should be fine as gasPlasticityImpact = 1.0 in standard mode, but this may be more efficient
            efficacy = syn.efficacy + learningrate 
                * (newneuronactivation[syn.toneuron]-0.5f)
                * (   neuronactivation[syn.fromneuron]-0.5f);
        }

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
            if (learningrate >= 0.0f)  // excitatory
                efficacy = MAX(0.0f, efficacy);
            if (learningrate < 0.0f)  // inhibitory
                efficacy = MIN(-1.e-10f, efficacy);
        }

		syn.efficacy = efficacy;
    }

    debugcheck( "after updating synapses" );

    float* saveneuronactivation = neuronactivation;
    neuronactivation = newneuronactivation;
    newneuronactivation = saveneuronactivation;

    if (Brain::config.gasnetsEnabled && !Brain::config.gasnetsGlobalFlatMode) 
    {
        int MAX = Brain::config.gasnetsGasChannelSize;
    
        int firstGasImpactedNeuronIndex = -1;
        int firstGasImpactedNeuronGasIndex = -1;
        
        for (int n=0; n< dims->numNeurons; n++) {
            for (int g=0; g < Brain::config.gasnetsNumGases; g++) {
                for( int t=0; t < (MAX - 1); t++) {
                    
                    if (Brain::config.gasnetsModulationTypeOrder == Brain::Configuration::ActPla  && Brain::config.gasnetsDebugMode > 3 && firstGasActivationImpactedNeuron == n) 
                    {
                        if (neurongasesconcentrations[n][g][t] != 0) 
                        {
                            firstGasImpactedNeuronIndex = n;
                            firstGasImpactedNeuronGasIndex = g;
                        }
                    }
                    
                    if (Brain::config.gasnetsModulationTypeOrder == Brain::Configuration::PlaAct && Brain::config.gasnetsDebugMode > 3 && firstGasPlasticityImpactedNeuron == n) 
                    {
                        if (neurongasesconcentrations[n][g][t] != 0) 
                        {
                            firstGasImpactedNeuronIndex = n;
                            firstGasImpactedNeuronGasIndex = g;
                        }
                    }
                    neurongasesconcentrations[n][g][t] = neurongasesconcentrations[n][g][t+1];
                }
                neurongasesconcentrations[n][g][MAX-1] = 0; //always start with emptyspot at end of channel
            }
        }
        
        //eventually we need to re-do this for two separate gases in case 4 gas mode is activated
        if (Brain::config.gasnetsDebugMode > 3) 
        {
            if (firstGasImpactedNeuronIndex != -1) 
            {
                cout << "     GasnetsDebugMode[4]: " ;
                cout << "GasChannel  n[" << firstGasImpactedNeuronIndex << "], g[" << firstGasImpactedNeuronGasIndex << "] : ";
                for( int t=0; t < (MAX - 1); t++) 
                {
                    cout << "["<< t << ": ";
                    cout  << neurongasesconcentrations[firstGasImpactedNeuronIndex][firstGasImpactedNeuronGasIndex][t];
                    cout << "]  ";
                }
                cout << "\n";
            }
        }
    }
}
