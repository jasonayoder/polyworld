#pragma once

#include "BaseNeuronModel.h"

// forward decls
class NervousSystem;

// note: activation levels are not maintained in the neuronstruct
// so that after the new activation levels are computed, the old
// and new blocks of memory can simply be repointered rather than
// copied.
struct FiringRateModel__Neuron
{
	float bias; //re-use as bias for gas emitting neurons as well
	float tau;
	long  startsynapses;
	long  endsynapses;
	
	long  startTargetGaschannels; 	//gasnets
	long  endTargetGaschannels; 	//gasnets
	long  startSourceGaschannels; 	//gasnets
	long  endSourceGaschannels; 	//gasnets
	
	float receptorStrength; //gasnets
	float emissionRate;		//gasnets
	float emissionRadius;	//gasnets
	float emissionMultiplier;	//gasnets
	int activatedByGas;		//gasnets
	
	//genome/sheets/SheetsGenomeSchema.cp::function setIE sets up the location of type in a sheet
	//These are defined to correlate with types in brain/sheets/SheetsBrain.cp::grow(...)
	int   type; //0 = E, 1 = I, 2 = EI, 3 = G1
};

struct FiringRateModel__NeuronAttrs
{
	float bias;  //genetically set bias "bi" from page 6: http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=19FB1EC6C36097BCD43DD35BF53564FF?doi=10.1.1.47.7254&rep=rep1&type=pdf
	float tau;
	float receptorStrength; //gasnets
	float emissionRate; 	//gasnets
	float emissionRadius;	//gasnets
	int type; 				//gasnets
	int activatedByGas;		//gasnets
};

struct FiringRateModel__Synapse
{
	float efficacy;   // > 0 for excitatory, < 0 for inhibitory //unless type is modulatory!
	float lrate;
	short fromneuron; //if the type of fromneuron is modulatory, then it doesn't activate anything... jasonayoder
	short toneuron;
};

struct FiringRateModel__GasChannel
{
	short fromneuron; //must be m-neuron
	short toneuron; //can be m-neuron or std-neuron
	short length; //timestep distance 
};

class FiringRateModel : public BaseNeuronModel<FiringRateModel__Neuron, FiringRateModel__NeuronAttrs, FiringRateModel__Synapse, FiringRateModel__GasChannel>
{
	typedef FiringRateModel__Neuron Neuron;
	typedef FiringRateModel__NeuronAttrs NeuronAttrs;
	typedef FiringRateModel__Synapse Synapse;
	typedef FiringRateModel__GasChannel GasChannel;

 public:
	FiringRateModel( NervousSystem *cns );
	virtual ~FiringRateModel();
	
	virtual void init_derived( float initial_activation );

	virtual void set_neuron( int index,
							 void *attributes,
							 int startsynapses,
							 int endsynapses,
							 int type,						//gasnets
							 int activatedByGas,			//gasnets
							 float receptorStrength, 		//gasnets
							 float emissionRate,			//gasnets
							 int startTargetGaschannels,    //gasnets
							 int endTargetGaschannels,		//gasnets
							 int startSourceGaschannels,	//gasnets
							 int endSourceGaschannels);		//gasnets
							 
	virtual void update( bool bprint );
	
	virtual float getGasActivationImpactAtNeuron(int neuronIndex);
	virtual float getGasPlasticityImpactAtNeuron(int neuronIndex);
	
    //GLOBAL GAS MODE
    float getGasConcentration(int index);                        		//gasnets
    void increaseGasConcentration(int index, float gasConcentration);   //gasnets

protected:
    float _gasConcentration [6]; //TODO cannot be set from Brain::config.gasnetsNumGases... might have to declare here and initialize otherwise
	
};


//GLOBAL gas regulating function
inline float FiringRateModel::getGasConcentration(int i) 
{ 
    return _gasConcentration[i]; 
}
//GLOBAL gas regulating function
inline void FiringRateModel::increaseGasConcentration(int i, float gasConcentration) 
{ 
    _gasConcentration[i] += gasConcentration; 
    if ( _gasConcentration[i] > 1 ) {
        _gasConcentration[i] = 1;
    }
}
