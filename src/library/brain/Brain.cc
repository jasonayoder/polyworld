/********************************************************************/
/* PolyWorld:  An Artificial Life Ecological Simulator              */
/* by Larry Yaeger                                                  */
/* Copyright Apple Computer 1990,1991,1992                          */
/********************************************************************/

// Self
#include "brain/Brain.h"

#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

// qt
//#include <qapplication.h>

// Local
#include "utils/AbstractFile.h"
#include "agent/agent.h"
#include "sim/debug.h"
#include "sim/globals.h"
#include "brain/groups/GroupsBrain.h"
#include "utils/misc.h"
#include "brain/NervousSystem.h"
#include "brain/NeuronModel.h"
#include "brain/sheets/SheetsBrain.h"
#include "sim/Simulation.h"

using namespace genome;


// Internal globals
#if DebugBrainGrow
	bool DebugBrainGrowPrint = true;
#endif

//===========================================================================
// Brain
//===========================================================================

Brain::Configuration Brain::config;

//---------------------------------------------------------------------------
// Brain::processWorldfile
//---------------------------------------------------------------------------
void Brain::processWorldfile( proplib::Document &doc )
{
	{
		string val = doc.get( "BrainArchitecture" );
		if( val == "Groups" )
			Brain::config.architecture = Brain::Configuration::Groups;
		else if( val == "Sheets" )
			Brain::config.architecture = Brain::Configuration::Sheets;
		else
			assert( false );
	}
	{
		string val = doc.get( "NeuronModel" );
		if( val == "F" )
			Brain::config.neuronModel = Brain::Configuration::FIRING_RATE;
		else if( val == "T" )
			Brain::config.neuronModel = Brain::Configuration::TAU;
		else if( val == "S" )
			Brain::config.neuronModel = Brain::Configuration::SPIKING;
		else
			assert( false );
	}

    Brain::config.Spiking.enableGenes = doc.get( "EnableSpikingGenes" );
	Brain::config.Spiking.aMinVal = doc.get( "SpikingAMin" );
	Brain::config.Spiking.aMaxVal = doc.get( "SpikingAMax" );
	Brain::config.Spiking.bMinVal = doc.get( "SpikingBMin" );
	Brain::config.Spiking.bMaxVal = doc.get( "SpikingBMax" );
	Brain::config.Spiking.cMinVal = doc.get( "SpikingCMin" );
	Brain::config.Spiking.cMaxVal = doc.get( "SpikingCMax" );
	Brain::config.Spiking.dMinVal = doc.get( "SpikingDMin" );
	Brain::config.Spiking.dMaxVal = doc.get( "SpikingDMax" );

	Brain::config.Tau.minVal = doc.get( "TauMin" );
	Brain::config.Tau.maxVal = doc.get( "TauMax" );
	Brain::config.Tau.seedVal = doc.get( "TauSeed" );

    Brain::config.maxbias = doc.get( "MaxBiasWeight" );
	Brain::config.maxneuron2energy = doc.get( "EnergyUseNeurons" );
	Brain::config.outputSynapseLearning = doc.get( "OutputSynapseLearning" );
	Brain::config.synapseFromOutputNeurons = doc.get( "SynapseFromOutputNeurons" );

	Brain::config.numPrebirthCycles = doc.get( "PreBirthCycles" );

    Brain::config.logisticSlope = doc.get( "LogisticSlope" );
    Brain::config.maxWeight = doc.get( "MaxSynapseWeight" );
    Brain::config.initMaxWeight = doc.get( "MaxSynapseWeightInitial" );

    Brain::config.minlrate = doc.get( "MinLearningRate" );
    Brain::config.maxlrate = doc.get( "MaxLearningRate" );

    Brain::config.maxsynapse2energy = doc.get( "EnergyUseSynapses" );
    Brain::config.decayRate = doc.get( "SynapseWeightDecayRate" );
    
    //Genetic Neuron Constraints
    Brain::config.geneticNeuronType = doc.get( "GeneticNeuronType" );
    Brain::config.geneticNeuronPosition = doc.get( "GeneticNeuronPosition" );
    
    Brain::config.allRandomAgents = doc.get( "AllRandomAgents" );
    
    
    // <Gasnets variables
    Brain::config.gasnetsEnabled = doc.get( "GasnetsEnabled" );
    if (Brain::config.gasnetsEnabled) {
        Brain::config.gasnetsDecayRate = doc.get( "GasnetsDecayRate" );
        {
            string val = doc.get( "GasnetsModulationTypeOrder" );
            if( val == "PlaAct" )
                Brain::config.gasnetsModulationTypeOrder = Brain::Configuration::PlaAct;
            else if( val == "ActPla" )
                Brain::config.gasnetsModulationTypeOrder = Brain::Configuration::ActPla;
            else
                assert( false );
        }
        Brain::config.gasnetsGlobalFlatMode = doc.get( "GasnetsGlobalFlatMode" );
        Brain::config.gasnetsConcentrationGradient = doc.get( "GasnetsConcentrationGradient" );
        Brain::config.gasnetsDiscreteEmissionBuildup = doc.get( "GasnetsDiscreteEmissionBuildup" );
        Brain::config.gasnetsFiringRateThresholdBiasBased = doc.get( "GasnetsFiringRateThresholdBiasBased" );
        Brain::config.gasnetsFiringRateThreshold = doc.get( "GasnetsFiringRateThreshold" );
        Brain::config.gasnetsRadiusMin = doc.get( "GasnetsRadiusMin" );
        Brain::config.gasnetsRadiusMax = doc.get( "GasnetsRadiusMax" );
        Brain::config.gasnetsReceptors = doc.get( "GasnetsReceptors" );
        Brain::config.gasnetsGasActivationPercentage = doc.get( "GasnetsGasActivationPercentage" );
        
        Brain::config.gasnetsMinEmissionRate = doc.get( "GasnetsMinEmissionRate" );
        Brain::config.gasnetsMaxEmissionRate = doc.get( "GasnetsMaxEmissionRate" );
        
        Brain::config.gasnetsDiscreteReceptorStrength = doc.get("DiscreteReceptorStrength");
        Brain::config.gasnetsNeuronGasGenerationMin = doc.get( "GasnetsNeuronGasGenerationMin" );
        Brain::config.gasnetsNeuronGasGenerationMax = doc.get( "GasnetsNeuronGasGenerationMax" );
        Brain::config.gasnetsNumGases = doc.get( "GasnetsNumGases" );    
        Brain::config.gasnetsDistancePerTimestep = doc.get( "GasnetsDistancePerTimestep" );    
        Brain::config.gasnetsGasChannelSize  = (int)(Brain::config.gasnetsRadiusMax / Brain::config.gasnetsDistancePerTimestep);
        Brain::config.gasnetsDebugMode  = doc.get( "GasnetsDebugMode" );
        
        
        if (Brain::config.gasnetsDebugMode > 0) {
            cout << " GasnetsDebugMode[1]: " ;
            cout << "Brain::config.gasnetsGasChannelSize: " << Brain::config.gasnetsGasChannelSize << "\n";
        }
        
    } else {
        Brain::config.gasnetsDebugMode = 0;
    }
    
    // Gasnets variables>

 	// Set up retina values
	Brain::config.minWin = doc.get( "RetinaWidth" );

	GroupsBrain::processWorldfile( doc );
	SheetsBrain::processWorldfile( doc );
}

//---------------------------------------------------------------------------
// Brain::init
//---------------------------------------------------------------------------
void Brain::init()
{
    Brain::config.retinaWidth = max(Brain::config.minWin, GroupsBrain::config.maxvisneurpergroup);
    
    if (Brain::config.retinaWidth & 1)
        Brain::config.retinaWidth++;  // keep it even for efficiency (so I'm told)
        
    Brain::config.retinaHeight = Brain::config.minWin;
    
    if (Brain::config.retinaHeight & 1)
        Brain::config.retinaHeight++;

	GroupsBrain::init();
	SheetsBrain::init();
}

//---------------------------------------------------------------------------
// Brain::Brain
//---------------------------------------------------------------------------
Brain::Brain( NervousSystem *cns )
: _cns(cns)
, _neuralnet(NULL)
, _renderer(NULL)
, _energyUse(0)
{
}


//---------------------------------------------------------------------------
// Brain::~Brain
//---------------------------------------------------------------------------
Brain::~Brain()
{
	delete _neuralnet;
	delete _renderer;
}

//---------------------------------------------------------------------------
// Brain::dumpAnatomical
//---------------------------------------------------------------------------
void Brain::dumpAnatomical( AbstractFile *file, long index, float fitness )
{
	// print the header, with index, fitness, and number of neurons
	file->printf( "brain %ld fitness=%g numneurons+1=%d maxWeight=%g maxBias=%g",
				  index, fitness, _dims.numNeurons+1, Brain::config.maxWeight, Brain::config.maxbias );

	_cns->dumpAnatomical( file );
	file->printf( "\n" );

	_neuralnet->dumpAnatomical( file );
}

//---------------------------------------------------------------------------
// Brain::startFunctional
//---------------------------------------------------------------------------
void Brain::startFunctional( AbstractFile *file, long index )
{
	file->printf( "version 1\n" );

	// print the header, with index (agent number)
	file->printf( "brainFunction %ld", index );	

	// print neuron count, number of input neurons, number of synapses
	_neuralnet->startFunctional( file );

	// print timestep born
	file->printf( " %ld", TSimulation::fStep );

	// print organs portion
	_cns->startFunctional( file );

	file->printf( "\n" );
}

//---------------------------------------------------------------------------
// Brain::endFunctional
//---------------------------------------------------------------------------
void Brain::endFunctional( AbstractFile *file, float fitness )
{
	file->printf( "end fitness = %g\n", fitness );
}

//---------------------------------------------------------------------------
// Brain::writeFunctional
//---------------------------------------------------------------------------
void Brain::writeFunctional( AbstractFile *file )
{
	_neuralnet->writeFunctional( file );
}

//---------------------------------------------------------------------------
// Brain::prebirth
//---------------------------------------------------------------------------
void Brain::prebirth()
{
    // now send some signals through the system
    // try pure noise for now...
    for( int i = 0; i < Brain::config.numPrebirthCycles; i++ )
    {
		_cns->prebirthSignal();

		update( false );
    }

    debugcheck( "after prebirth cycling" );
}

//---------------------------------------------------------------------------
// Brain::update
//---------------------------------------------------------------------------
void Brain::update( bool bprint )
{	 
	_neuralnet->update( bprint );
}

//---------------------------------------------------------------------------
// Brain::getRenderer
//---------------------------------------------------------------------------
NeuralNetRenderer *Brain::getRenderer()
{
	return _renderer;
}
