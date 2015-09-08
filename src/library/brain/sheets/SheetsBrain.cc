#include "brain/sheets/SheetsBrain.h"

#include <assert.h>

#include "brain/FiringRateModel.h"
#include "utils/misc.h"
#include "brain/NervousSystem.h"
#include "genome/sheets/SheetsGenome.h"
#include "brain/SpikingModel.h"

using namespace genome;
using namespace sheets;

SheetsBrain::Configuration SheetsBrain::config;

//---------------------------------------------------------------------------
// SheetsBrain::processWorldfile
//---------------------------------------------------------------------------
void SheetsBrain::processWorldfile( proplib::Document &doc )
{
	proplib::Property &sheets = doc.get( "Sheets" );

	SheetsBrain::config.minBrainSize.x = sheets.get( "MinBrainSize" ).get( "X" );
	SheetsBrain::config.minBrainSize.y = sheets.get( "MinBrainSize" ).get( "Y" );
	SheetsBrain::config.minBrainSize.z = sheets.get( "MinBrainSize" ).get( "Z" );

	SheetsBrain::config.maxBrainSize.x = sheets.get( "MaxBrainSize" ).get( "X" );
	SheetsBrain::config.maxBrainSize.y = sheets.get( "MaxBrainSize" ).get( "Y" );
	SheetsBrain::config.maxBrainSize.z = sheets.get( "MaxBrainSize" ).get( "Z" );

	SheetsBrain::config.minSynapseProbabilityX = sheets.get( "MinSynapseProbabilityX" );
	SheetsBrain::config.maxSynapseProbabilityX = sheets.get( "MaxSynapseProbabilityX" );

	SheetsBrain::config.minLearningRate = sheets.get( "MinLearningRate" );
	SheetsBrain::config.maxLearningRate = sheets.get( "MaxLearningRate" );

    SheetsBrain::config.minVisionNeuronsPerSheet = sheets.get( "MinVisionNeuronsPerSheet" );
    SheetsBrain::config.maxVisionNeuronsPerSheet = sheets.get( "MaxVisionNeuronsPerSheet" );

    SheetsBrain::config.minInternalSheetsCount = sheets.get( "MinInternalSheetsCount" );
    SheetsBrain::config.maxInternalSheetsCount = sheets.get( "MaxInternalSheetsCount" );

    SheetsBrain::config.minInternalSheetSize = sheets.get( "MinInternalSheetSize" );
    SheetsBrain::config.maxInternalSheetSize = sheets.get( "MaxInternalSheetSize" );

    SheetsBrain::config.minInternalSheetNeuronCount = sheets.get( "MinInternalSheetNeuronCount" );
    SheetsBrain::config.maxInternalSheetNeuronCount = sheets.get( "MaxInternalSheetNeuronCount" );
}

//---------------------------------------------------------------------------
// SheetsBrain::init
//---------------------------------------------------------------------------
void SheetsBrain::init()
{
}

//---------------------------------------------------------------------------
// SheetsBrain::SheetsBrain
//---------------------------------------------------------------------------
SheetsBrain::SheetsBrain( NervousSystem *cns, SheetsGenome *genome, SheetsModel *model )
: Brain( cns )
, _numInternalSheets( 0 )
, _numInternalNeurons( 0 )
{
	memset( _numSynapses, 0, sizeof(_numSynapses) );
	
	memset( _numGasChannels, 0, sizeof(_numGasChannels) ); //gaschannels

	grow( genome, model );
	delete model;
}

//---------------------------------------------------------------------------
// SheetsBrain::~SheetsBrain
//---------------------------------------------------------------------------
SheetsBrain::~SheetsBrain()
{
}

//---------------------------------------------------------------------------
// SheetsBrain::grow
//---------------------------------------------------------------------------
void SheetsBrain::grow( SheetsGenome *genome, SheetsModel *model )
{

	// ---
	// --- Configure Neuron Count
	// ---
	_dims.numNeurons = (int)model->getNeurons().size();

	// ---
	// --- Configure Synapse Count
	// ---
	itfor( NeuronVector, model->getNeurons(), it )
		_dims.numSynapses += (int)(*it)->synapsesOut.size();
		
		
	// ---
	// --- Configure Input/Output Neurons/Nerves
	// ---
	{
		int numInputSheets = (int)model->getSheets(Sheet::Input).size();
		int numOutputSheets = (int)model->getSheets(Sheet::Output).size();
		int numInternalSheets = (int)model->getSheets(Sheet::Internal).size();
		int numSheets = numInputSheets + numOutputSheets + numInternalSheets;
		int sheetNeuronCount[ numSheets ];
		memset( sheetNeuronCount, 0, sizeof(sheetNeuronCount) );

		itfor( NeuronVector, model->getNeurons(), it )
		{
			Neuron *neuron = *it;
			
			
			switch( neuron->sheet->getType() )
			{
			case Sheet::Input:
				_dims.numInputNeurons++;
				break;
			case Sheet::Output:
				_dims.numOutputNeurons++;
				break;
			case Sheet::Internal:
			
                switch (neuron->attrs.type) 
                {
                case Neuron::Attributes::E:
                case Neuron::Attributes::I:
                case Neuron::Attributes::EI:
                  //do nothing not gas producing
                  break;
                case Neuron::Attributes::G1:
                  if (Brain::config.gasnetsDebugMode > 2) { 
                      cout << "   GasnetsDebugMode[3]: " ;
                      cout << " n[" << neuron->id  << "] is a " ;
                      cout << "G1 Neuron having GasChannels added \n";
                  }
                  //MUST ADD TO DIMS OR YOU WILL GET SEGFAULTS
                  _dims.numGasChannels += neuron->sheet->addGasChannels(neuron, model);
                  break;
                case Neuron::Attributes::G2:
                  if (Brain::config.gasnetsDebugMode > 2) { 
                      cout << "   GasnetsDebugMode[3]: " ;
                      cout << " n[" << neuron->id  << "] is a " ;
                      cout << "G2 Neuron having GasChannels added \n";
                  }
                  //MUST ADD TO DIMS OR YOU WILL GET SEGFAULTS
                  
                  _dims.numGasChannels += neuron->sheet->addGasChannels(neuron, model);
                  
                  break;
                case Neuron::Attributes::G3:
                  if (Brain::config.gasnetsDebugMode > 2) { 
                      cout << "   GasnetsDebugMode[3]: " ;
                      cout << " n[" << neuron->id  << "] is a " ;
                      cout << "G3 Neuron having GasChannels added \n";
                  }
                  //MUST ADD TO DIMS OR YOU WILL GET SEGFAULTS
                  _dims.numGasChannels += neuron->sheet->addGasChannels(neuron, model);
                  break;
                case Neuron::Attributes::G4:
                  if (Brain::config.gasnetsDebugMode > 2) { 
                      cout << "   GasnetsDebugMode[3]: " ;
                      cout << " n[" << neuron->id  << "] is a " ;
                      cout << "G4 Neuron having GasChannels added \n";
                  }
                  //MUST ADD TO DIMS OR YOU WILL GET SEGFAULTS
                  _dims.numGasChannels += neuron->sheet->addGasChannels(neuron, model);
                  break;
                case Neuron::Attributes::G5:
                  if (Brain::config.gasnetsDebugMode > 2) { 
                      cout << "   GasnetsDebugMode[3]: " ;
                      cout << " n[" << neuron->id  << "] is a " ;
                      cout << "G5 Neuron having GasChannels added \n";
                  }
                  //MUST ADD TO DIMS OR YOU WILL GET SEGFAULTS
                  _dims.numGasChannels += neuron->sheet->addGasChannels(neuron, model);
                  break;
                case Neuron::Attributes::G6:
                  if (Brain::config.gasnetsDebugMode > 2) { 
                      cout << "   GasnetsDebugMode[3]: " ;
                      cout << " n[" << neuron->id  << "] is a " ;
                      cout << "G6 Neuron having GasChannels added \n";
                  }
                  //MUST ADD TO DIMS OR YOU WILL GET SEGFAULTS
                  _dims.numGasChannels += neuron->sheet->addGasChannels(neuron, model);
                  break;
                default:
                  assert(false); //should not be M - this is a generic set - see comments in SheetsGenomeSchema.cp
                }
                ////////////
			
				// no-op
				break;
			default:
				assert( false );
			}

			sheetNeuronCount[ neuron->sheet->getId() ]++;
		}

		int neuronIndex = 0;

		for( int sheetId = 0; sheetId < (numInputSheets + numOutputSheets); sheetId++ )
		{
			Sheet *sheet = model->getSheet( sheetId );
			Nerve *nerve = _cns->getNerve( sheet->getName() );
		
			int neuronCount = sheetNeuronCount[sheetId];
			nerve->config( neuronCount, neuronIndex );
			neuronIndex += neuronCount;
		}

		for( int sheetId = numInputSheets + numOutputSheets; sheetId < numSheets; sheetId++ )
		{
			if( sheetNeuronCount[sheetId] )
			{
				_numInternalSheets++;
				_numInternalNeurons += sheetNeuronCount[sheetId];

			}
		}
	}

	// ---
	// --- Instantiate Neural Net 
	// ---
	{
		switch( Brain::config.neuronModel )
		{
		case Brain::Configuration::SPIKING:
			{
				SpikingModel *spiking = new SpikingModel( _cns,
														  genome->get("ScaleLatestSpikes") );
				_neuralnet = spiking;
			}
			break;
		case Brain::Configuration::FIRING_RATE:
		case Brain::Configuration::TAU:
			{
				FiringRateModel *firingRate = new FiringRateModel( _cns );
				_neuralnet = firingRate;
			}
			break;
		default:
			assert(false);
		}

		_neuralnet->init( &_dims, 0.0f );
	}

	// ---
	// --- Configure Neural Net 
	// ---
	{
		int synapseIndex = 0;
		
		//gaschannel
		int targetGaschannelIndex = 0;
		int sourceGaschannelIndex = 0;
		
		

		itfor( NeuronVector, model->getNeurons(), it )
		{
			Neuron *neuron = *it;
			
			//This code is vitally important to giving access to the type of neuron from within the FiringRateModel.cp file
			// M exists, but only as a generic and is unreferenced other than for the purpose of receptive fields
			//so we want to assert false if somehow it ever creeped into actual neuron structs
			int neuronType = 0;
			if (neuron->attrs.type == Neuron::Attributes::E ) {
			    neuronType = 0;
			} else if (neuron->attrs.type == Neuron::Attributes::I ) {
			    neuronType = 1;
			} else if (neuron->attrs.type == Neuron::Attributes::EI ) {
			    neuronType = 2;
			} else if (neuron->attrs.type == Neuron::Attributes::G1 ) {
			    neuronType = 3;
			} else if (neuron->attrs.type == Neuron::Attributes::G2 ) {
			    neuronType = 4;
			} else if (neuron->attrs.type == Neuron::Attributes::G3 ) {
			    neuronType = 5;			    
			} else if (neuron->attrs.type == Neuron::Attributes::G4 ) {
			    neuronType = 6;			    
			} else if (neuron->attrs.type == Neuron::Attributes::G5 ) {
			    neuronType = 7;			    
			} else if (neuron->attrs.type == Neuron::Attributes::G6 ) {
			    neuronType = 8;			    			    
            } else {
                assert(false);
            }
            // end gasnets
            
			_neuralnet->set_neuron( neuron->id,
									&(neuron->attrs.neuronModel),
									synapseIndex,
									-1, 			
									neuronType , 			//gasnets,
									-1,						//activatedByGas
									-1, 					//receptorStrength
									-1, 					//emissionRate
									targetGaschannelIndex, 
									-1,						//targetendgaschannels
									sourceGaschannelIndex,
									-1);					//sourceendgaschannels
            
			itfor( SynapseMap, neuron->synapsesIn, it_synapse )
			{
				Synapse *synapse = it_synapse->second;

				_neuralnet->set_synapse( synapseIndex++,
										 synapse->from->id,
										 synapse->to->id,
										 synapse->attrs.weight,
										 synapse->attrs.lrate );

                //segfault hits here, because synapse's neuron (from or to) must not have a sheet and/or type
				if (synapse->from->sheet->getType() < 0 || synapse->to->sheet->getType() < 0 || synapse->from->sheet->getType() > 2 || synapse->to->sheet->getType() > 2) {
                    cout << "You are going to hit a segfault because you are using a bad index! The can happen when you don't allocate enough space for GasChannels";
                    assert(false);
				}

				_numSynapses[ synapse->from->sheet->getType() ][ synapse->to->sheet->getType() ] += 1;
				
			}

			_neuralnet->set_neuron_endsynapses( neuron->id, synapseIndex );
			
			// setup target gaschannels
			itfor( GasChannelMap, neuron->gasChannelsIn, it_gasChannel )
			{
				GasChannel *gasChannel = it_gasChannel->second;  //accesses the second of the pair
				
				int maxLength = Brain::config.gasnetsGasChannelSize - 1;
				float maxDistance = Brain::config.gasnetsRadiusMax;	//maximum gas radius
				
				int length = (int)(1 +  (maxLength /  (maxDistance / gasChannel->from->absPosition.distance( gasChannel->to->absPosition))));
				
				if (Brain::config.gasnetsDebugMode > 2) {
				    cout << "   GasnetsDebugMode[3]: " ;
    				cout << "absolute distance from n[" << gasChannel->from->id << "] to n[" << gasChannel->to->id << "]: " << gasChannel->from->absPosition.distance( gasChannel->to->absPosition) << " ";
    				cout << "length of gas channel in timesteps: " << length << "\n";
                }
				
				_neuralnet->set_target_gaschannel( targetGaschannelIndex++,
										 gasChannel->from->id,
										 gasChannel->to->id,
										 length);

				_numGasChannels[ gasChannel->from->sheet->getType() ][ gasChannel->to->sheet->getType() ] += 1;
				
			}
			// set target gaschannel endindex			
			_neuralnet->set_neuron_end_target_gaschannels( neuron->id, targetGaschannelIndex );

			
			// setup source gaschannels
			itfor( GasChannelMap, neuron->gasChannelsOut, it_gasChannel )
			{
				GasChannel *gasChannel = it_gasChannel->second;  //accesses the second of the pair
				
				int maxLength = Brain::config.gasnetsGasChannelSize - 1;
				float maxDistance = Brain::config.gasnetsRadiusMax;	//maximum gas radius
				
				int length = (int)(1 +  (maxLength /  (maxDistance / gasChannel->from->absPosition.distance( gasChannel->to->absPosition))));
				
				_neuralnet->set_source_gaschannel( sourceGaschannelIndex++,
										 gasChannel->from->id,
										 gasChannel->to->id,
										 length);

			}
			// set gaschannel endchannels
			_neuralnet->set_neuron_end_source_gaschannels( neuron->id, sourceGaschannelIndex );
			
			
			
		}
	}
}
