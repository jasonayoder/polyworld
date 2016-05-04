#include "genome/sheets/SheetsGenomeSchema.h"

#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "agent/agent.h"
#include "brain/sheets/SheetsBrain.h"
#include "genome/sheets/SheetsGenome.h"

using namespace std;
using namespace genome;
using namespace sheets;

static const float InputOutputSheetHeight = 0.1;
static const float InputOutputSheetDefaultWidth = 0.1;


SheetsGenomeSchema::Configuration SheetsGenomeSchema::config;

void SheetsGenomeSchema::processWorldfile( proplib::Document &doc )
{
	proplib::Property &sheets = doc.get( "Sheets" );

	{
		proplib::Property &cross = sheets.get( "CrossoverProbability" );

		SheetsGenomeSchema::config.crossoverProbability.sheet = cross.get( "Sheet" );
		SheetsGenomeSchema::config.crossoverProbability.receptiveField = cross.get( "ReceptiveField" );
		SheetsGenomeSchema::config.crossoverProbability.neuronAttr = cross.get( "NeuronAttr" );
		SheetsGenomeSchema::config.crossoverProbability.gene = cross.get( "Gene" );
	}

	{
		string val = sheets.get( "NeuronAttrEncoding" );
		if( val == "Sheet" )
			SheetsGenomeSchema::config.neuronAttrEncoding = SheetsGenomeSchema::Configuration::SheetNeuronAttr;
		else if( val == "Neuron" )
			SheetsGenomeSchema::config.neuronAttrEncoding = SheetsGenomeSchema::Configuration::PerNeuronAttr;
		else
			assert( false );
	}

	{
		string val = sheets.get( "SynapseAttrEncoding" );
		if( val == "Field" )
			SheetsGenomeSchema::config.synapseAttrEncoding = SheetsGenomeSchema::Configuration::FieldSynapseAttr;
		else
			assert( false );
	}

	{
		string val = sheets.get( "ReceptiveFieldEncoding" );
		if( val == "Permutations" )
			SheetsGenomeSchema::config.receptiveFieldEncoding = SheetsGenomeSchema::Configuration::Permutations;
		else if( val == "ExplicitVector" )
		{
			SheetsGenomeSchema::config.receptiveFieldEncoding = SheetsGenomeSchema::Configuration::ExplicitVector;
			SheetsGenomeSchema::config.minExplicitVectorSize = sheets.get( "MinExplicitVectorSize" );
			SheetsGenomeSchema::config.maxExplicitVectorSize = sheets.get( "MaxExplicitVectorSize" );
		}
		else
			assert( false );
	}

	SheetsGenomeSchema::config.enableReceptiveFieldCurrentRegion = sheets.get( "EnableReceptiveFieldCurrentRegion" );
	SheetsGenomeSchema::config.enableReceptiveFieldOtherRegion = sheets.get( "EnableReceptiveFieldOtherRegion" );
}


SheetsGenomeSchema::SheetsGenomeSchema()
{
}

SheetsGenomeSchema::~SheetsGenomeSchema()
{
}

Genome *SheetsGenomeSchema::createGenome( GenomeLayout *layout )
{
	return new SheetsGenome( this, layout );
}

SheetsCrossover *SheetsGenomeSchema::getCrossover()
{
	return &_crossover;
}

// Creates either mutable or immutable scalar gene, dependent on MIN == MAX
#define __SCALAR( CONT, NAME, MIN, MAX, ROUND )						\
	if( (MIN) != (MAX) )												\
		CONT->add( new MutableScalarGene(NAME, MIN, MAX, __InterpolatedGene::ROUND) ); \
	else																\
		CONT->add( new ImmutableScalarGene(NAME, MIN) )

#define SCALAR(CONT, NAME, MIN, MAX) __SCALAR(CONT, NAME, MIN, MAX, ROUND_INT_FLOOR)
#define SCALARJ(CONT, NAME, MIN, MAX) __SCALAR(CONT, NAME, MIN, MAX, ROUND_NONE)
#define SCALARV(CONT, NAME, MINMAX ) SCALAR(CONT, NAME, (MINMAX).a, (MINMAX).b)
#define INDEX(CONT, NAME, MIN, MAX) __SCALAR(CONT, NAME, MIN, MAX, ROUND_INT_BIN)
#define INDEXV(CONT, NAME, MINMAX) INDEX(CONT, NAME, (MINMAX).a, (MINMAX).b )
#define CONST(CONT, NAME, VAL) CONT->add( new ImmutableScalarGene(NAME, VAL) )


template<typename TSheetDef>
static void layoutInputOutputSheets( vector<TSheetDef> &defs )
{
	float centerB = 0.5;
	list<TSheetDef *> rowMembers;

	function<void ()> nextRow =
		[&centerB, &rowMembers]()
		{
			rowMembers.clear();
			if( centerB >= 0.5 )
				// Increase offset, opposite side of center.
				centerB = 1.0 - centerB - InputOutputSheetHeight;
			else
				// Same offset, opposite side of center.
				centerB = 1.0 - centerB;
		};

	for( size_t i = 0; i < defs.size(); i++ )
	{
		TSheetDef *def = &(defs[i]);

		if( !rowMembers.empty() )
		{
			if( (rowMembers.size() % 2) == 0 )
			{
				TSheetDef *end = rowMembers.back();
				if( (end->center.a + (end->maxWidth/2)) + def->maxWidth > 1.0 )
					nextRow();
				else
				{
					def->center.a = (end->center.a + (end->maxWidth/2)) + (def->maxWidth/2);
					def->center.b = centerB;
					rowMembers.push_back( def );
				}
			}
			else
			{
				TSheetDef *end = rowMembers.front();
				if( (end->center.a - (end->maxWidth/2)) - def->maxWidth < 0.0 )
					nextRow();
				else
				{
					def->center.a = (end->center.a - (end->maxWidth/2)) - (def->maxWidth/2);
					def->center.b = centerB;
					rowMembers.push_front( def );
				}
			}
		}

		if( rowMembers.empty() )
		{
			rowMembers.push_back( def );
			def->center.a = 0.5;
			def->center.b = centerB;
		}

		//cout << def->name << " @ " << def->center << endl;
	}
}

void SheetsGenomeSchema::define()
{
	// ---
	// --- Common Genes (e.g. MutationRate, Strength)
	// ---
	GenomeSchema::define();

	_crossover.PhysicalLevel.defineRange( getAll().front(), getAll().back() );

	int sheetId = 0;

	// ---
	// --- Model Properties
	// ---
	SCALAR( this, "SizeX", SheetsBrain::config.minBrainSize.x, SheetsBrain::config.maxBrainSize.x );
	SCALAR( this, "SizeY", SheetsBrain::config.minBrainSize.y, SheetsBrain::config.maxBrainSize.y );
	SCALAR( this, "SizeZ", SheetsBrain::config.minBrainSize.z, SheetsBrain::config.maxBrainSize.z );

	SCALAR( this, "InternalSheetsCount", SheetsBrain::config.minInternalSheetsCount, SheetsBrain::config.maxInternalSheetsCount );

	SCALAR( this, "SynapseProbabilityX", SheetsBrain::config.minSynapseProbabilityX, SheetsBrain::config.maxSynapseProbabilityX );

	SCALAR( this, "LearningRate", SheetsBrain::config.minLearningRate, SheetsBrain::config.maxLearningRate );
	

	// ---
	// --- Input Sheets
	// ---
	ContainerGene *inputSheets = new ContainerGene( "InputSheets" );
	{
		struct InputSheetDef
		{
			const char *name;
			float maxWidth;
			IntMinMax neuronCount;
			bool enabled;

			Vector2f center;
		};
		const IntMinMax VisionNeuronCount( SheetsBrain::config.minVisionNeuronsPerSheet,
										   SheetsBrain::config.maxVisionNeuronsPerSheet );

		// If you are adding a new input, please append to this list in order to maintain stable layout.
		vector<InputSheetDef> defs =
			{
				// Name					MaxWidth						NeuronCount        Enabled
				{ "Red",				1.0,							VisionNeuronCount, true },
				{ "Green",				1.0,							VisionNeuronCount, true },
				{ "Blue",				1.0,							VisionNeuronCount, true },
				{ "Random",				InputOutputSheetDefaultWidth,	1,			       true },
				{ "Energy",				InputOutputSheetDefaultWidth,	1,			       true },
				{ "MateWaitFeedback",	InputOutputSheetDefaultWidth,	1,			       agent::config.enableMateWaitFeedback },
				{ "SpeedFeedback",		InputOutputSheetDefaultWidth,	1,			       agent::config.enableSpeedFeedback },
				{ "Carrying",			InputOutputSheetDefaultWidth,	1,			       agent::config.enableCarry },
				{ "BeingCarried",		InputOutputSheetDefaultWidth,	1,			       agent::config.enableCarry }
			};

		layoutInputOutputSheets( defs );

		for( InputSheetDef def : defs )
		{
			if( def.enabled )
			{
				inputSheets->add(
								 defineSheet( /* name */ def.name,
											  /* id */ sheetId++,
											  /* type*/ Sheet::Input,
											  /* orientation */ (int)PlaneZY,
											  /* slot */ 0.0f,
											  /* centerA */ def.center.a,
											  /* centerB */ def.center.b,
											  /* sizeA */ def.maxWidth,
											  /* sizeB */ InputOutputSheetHeight,
											  /* neuronCountA */ def.neuronCount,
											  /* neuronCountB */ 1,
											  /* scaleToNeuronCountA */ def.neuronCount.a != def.neuronCount.b ) );
			}
		}

		add( inputSheets );
	}

	// ---
	// --- Output Sheets
	// ---
	ContainerGene *outputSheets = new ContainerGene( "OutputSheets" );
	{
		struct OutputSheetDef
		{
			const char *name;
			bool enabled;
			float maxWidth;
			Vector2f center;
		};

		// If you are adding a new output, please append to this list in order to maintain stable layout.
		vector<OutputSheetDef> defs =
			{
				// Name				Enabled
				{ "Eat",			agent::config.enableEat },
				{ "Mate",			agent::config.enableMate },
				{ "Fight",			agent::config.enableFight },
				{ "Speed",			agent::config.enableSpeed },
				{ "Yaw",			agent::config.enableYaw },
				{ "YawOppose",		agent::config.yawEncoding == agent::YE_OPPOSE },
				{ "Light",			agent::config.enableLight },
				{ "Focus",			agent::config.enableFocus },
				{ "VisionPitch",	agent::config.enableVisionPitch },
				{ "VisionYaw",		agent::config.enableVisionYaw },
				{ "Give",			agent::config.enableGive },
				{ "Pickup",			agent::config.enableCarry },
				{ "Drop",			agent::config.enableCarry }
			};

		for( OutputSheetDef &def : defs )
			def.maxWidth = InputOutputSheetDefaultWidth;

		layoutInputOutputSheets( defs );

		for( OutputSheetDef def : defs )
		{
			if( def.enabled )
			{
				outputSheets->add(
								  defineSheet( /* name */ def.name,
											   /* id */ sheetId++,
											   /* type */ Sheet::Output,
											   /* orientation */ (int)PlaneZY,
											   /* slot */ 1.0f,
											   /* centerA */ def.center.a,
											   /* centerB */ def.center.b,
											   /* sizeA */ def.maxWidth,
											   /* sizeB */ InputOutputSheetHeight,
											   /* neuronCountA */ 1,
											   /* neuronCountB */ 1 ) );
			}
		}

		add( outputSheets );
	}

	// ---
	// --- Internal Sheets
	// ---
	ContainerGene *internalSheets = new ContainerGene( "InternalSheets" );
	{
		for( int i = 0; i < SheetsBrain::config.maxInternalSheetsCount; i++ )
		{
			ContainerGene *sheet =
				defineSheet( /* name */ itoa(i),
							 /* id */ sheetId++,
							 /* type */ Sheet::Internal,
							 /* orientation */ IntMinMax( 0, 2 ),
							 /* slot */ FloatMinMax( 0.0f, 1.0f ),
							 /* centerA */ FloatMinMax( 0.0, 1.0 ),
							 /* centerB */ FloatMinMax( 0.0, 1.0 ),
							 /* sizeA */ FloatMinMax( SheetsBrain::config.minInternalSheetSize,
													  SheetsBrain::config.maxInternalSheetSize ),
							 /* sizeB */ FloatMinMax( SheetsBrain::config.minInternalSheetSize,
													  SheetsBrain::config.maxInternalSheetSize ),
							 /* neuronCountA */ IntMinMax( SheetsBrain::config.minInternalSheetNeuronCount,
														   SheetsBrain::config.maxInternalSheetNeuronCount ),
							 /* neuronCountB */ IntMinMax( SheetsBrain::config.minInternalSheetNeuronCount,
														   SheetsBrain::config.maxInternalSheetNeuronCount ) );

			internalSheets->add( sheet );
		}

		add( internalSheets );
	}

	// ---
	// --- Receptive Fields
	// ---
	switch( SheetsGenomeSchema::config.receptiveFieldEncoding )
	{
	case SheetsGenomeSchema::Configuration::Permutations:
		{
			defineReceptiveFields( inputSheets, internalSheets );
			defineReceptiveFields( inputSheets, outputSheets );

			defineReceptiveFields( internalSheets, internalSheets );
			defineReceptiveFields( internalSheets, outputSheets );

			defineReceptiveFields( outputSheets, internalSheets );
			defineReceptiveFields( outputSheets, outputSheets );
		}
		break;
	case SheetsGenomeSchema::Configuration::ExplicitVector:
		{
			int ninputSheets = inputSheets->getAll().size();
			int nsheets = ninputSheets + outputSheets->getAll().size() + internalSheets->getAll().size();

			defineReceptiveFieldVector( inputSheets, Sheet::Target, ninputSheets, nsheets );

			defineReceptiveFieldVector( outputSheets, Sheet::Source, ninputSheets, nsheets );
			defineReceptiveFieldVector( outputSheets, Sheet::Target, ninputSheets, nsheets );

			defineReceptiveFieldVector( internalSheets, Sheet::Source, ninputSheets, nsheets );
			defineReceptiveFieldVector( internalSheets, Sheet::Target, ninputSheets, nsheets );
		}
		break;
	default:
		assert( false );
	}
}

void SheetsGenomeSchema::complete( int offset )
{
	GenomeSchema::complete();

	_crossover.GeneLevel.defineRange( NULL, NULL );

	_crossover.SheetLevel.probability = SheetsGenomeSchema::config.crossoverProbability.sheet;
	_crossover.ReceptiveFieldLevel.probability = SheetsGenomeSchema::config.crossoverProbability.receptiveField;
	_crossover.NeuronAttrLevel.probability = SheetsGenomeSchema::config.crossoverProbability.neuronAttr;
	_crossover.GeneLevel.probability = SheetsGenomeSchema::config.crossoverProbability.gene;

	_crossover.complete( this );

#if false
	{
		SheetsCrossover::Segment segment;

		cout << "--> crossover" << endl;
		while( _crossover.nextSegment(segment) )
		{
			cout << "[" << segment.start << "," << segment.end << "]";
			if( segment.level ) cout << " " << segment.level->name;
			cout << endl;
		}
		exit(0);
	}
#endif
}

Neuron::Attributes::Type SheetsGenomeSchema::getNeuronType( SynapseType synapseType,
															Sheet::ReceptiveFieldNeuronRole role )
{
	switch( role )
	{
	case Sheet::From:
		switch( synapseType )
		{
		case EM:
		case EE:
		case EI:
			return Neuron::Attributes::E;
        case IM:
		case IE:
		case II:
			return Neuron::Attributes::I;
		case ME:
		case MI:
		    //whenever this is returned, we should get a false
			return Neuron::Attributes::M; //gasnets
		default:
			assert( false );
			break;
		}
		break;
	case Sheet::To:
		switch( synapseType )
		{
		case ME:
		case EE:
		case IE:
			return Neuron::Attributes::E;
        case MI:
		case EI:
		case II:
			return Neuron::Attributes::I;
		case EM:
		case IM:
			return Neuron::Attributes::M; //gasnets
		default:
			assert( false );
			break;
		}
		break;
	default:
		assert( false );
	}
}

ContainerGene *SheetsGenomeSchema::defineSheet( const char *name,
												int id,
												Sheet::Type sheetType,
												IntMinMax orientation,
												FloatMinMax slot,
												FloatMinMax centerA,
												FloatMinMax centerB,
												FloatMinMax sizeA,
												FloatMinMax sizeB,
												IntMinMax neuronCountA,
												IntMinMax neuronCountB,
												bool scaleToNeuronCountA )
{
	ContainerGene *sheet = new ContainerGene( name );

	CONST( sheet, "Id", id );
	CONST( sheet, "Type", sheetType );
	INDEXV( sheet, "Orientation", orientation );
	SCALARV( sheet, "Slot", slot );
	SCALARV( sheet, "CenterA", centerA );
	SCALARV( sheet, "CenterB", centerB );
	SCALARV( sheet, "SizeA", sizeA );
	SCALARV( sheet, "SizeB", sizeB );
	SCALARV( sheet, "NeuronCountA", neuronCountA );
	SCALARV( sheet, "NeuronCountB", neuronCountB );
	CONST( sheet, "ScaleToNeuronCountA", (int)scaleToNeuronCountA );

	if( sheetType != Sheet::Input )
	{
		ContainerGene *attrs = new ContainerGene( "NeuronAttrs" );

		if( SheetsGenomeSchema::config.neuronAttrEncoding == Configuration::SheetNeuronAttr )
		{
			defineNeuronAttrs( attrs );
		}
		else if( SheetsGenomeSchema::config.neuronAttrEncoding == Configuration::PerNeuronAttr )
		{
			for( int a = 0; a < neuronCountA.b; a++ )
			{
				char aname[8];
				sprintf( aname, "a%d", a );
				ContainerGene *agene = new ContainerGene( aname );

				for( int b = 0; b < neuronCountB.b; b++ )
				{
					char bname[8];
					sprintf( bname, "b%d", b );
					ContainerGene *bgene = new ContainerGene( bname );

					defineNeuronAttrs( bgene );

					agene->add( bgene );
				}

				attrs->add( agene );
			}
		}

		sheet->add( attrs );

		_crossover.NeuronAttrLevel.defineBoundaryPoints( attrs );
	}

	_crossover.SheetLevel.defineBoundaryPoints( sheet );

	return sheet;
}

void SheetsGenomeSchema::defineNeuronAttrs( ContainerGene *attrs )
{
	SCALAR( attrs, "Bias", -Brain::config.maxbias, Brain::config.maxbias );
	
    if (Brain::config.gasnetsEnabled) {
        
        SCALARJ( attrs, "EmissionRate", Brain::config.gasnetsMinEmissionRate, Brain::config.gasnetsMaxEmissionRate );
        
       	if (Brain::config.gasnetsReceptors) {
       	  SCALARJ( attrs, "ReceptorStrength", 0.0, 1.0 );
       	} else {
       	  SCALARJ( attrs, "ReceptorStrength", 1.0, 1.0 );
       	}
       	
       	//if less than GasnetsGasActivationPercentage it will be a activated by a gas
        SCALARJ( attrs, "ActivationBySynapses", 0.0, 1.0 );
      	SCALARJ( attrs, "EmissionRadius", Brain::config.gasnetsRadiusMin, Brain::config.gasnetsRadiusMax );
	}
	
	
	//normally we have types E or I only 0 or 1
	int numNeuronTypes = 2;
	
	//for each gas there is one more type
	if (Brain::config.gasnetsEnabled) {
	    numNeuronTypes += Brain::config.gasnetsNumGases;
	}
	
	if (Brain::config.geneticNeuronType) {
	  INDEX( attrs, "GeneticType", 0, numNeuronTypes - 1  );
	}
	
	if (Brain::config.geneticNeuronPosition) {
	  //TODO - ADD NEURON POSITION HERE 
	}

	switch( Brain::config.neuronModel )
	{
	case Brain::Configuration::FIRING_RATE:
		// no-op
		break;
	case Brain::Configuration::TAU:
		SCALAR( attrs, "Tau", Brain::config.Tau.minVal, Brain::config.Tau.maxVal );
		break;
	case Brain::Configuration::SPIKING:
		if( Brain::config.Spiking.enableGenes )
		{
			SCALAR( attrs, "SpikingParameterA", Brain::config.Spiking.aMinVal, Brain::config.Spiking.aMaxVal );
			SCALAR( attrs, "SpikingParameterB", Brain::config.Spiking.bMinVal, Brain::config.Spiking.bMaxVal );
			SCALAR( attrs, "SpikingParameterC", Brain::config.Spiking.cMinVal, Brain::config.Spiking.cMaxVal );
			SCALAR( attrs, "SpikingParameterD", Brain::config.Spiking.dMinVal, Brain::config.Spiking.dMaxVal );
		}
		break;
	default:
		assert( false );
	}
}

void SheetsGenomeSchema::defineReceptiveFields( ContainerGene *fromSheets,
												ContainerGene *toSheets )
{
	citfor( GeneVector, fromSheets->getAll(), it_from )
	{
		ContainerGene *from = GeneType::to_Container( *it_from );

		citfor( GeneVector, toSheets->getAll(), it_to )
		{
			ContainerGene *to = GeneType::to_Container( *it_to );
			{
				IntMinMax toId = (int)to->getConst("Id");
				defineReceptiveField( from, toId, Sheet::Target, EE );
				defineReceptiveField( from, toId, Sheet::Target, EI );
				defineReceptiveField( from, toId, Sheet::Target, IE );
				defineReceptiveField( from, toId, Sheet::Target, II );
				
				//only synapses GOING TO M-NEURONS
				if (Brain::config.gasnetsEnabled) {
    				defineReceptiveField( from, toId, Sheet::Target, EM );
	    			defineReceptiveField( from, toId, Sheet::Target, IM );
                }
								
			}

			{
				IntMinMax fromId = (int)from->getConst("Id");
				defineReceptiveField( to, fromId, Sheet::Source, EE );
				defineReceptiveField( to, fromId, Sheet::Source, EI );
				defineReceptiveField( to, fromId, Sheet::Source, IE );
				defineReceptiveField( to, fromId, Sheet::Source, II );
				
				//We only want E->M or I->M synapses
				if (Brain::config.gasnetsEnabled) {
    				defineReceptiveField( to, fromId, Sheet::Source, EM );
    				defineReceptiveField( to, fromId, Sheet::Source, IM );
                }
                
			}
		}
	}
}

void SheetsGenomeSchema::defineReceptiveFieldVector( ContainerGene *sheets,
													 Sheet::ReceptiveFieldRole role,
													 int ninputSheets,
													 int nsheets )
{
	IntMinMax otherSheetId;
	const char *vectorSizeName;
	switch( role )
	{
	case Sheet::Source:
		otherSheetId.set( 0, nsheets - 1 );
		vectorSizeName = "SourceReceptiveFieldsCount";
		break;
	case Sheet::Target:
		otherSheetId.set( ninputSheets, nsheets - 1 );
		vectorSizeName = "TargetReceptiveFieldsCount";
		break;
	default:
		assert( false );
	}

	IntMinMax synapseType( 0, __NumSynapseTypes - 1 );

	for( Gene *sheetGene_ : sheets->getAll() )
	{
		ContainerGene *sheetGene = GeneType::to_Container( sheetGene_ );
		SCALAR( sheetGene,
				vectorSizeName,
				SheetsGenomeSchema::config.minExplicitVectorSize,
				SheetsGenomeSchema::config.maxExplicitVectorSize );

		for( int i = 0; i < SheetsGenomeSchema::config.maxExplicitVectorSize; i++ )
		{
			defineReceptiveField( sheetGene,
								  otherSheetId,
								  role,
								  synapseType );
		}
		
	}
}

void SheetsGenomeSchema::defineReceptiveField( ContainerGene *currentSheet,
											   IntMinMax otherSheetId,
											   Sheet::ReceptiveFieldRole role,
											   IntMinMax synapseType )
{
	const char *fieldArrayName;
	switch( role )
	{
	case Sheet::Target:
		fieldArrayName = "TargetReceptiveFields";
		break;
	case Sheet::Source:
		fieldArrayName = "SourceReceptiveFields";
		break;
	}
	ContainerGene *fieldArray = GeneType::to_Container( currentSheet->gene(fieldArrayName) );
	if( !fieldArray )
	{
		fieldArray = new ContainerGene( fieldArrayName );
		currentSheet->add( fieldArray );
	}
	
	ContainerGene *field = new ContainerGene( itoa(fieldArray->getAll().size()) );

	SCALARV( field, "OtherSheetId", otherSheetId );
	CONST( field, "Role", (int)role );
	SCALARV( field, "SynapseType", synapseType );
	if( SheetsGenomeSchema::config.enableReceptiveFieldCurrentRegion )
	{
		SCALAR( field, "CurrentCenterA", 0.0, 1.0 );
		SCALAR( field, "CurrentCenterB", 0.0, 1.0 );
		SCALAR( field, "CurrentSizeA", 0.0, 1.0 );
		SCALAR( field, "CurrentSizeB", 0.0, 1.0 );
	}
	else
	{
		SCALAR( field, "CurrentCenterA", 0.5, 0.5 );
		SCALAR( field, "CurrentCenterB", 0.5, 0.5 );
		SCALAR( field, "CurrentSizeA", 1.0, 1.0 );
		SCALAR( field, "CurrentSizeB", 1.0, 1.0 );
	}
	if( SheetsGenomeSchema::config.enableReceptiveFieldOtherRegion )
	{
		SCALAR( field, "OtherCenterA", 0.0, 1.0 );
		SCALAR( field, "OtherCenterB", 0.0, 1.0 );
		SCALAR( field, "OtherSizeA", 0.0, 1.0 );
		SCALAR( field, "OtherSizeB", 0.0, 1.0 );
	}
	else
	{
		SCALAR( field, "OtherCenterA", 0.5, 0.5 );
		SCALAR( field, "OtherCenterB", 0.5, 0.5 );
		SCALAR( field, "OtherSizeA", 1.0, 1.0 );
		SCALAR( field, "OtherSizeB", 1.0, 1.0 );
	}
	SCALAR( field, "OffsetA", -1.0f, 1.0f );
	SCALAR( field, "OffsetB", -1.0f, 1.0f );
	SCALAR( field, "SizeA", 0.0f, 1.0f );
	SCALAR( field, "SizeB", 0.0f, 1.0f );

	if( SheetsGenomeSchema::config.synapseAttrEncoding == SheetsGenomeSchema::Configuration::FieldSynapseAttr )
		defineSynapseAttrs( field );

	fieldArray->add( field );

	_crossover.ReceptiveFieldLevel.defineBoundaryPoints( field );
}

void SheetsGenomeSchema::defineSynapseAttrs( ContainerGene *container )
{
	SCALAR( container, "Weight", 0.0f, Brain::config.initMaxWeight );
	SCALAR( container, "LearningRate", 0.0f, 1.0f );
}

SheetsModel *SheetsGenomeSchema::createSheetsModel( SheetsGenome *genome )
{
	Vector3f size( genome->get("SizeX"),
				   genome->get("SizeY"),
				   genome->get("SizeZ") );

	float synapseProbabilityX = genome->get( "SynapseProbabilityX" );

	SheetsModel *model = new SheetsModel( size, synapseProbabilityX );
	
	createSheets( model, genome, genome->gene("InputSheets") );
	createSheets( model, genome, genome->gene("OutputSheets") );
	createSheets( model, genome, genome->gene("InternalSheets"), genome->get("InternalSheetsCount") );

	createReceptiveFields( model, genome, genome->gene("InputSheets") );
	createReceptiveFields( model, genome, genome->gene("OutputSheets") );
	createReceptiveFields( model, genome, genome->gene("InternalSheets") );
	
	model->cull();

	return model;
}

void SheetsGenomeSchema::createSheets( SheetsModel *model, SheetsGenome *g, Gene *sheetsGene, int numSheets )
{
	ContainerGene *container = GeneType::to_Container( sheetsGene );
	if( numSheets == -1 )
		numSheets = (int)container->getAll().size();

	for( int i = 0; i < numSheets; i++ )
	{
		ContainerGene *sheetGene = GeneType::to_Container( container->getAll()[i] );

		createSheet( model, g, sheetGene );
	}
}

void SheetsGenomeSchema::createSheet( SheetsModel *model, SheetsGenome *g, ContainerGene *sheetGene )
{
	// ---
	// --- Fetch sheet attributes
	// ---
	string name = sheetGene->name;
	int id = g->get( sheetGene->gene("Id") );
	Sheet::Type sheetType = (Sheet::Type)(int)g->get( sheetGene->gene("Type") );
	int orientation = g->get( sheetGene->gene("Orientation") );
	float slot = g->get( sheetGene->gene("Slot") );
	Vector2f center( g->get( sheetGene->gene("CenterA") ),
					 g->get( sheetGene->gene("CenterB") ) );
	Vector2f size( g->get( sheetGene->gene("SizeA") ),
					 g->get( sheetGene->gene("SizeB") ) );
	Vector2i neuronCount( g->get( sheetGene->gene("NeuronCountA") ),
						  g->get( sheetGene->gene("NeuronCountB") ) );

	bool scaleToNeuronCountA = (int)sheetGene->getConst( "ScaleToNeuronCountA" ) != 0;
	if( scaleToNeuronCountA )
	{
		int maxNeurons = GeneType::to_MutableScalar( sheetGene->gene("NeuronCountA") )->getMax();
		float fraction = float(neuronCount.a) / maxNeurons;
		size.a *= fraction;
	}

	// ---
	// --- Create neuron attribute configuration function
	// ---
	function<void (Neuron *)> neuronCreated;

	//created this function instead of lots of bad alternative ways of coding
	function<Neuron::Attributes::Type (Neuron *, int num)> convertIntToNeuronTypeEnum = [](Neuron *neuron, int num) {
                    switch (num)
                    {
                        //we skip the IE type because we don't make internal neurons that type
                        case 0:
                          return neuron->attrs.type = Neuron::Attributes::E;
                        case 1:
                          return neuron->attrs.type = Neuron::Attributes::I;
                        case 2:
                          return neuron->attrs.type = Neuron::Attributes::G1;
                        case 3:
                          return neuron->attrs.type = Neuron::Attributes::G2;
                        case 4:
                          return neuron->attrs.type = Neuron::Attributes::G3;
                        case 5:
                          return neuron->attrs.type = Neuron::Attributes::G4;
                        case 6:
                          return neuron->attrs.type = Neuron::Attributes::G5;
                        case 7:
                          return neuron->attrs.type = Neuron::Attributes::G6;
                        default:
                          assert(false);
                    }
                };



	function<void (Neuron *)> setIE =
		[sheetType, neuronCount, convertIntToNeuronTypeEnum] ( Neuron *neuron ) //added neuronCount to capture list, also added function to convert numbers
		{
		
			switch( sheetType )
			{
			case Sheet::Input:
			case Sheet::Output:
				neuron->attrs.type = Neuron::Attributes::EI;
				
				if (Brain::config.gasnetsEnabled && Brain::config.neuronModel == Brain::Configuration::TAU) {
				  neuron->attrs.receptorStrength = neuron->attrs.neuronModel.tau.receptorStrength;
				  //neuron->attrs.activatedByGas = neuron->attrs.neuronModel.tau.activatedByGas;
				}
				
				
				break;
			case Sheet::Internal:
			    
			    //Set TAU variables
                if (Brain::config.neuronModel == Brain::Configuration::TAU  ) 
                {
                    //Gasnets
                    if (Brain::config.gasnetsEnabled) {
                        neuron->attrs.emissionRadius = neuron->attrs.neuronModel.tau.emissionRadius;
                        neuron->attrs.receptorStrength = neuron->attrs.neuronModel.tau.receptorStrength;
                        neuron->attrs.emissionRate = neuron->attrs.neuronModel.tau.emissionRate;
                        neuron->attrs.activatedByGas = neuron->attrs.neuronModel.tau.activatedByGas;
                    }

                    //Genetic Type
                    if (Brain::config.geneticNeuronType) {
                      neuron->attrs.type = convertIntToNeuronTypeEnum( neuron, (int)neuron->attrs.neuronModel.tau.type );
                    }
                
                //Set FIRINGRATE variables
                } else if (Brain::config.neuronModel == Brain::Configuration::FIRING_RATE ) {
                
                    //Gasnets
                    if (Brain::config.gasnetsEnabled) {
                        neuron->attrs.emissionRadius = neuron->attrs.neuronModel.firingRate.emissionRadius;
                        neuron->attrs.receptorStrength = neuron->attrs.neuronModel.firingRate.receptorStrength;
                        neuron->attrs.emissionRate = neuron->attrs.neuronModel.firingRate.emissionRate;
                        neuron->attrs.activatedByGas = neuron->attrs.neuronModel.firingRate.activatedByGas;
                    }
                    
                    //Genetic Type
                    if (Brain::config.geneticNeuronType) {
                        neuron->attrs.type = convertIntToNeuronTypeEnum( neuron, (int)neuron->attrs.neuronModel.firingRate.type );
                      
                    }
                    
                //Not prepared for any other neuronModel
                } else {
                    assert(false);
                }
                
                //for modular arthimetic based neuron types
                if (!Brain::config.geneticNeuronType) {
                    
                    
                    if (!Brain::config.gasnetsEnabled) 
                    {
                        //normal non-genetic non-gas version
                        neuron->attrs.type = (neuron->sheetIndex.a % 2) == (neuron->sheetIndex.b % 2)
                            ? Neuron::Attributes::E
                            : Neuron::Attributes::I;
                    }
                    else
                    {
                        //TODO - redo code here without a nested switch statement
                        if ( Brain::config.gasnetsNumGases == 2 )  
                        {
                            switch ( (neuron->sheetIndex.a + neuron->sheetIndex.b*neuronCount.a) % 6) 
                            {
                                case 0:
                                case 4:
                                  neuron->attrs.type = Neuron::Attributes::E;
                                  break;
                                case 2:
                                case 5:
                                  neuron->attrs.type = Neuron::Attributes::I;
                                  break;
                                case 1:
                                  neuron->attrs.type = Neuron::Attributes::G1;
                                  break;
                                case 3:
                                  neuron->attrs.type = Neuron::Attributes::G2;
                                  break;
                                default:
                                  assert(false);
                              }
                              
                              if (Brain::config.gasnetsDebugMode > 2) {
                                cout << "   GasnetsDebugMode[3]: Neuron created (type: "<<neuron->attrs.type <<") @ (" ;
                                cout << neuron->sheetIndex.a << "," << neuron->sheetIndex.b << ")\n";
                              }
                              
                          } else if(Brain::config.gasnetsNumGases == 4) {
                              
                              switch ( (neuron->sheetIndex.a + neuron->sheetIndex.b*neuronCount.a) % 6) 
                              {
                                case 0:
                                  neuron->attrs.type = Neuron::Attributes::G3;
                                  break;
                                case 4:
                                  neuron->attrs.type = Neuron::Attributes::E;
                                  break;
                                case 2:
                                  neuron->attrs.type = Neuron::Attributes::G4;
                                  break;
                                case 5:
                                  neuron->attrs.type = Neuron::Attributes::I;
                                  break;
                                case 1:
                                  neuron->attrs.type = Neuron::Attributes::G1;
                                  break;
                                case 3:
                                  neuron->attrs.type = Neuron::Attributes::G2;
                                  break;
                                default:
                                  assert(false);
                              }
                              
                              if (Brain::config.gasnetsDebugMode > 2) {
                                cout << "   GasnetsDebugMode[3]: Neuron created (type: "<<neuron->attrs.type <<") @ (" ;
                                cout << neuron->sheetIndex.a << "," << neuron->sheetIndex.b << ")\n";
                              }
                          
                          } else {
                              cout <<"MORE THAN TWO gases????";
                              assert(false);
                          }
                    }
			    } else {

			        //When the neuron type is genetically specified!
			        if (Brain::config.gasnetsDebugMode > 2) {
                                      cout << "   GasnetsDebugMode[3]: Neuron created (type: " << neuron->attrs.type << ") @ " ;
                                      cout << " (" << neuron->sheetIndex.a << "," << neuron->sheetIndex.b << ")\n";
                    }
                    
                    if (neuron->attrs.type < 0 || neuron->attrs.type > 2 + Brain::config.gasnetsNumGases) {
                      assert(false);
                    }
                    
                    
			    }
			    
			    
			    
				break;
			default:
				assert( false );
			}
		};

	if( sheetType == Sheet::Input )
	{
		neuronCreated = [setIE]( Neuron *neuron )
			{
				memset( &neuron->attrs, 0, sizeof(neuron->attrs) );
				setIE( neuron );
			};
	}
	else if( SheetsGenomeSchema::config.neuronAttrEncoding == Configuration::SheetNeuronAttr )
	{
		ContainerGene *attrsGene = GeneType::to_Container( sheetGene->gene( "NeuronAttrs" ) );
		Neuron::Attributes attrs = decodeNeuronAttrs( g, attrsGene );
		
		neuronCreated = [setIE, attrs]( Neuron *neuron )
		{
			neuron->attrs = attrs;
			setIE( neuron );
		};
	}
	else if( SheetsGenomeSchema::config.neuronAttrEncoding == Configuration::PerNeuronAttr )
	{
		ContainerGene *attrsGene = GeneType::to_Container( sheetGene->gene( "NeuronAttrs" ) );

		neuronCreated = [this, g, setIE, attrsGene]( Neuron *neuron )
		{
			ContainerGene *agene = GeneType::to_Container( attrsGene->getAll()[neuron->sheetIndex.a] );
			ContainerGene *bgene = GeneType::to_Container( agene->getAll()[neuron->sheetIndex.b] );
			neuron->attrs = decodeNeuronAttrs( g, bgene );
			setIE( neuron );
		};
		
	}
	else
	{
		assert( false );
	}

	// ---
	// --- Create Sheet
	// ---
	model->createSheet( name,
						id,
						(Sheet::Type)sheetType,
						(Orientation)orientation,
						slot,
						center,
						size,
						neuronCount,
						neuronCreated );
}

Neuron::Attributes SheetsGenomeSchema::decodeNeuronAttrs( SheetsGenome *g,
														  ContainerGene *attrsGene )
{
	Neuron::Attributes attrs;
	
	


	
	switch( Brain::config.neuronModel )
	{
	case Brain::Configuration::FIRING_RATE:
		attrs.neuronModel.firingRate.bias = g->get( attrsGene->gene("Bias") );
		if (Brain::config.gasnetsEnabled) {
		    attrs.neuronModel.firingRate.receptorStrength = g->get( attrsGene->gene("ReceptorStrength") ); 	//gasnets
		    
		    if (Brain::config.gasnetsDiscreteReceptorStrength) {
		        if (attrs.neuronModel.firingRate.receptorStrength < 0.33) {
		            attrs.neuronModel.firingRate.receptorStrength = 0;
		        } else if (attrs.neuronModel.firingRate.receptorStrength < 0.66) {
		            attrs.neuronModel.firingRate.receptorStrength = 0.5;
		        } else {
		            attrs.neuronModel.firingRate.receptorStrength = 1.0;
                }
		    }
		    
		    
		    if (   (float)g->get( attrsGene->gene("ActivationBySynapses") )  < 1.0 - Brain::config.gasnetsGasActivationPercentage ) {
		        attrs.neuronModel.firingRate.activatedByGas = 0;
		    } else {
		    
		        float distributionLevel =  Brain::config.gasnetsGasActivationPercentage / Brain::config.gasnetsNumGases ;
		        //starting at gas 1
		        attrs.neuronModel.firingRate.activatedByGas = 1 +
		            floor(   (  (float)g->get(attrsGene->gene("ActivationBySynapses")) - (1.0 - Brain::config.gasnetsGasActivationPercentage)  )   /  distributionLevel);

		        //catch for the case when an EXACT 1 is returned by the gene
		        if (attrs.neuronModel.firingRate.activatedByGas >  Brain::config.gasnetsNumGases) {
		          attrs.neuronModel.firingRate.activatedByGas -= 1;
		        }
		        
		    }
		    
		    
		    
		    attrs.neuronModel.firingRate.emissionRate = g->get( attrsGene->gene("EmissionRate") ); 			//gasnets
		    attrs.neuronModel.firingRate.emissionRadius = g->get( attrsGene->gene("EmissionRadius") ); 		//gasnets
		}
		
		if (Brain::config.geneticNeuronType) {
		    attrs.neuronModel.firingRate.type = g->get( attrsGene->gene("GeneticType") );
		    
		}
		
		break;
	case Brain::Configuration::TAU:

	    if (Brain::config.gasnetsEnabled) {
    	    attrs.neuronModel.tau.receptorStrength = g->get( attrsGene->gene("ReceptorStrength") ); 		//gasnets
	        
	        if (Brain::config.gasnetsDiscreteReceptorStrength) {
		        if (attrs.neuronModel.tau.receptorStrength < 0.33) {
		            attrs.neuronModel.tau.receptorStrength = 0;
		        } else if (attrs.neuronModel.tau.receptorStrength < 0.66) {
		            attrs.neuronModel.tau.receptorStrength = 0.5;
		        } else {
		            attrs.neuronModel.tau.receptorStrength = 1.0;
                }
		    }
		    
		    //determine if a neuron's gene makes it activate by standard neuron or by a gas
		    //This seems reasonable for some values, but others (at extremes) there seem to be too many
		    //outliers, for instance with 1%, there are some agents with many gas activated neurons- should be much more rare
		    
		    //This current configuration has opposite names so that the lower percentage of gasActivation will be reasonable
		    //Why this isn't symmetric??????
		    //TODO - investigate the reason for the asymmetric distribution of probabilities with genes
		    if (   (float)g->get( attrsGene->gene("ActivationBySynapses") )  < 1.0 - Brain::config.gasnetsGasActivationPercentage ) {
		        attrs.neuronModel.tau.activatedByGas = 0;
		    } else {
		    
		        float distributionLevel =  Brain::config.gasnetsGasActivationPercentage / Brain::config.gasnetsNumGases ;
		        //starting at gas 1
		        attrs.neuronModel.tau.activatedByGas = 1 +
		            floor(   (  (float)g->get(attrsGene->gene("ActivationBySynapses")) - (1.0 - Brain::config.gasnetsGasActivationPercentage)  )   /  distributionLevel);

		        //catch for the case when an EXACT 1 is returned by the gene
		        if (attrs.neuronModel.tau.activatedByGas >  Brain::config.gasnetsNumGases) {
		          attrs.neuronModel.tau.activatedByGas -= 1;
		        }
		        
		    }
	        
	        attrs.neuronModel.tau.emissionRate = g->get( attrsGene->gene("EmissionRate") ); 				//gasnets
	        attrs.neuronModel.tau.emissionRadius = g->get( attrsGene->gene("EmissionRadius") ); 			//gasnets
        }
        
        if (Brain::config.geneticNeuronType) {
            attrs.neuronModel.tau.type = g->get( attrsGene->gene("GeneticType") ) ; 			//gasnets
            
        }
        
		attrs.neuronModel.tau.bias = g->get( attrsGene->gene("Bias") );
		attrs.neuronModel.tau.tau = g->get( attrsGene->gene("Tau") );
		break;
	case Brain::Configuration::SPIKING:
		if( Brain::config.Spiking.enableGenes )
		{
		    if (Brain::config.gasnetsEnabled) {
    		    attrs.neuronModel.spiking.receptorStrength = g->get( attrsGene->gene("ReceptorStrength") ); //gasnets
    		    
    		    if (Brain::config.gasnetsDiscreteReceptorStrength) {
		          if (attrs.neuronModel.spiking.receptorStrength < 0.33) {
    		            attrs.neuronModel.spiking.receptorStrength = 0;
	    	        } else if (attrs.neuronModel.spiking.receptorStrength < 0.66) {
		                attrs.neuronModel.spiking.receptorStrength = 0.5;
                    } else {
		                attrs.neuronModel.spiking.receptorStrength = 1.0;
                    }
                }

                if (   (float)g->get( attrsGene->gene("ActivationBySynapses") )  < 1.0 - Brain::config.gasnetsGasActivationPercentage ) {
		          attrs.neuronModel.spiking.activatedByGas = 0;
                } else {
                
                    float distributionLevel =  Brain::config.gasnetsGasActivationPercentage / Brain::config.gasnetsNumGases ;
                    //starting at gas 1
                    attrs.neuronModel.spiking.activatedByGas = 1 +
                        floor(   (  (float)g->get(attrsGene->gene("ActivationBySynapses")) - (1.0 - Brain::config.gasnetsGasActivationPercentage)  )   /  distributionLevel);

                    //catch for the case when an EXACT 1 is returned by the gene
                    if (attrs.neuronModel.spiking.activatedByGas >  Brain::config.gasnetsNumGases) {
                      attrs.neuronModel.spiking.activatedByGas -= 1;
                    }
                }
    		    
	    	    attrs.neuronModel.spiking.emissionRate = g->get( attrsGene->gene("EmissionRate") ); 		//gasnets
	    	    attrs.neuronModel.spiking.emissionRadius = g->get( attrsGene->gene("EmissionRadius") ); 	//gasnets
	    	    
            }
			attrs.neuronModel.spiking.bias = g->get( attrsGene->gene("Bias") );
			attrs.neuronModel.spiking.SpikingParameter_a = g->get( attrsGene->gene("SpikingParameterA") );
			attrs.neuronModel.spiking.SpikingParameter_b = g->get( attrsGene->gene("SpikingParameterB") );
			attrs.neuronModel.spiking.SpikingParameter_c = g->get( attrsGene->gene("SpikingParameterC") );
			attrs.neuronModel.spiking.SpikingParameter_d = g->get( attrsGene->gene("SpikingParameterD") );
		}

		if (Brain::config.geneticNeuronType) {
            attrs.neuronModel.spiking.type = g->get( attrsGene->gene("GeneticType") ) ;
		}
		
		break;
	default:
		assert( false );
	}			

	return attrs;
}

void SheetsGenomeSchema::createReceptiveFields( SheetsModel *model, SheetsGenome *g, Gene *sheetsGene )
{
	ContainerGene *container = GeneType::to_Container( sheetsGene );

	citfor( GeneVector, container->getAll(), it )
	{
		ContainerGene *sheetGene = GeneType::to_Container( *it );

		int id = g->get( sheetGene->gene("Id") );
		Sheet *sheet = model->getSheet( id );
		
		if( sheet == NULL )
			continue;
			

		vector<const char*> vectorNames = { "SourceReceptiveFields", "TargetReceptiveFields" };

		for( const char *vectorName : vectorNames )
		{
			ContainerGene *fieldsGene = GeneType::to_Container(sheetGene->gene(vectorName));
			if( fieldsGene )
			{
				int nfields;
				switch( SheetsGenomeSchema::config.receptiveFieldEncoding )
				{
				case SheetsGenomeSchema::Configuration::Permutations:
					nfields = fieldsGene->getAll().size();
					break;
				case SheetsGenomeSchema::Configuration::ExplicitVector:
					nfields = g->get( sheetGene->gene(string(vectorName) + "Count") );
					break;
				default:
					assert( false );
				}

				for( int i = 0; i < nfields; i++ )
				{
					ContainerGene *fieldGene = GeneType::to_Container( fieldsGene->getAll()[i] );
		
					createReceptiveField( model, sheet, g, fieldGene );
				}
			}
		}
	}
	
}

void SheetsGenomeSchema::createReceptiveField( SheetsModel *model,
											   Sheet *sheet,
											   SheetsGenome *g,
											   ContainerGene *fieldGene )
{
	// ---
	// --- Fetch field attributes
	// ---
	int otherSheetId = g->get( fieldGene->gene("OtherSheetId") );
	Sheet *otherSheet = model->getSheet( otherSheetId );
	if( otherSheet == NULL )
		return;

	Sheet::ReceptiveFieldRole role = (Sheet::ReceptiveFieldRole)(int)g->get( fieldGene->gene("Role") );
	SynapseType synapseType = (SynapseType)(int)g->get( fieldGene->gene("SynapseType") );
	Vector2f currentCenter( g->get( fieldGene->gene("CurrentCenterA") ),
							g->get( fieldGene->gene("CurrentCenterB") ) );
	Vector2f currentSize( g->get( fieldGene->gene("CurrentSizeA") ),
						  g->get( fieldGene->gene("CurrentSizeB") ) );
	Vector2f otherCenter( g->get( fieldGene->gene("OtherCenterA") ),
							g->get( fieldGene->gene("OtherCenterB") ) );
	Vector2f otherSize( g->get( fieldGene->gene("OtherSizeA") ),
						  g->get( fieldGene->gene("OtherSizeB") ) );
	Vector2f offset( g->get( fieldGene->gene("OffsetA") ),
					 g->get( fieldGene->gene("OffsetB") ) );
	Vector2f size( g->get( fieldGene->gene("SizeA") ),
				   g->get( fieldGene->gene("SizeB") ) );

	// ---
	// --- Create neuron predicate function
	// ---
	function<bool (Neuron *, Sheet::ReceptiveFieldNeuronRole)> neuronPredicate =
		[this, synapseType]( Neuron *neuron, Sheet::ReceptiveFieldNeuronRole role )
		{
		    if (Brain::config.gasnetsEnabled) 
		    {
		        
		        if( neuron->attrs.type == Neuron::Attributes::EI && (getNeuronType( synapseType, role ) == Neuron::Attributes::EI
			                                                    || getNeuronType( synapseType, role ) == Neuron::Attributes::I 
			                                                    || getNeuronType( synapseType, role ) == Neuron::Attributes::E)  ) //gasnets
                    return true;
            
                //if the type is not a standard neuron, then we want to say it matches the generic type
                if (neuron->attrs.type != Neuron::Attributes::EI 
                        && neuron->attrs.type != Neuron::Attributes::E 
                        && neuron->attrs.type != Neuron::Attributes::I 
                        && getNeuronType( synapseType, role ) == Neuron::Attributes::M ) 
                {
                    return true;  
                }

            } 
            else 
            {
                
                if( neuron->attrs.type == Neuron::Attributes::EI )
                                              return true;	
            }
            return getNeuronType( synapseType, role ) == neuron->attrs.type;
		};

	// ---
	// --- Create synapse attribute configuration function
	// ---
	function<void (Synapse *)> synapseCreated;
	Synapse::Attributes attrs;
	if( SheetsGenomeSchema::config.synapseAttrEncoding == SheetsGenomeSchema::Configuration::FieldSynapseAttr )
	{
		WARN_ONCE( "IMPLEMENT INDEX-BASED WEIGHT" ); 
		attrs = decodeSynapseAttrs( g, fieldGene );

		if( getNeuronType(synapseType, Sheet::From) == Neuron::Attributes::I )
		{
			attrs.weight *= -1;
			attrs.lrate *= -1;
		}
		
		synapseCreated =
			[attrs] ( Synapse *synapse )
			{
				synapse->attrs = attrs;
			};
	}
	else
	{
		assert( false );
	}

	sheet->addReceptiveField( role,
							  currentCenter,
							  currentSize,
							  otherCenter,
							  otherSize,
							  offset,
							  size,
							  otherSheet,
							  neuronPredicate,
							  synapseCreated );
							  
}

Synapse::Attributes SheetsGenomeSchema::decodeSynapseAttrs( SheetsGenome *g,
															ContainerGene *container )
{
	Synapse::Attributes attrs;

	attrs.weight = g->get( container->gene("Weight") );
	attrs.lrate = g->get( "LearningRate" );
	attrs.lrate *= (float)g->get( container->gene("LearningRate") );

	return attrs;
}

