#include "brain/sheets/SheetsModel.h"

#include <stdlib.h>

#include "utils/misc.h"

using namespace sheets;
using namespace std;

#define trc(x...)
//#define trc(x...) cout << x << endl;

// verbose trace
#define vtrc(x...)
//#define vtrc(x...) cout << x << endl;

namespace sheets
{
	std::ostream &operator<<( std::ostream &out, const NeuronSubset &n )
	{
		return out << "{" << n._begin << " --> " << n._end << "}";
	}
}


//===========================================================================
// Neuron
//===========================================================================

bool NeuronKeyCompare::operator()( Neuron *x, Neuron *y ) 
{
	return x->nonCulledId < y->nonCulledId;
}

//---------------------------------------------------------------------------
// Neuron::~Neuron
//---------------------------------------------------------------------------
Neuron::~Neuron()
{
	itfor( SynapseMap, synapsesOut, it )
		delete it->second;
		
	itfor( GasChannelMap, gasChannelsOut, it )
	        delete it->second;
	
}

//===========================================================================
// NeuronSubset
//===========================================================================

NeuronSubset::Iterator::Iterator( NeuronSubset *subset, const Vector2i &index ) : _subset(subset), _index(index)
{
}

Vector2i &NeuronSubset::Iterator::operator*()
{
	return _index;
}

NeuronSubset::Iterator &NeuronSubset::Iterator::operator++()
{
	_index.b++;
	if( _index.b > _subset->_end.b )
	{
		_index.a++;
		_index.b = _subset->_begin.b;
	}

	return *this;
}

bool NeuronSubset::Iterator::operator!=( const Iterator &other ) const
{
	return _index.a != other._index.a || _index.b != other._index.b;
}

size_t NeuronSubset::size()
{
	return max(0, _end.a - _begin.a + 1) * max(0, _end.b - _begin.b + 1);
}

//===========================================================================
// Sheet
//===========================================================================

const char *Sheet::getName( Type type )
{
	switch( type )
	{
	case Input: return "Input";
	case Output: return "Output";
	case Internal: return "Internal";
	default: assert( false );
	}
}

//---------------------------------------------------------------------------
// Sheet::Sheet
//---------------------------------------------------------------------------
Sheet::Sheet( SheetsModel *model,
			  string name,
			  int id,
			  Type type,
			  Orientation orientation,
			  float slot,
			  const Vector2f &center,
			  const Vector2f &size,
			  const Vector2i &neuronCount,
			  function<void (Neuron*)> &neuronCreated )
: _sheetsModel( model )
, _name( name )
, _id( id )
, _type( type )
, _orientation( orientation )
, _slot( slot )
, _center( center )
, _size( size )
, _neuronCount( neuronCount )
, _neuronSpacing( 1.0f / _neuronCount.a,
				  1.0f / _neuronCount.b )
, _neuronInsets( _neuronSpacing / 2 )
, _nneurons( neuronCount.a * neuronCount.b )
{
	trc( "Sheet [orig] Center=" << _center << " Size=" << _size << " --> Corners=[ " << (_center - (_size / 2)) << " , " << (_center + (_size / 2)) << " ]" );

	// If the sheet exceeds the model boundaries, then we must scale it down and translate it to fit.
	trim( _center.a, _size.a );
	trim( _center.b, _size.b );

	trc( "Sheet [trim] Center=" << _center << " Size=" << _size << " --> Corners[ " << (_center - (_size / 2)) << " , " << (_center + (_size / 2)) << " ]" );

	createNeurons( neuronCreated );
}

//---------------------------------------------------------------------------
// Sheet::~Sheet
//---------------------------------------------------------------------------
Sheet::~Sheet()
{
	delete [] _neurons;
}

//---------------------------------------------------------------------------
// Sheet::getName
//---------------------------------------------------------------------------
string Sheet::getName()
{
	return _name;
}

//---------------------------------------------------------------------------
// Sheet::getId
//---------------------------------------------------------------------------
int Sheet::getId()
{
	return _id;
}

//---------------------------------------------------------------------------
// Sheet::getType
//---------------------------------------------------------------------------
Sheet::Type Sheet::getType()
{
	return _type;
}

//---------------------------------------------------------------------------
// Sheet::getNeuronCount
//---------------------------------------------------------------------------
const Vector2i &Sheet::getNeuronCount()
{
	return _neuronCount;
}

//---------------------------------------------------------------------------
// Sheet::getNeuron
//---------------------------------------------------------------------------
Neuron *Sheet::getNeuron( int a, int b )
{
	//cout << "a: " << a <<"\n";
	//cout << "_neuronCount.a: " << _neuronCount.a <<"\n";
	

	assert( a < _neuronCount.a );
	assert( a >= 0 );
	assert( b < _neuronCount.b );
	assert( b >= 0 );

	int offset = (a + (b * _neuronCount.a));
	assert( (offset >= 0) && (offset < _nneurons) );

	return _neurons + offset;
}

//---------------------------------------------------------------------------
// Sheet::getNeuron
//---------------------------------------------------------------------------
Neuron *Sheet::getNeuron( Vector2i index )
{
	return getNeuron( index.a, index.b );
}

//---------------------------------------------------------------------------
// Sheet::addReceptiveField
//---------------------------------------------------------------------------
void Sheet::addReceptiveField( ReceptiveFieldRole role,
							   Vector2f currentCenter,
							   Vector2f currentSize,
							   Vector2f otherCenter,
							   Vector2f otherSize,
							   Vector2f fieldOffset,
							   Vector2f fieldSize,
							   Sheet *other,
							   function<bool (Neuron *, ReceptiveFieldNeuronRole)> neuronPredicate,
							   function<void (Synapse *)> synapseCreated )
{
	currentSize.a = max( currentSize.a, _neuronSpacing.a );
	currentSize.b = max( currentSize.b, _neuronSpacing.b );

	fieldSize.a = max( fieldSize.a, other->_neuronSpacing.a );
	fieldSize.b = max( fieldSize.b, other->_neuronSpacing.b );

	ReceptiveFieldNeuronRole currentNeuronRole;
	ReceptiveFieldNeuronRole otherNeuronRole;
	switch( role )
	{
	case Source:
		currentNeuronRole = To;
		otherNeuronRole = From;
		break;
	case Target:
		currentNeuronRole = From;
		otherNeuronRole = To;
		break;
	default:
		assert( false );
	}

	NeuronSubset currentNeurons = findNeurons( currentCenter, currentSize );
	NeuronSubset otherNeurons = other->findNeurons( otherCenter, otherSize );

	// Iterate over each neuron in this sheet's region.
	for( Vector2i currentNeuronIndex : currentNeurons )
	{
		Neuron *currentNeuron = getNeuron( currentNeuronIndex );

		if( !neuronPredicate(currentNeuron, currentNeuronRole) ) 
			continue;
		
		
		// ---
		// --- Find all neurons within the receptive field.
		// ---
		NeuronSubset allReceptiveFieldNeurons =
			other->findReceptiveFieldNeurons( currentNeuron->sheetPosition,
											  fieldOffset,
											  fieldSize );

		// ---
		// --- Constrain by other center/size
		// ---
		NeuronSubset constrainedReceptiveFieldNeurons;
		constrainedReceptiveFieldNeurons._begin.a = max( allReceptiveFieldNeurons._begin.a,
														 otherNeurons._begin.a );
		constrainedReceptiveFieldNeurons._begin.b = max( allReceptiveFieldNeurons._begin.b,
														 otherNeurons._begin.b );
		constrainedReceptiveFieldNeurons._end.a = min( allReceptiveFieldNeurons._end.a,
													   otherNeurons._end.a );
		constrainedReceptiveFieldNeurons._end.b = min( allReceptiveFieldNeurons._end.b,
													   otherNeurons._end.b );

		//cout << "otherCount=" << other->_neuronCount << ", otherCenter=" << otherCenter << ", otherSize=" << otherSize << ", otherNeurons=" << otherNeurons << ", allRF=" << allReceptiveFieldNeurons << ", constrainRF=" << constrainedReceptiveFieldNeurons << ", size=" << constrainedReceptiveFieldNeurons.size() << endl;

		if( constrainedReceptiveFieldNeurons.size() < 1 )
			continue;
		
		// ---
		// --- Find which neurons pass the predicate, and determine the mean distance.
		// ---

		Neuron *fieldNeurons[ constrainedReceptiveFieldNeurons.size() ];
		int nfieldNeurons = 0;
		float totalDistance = 0;
		
		for( Vector2i otherNeuronIndex : constrainedReceptiveFieldNeurons )
		{
			Neuron *otherNeuron = other->getNeuron( otherNeuronIndex );

			// Ignore self-synapse
			if( currentNeuron == otherNeuron )
				continue;

			if( !neuronPredicate(otherNeuron, otherNeuronRole) ) 
				continue;
			
			fieldNeurons[ nfieldNeurons++ ] = otherNeuron;
			totalDistance += currentNeuron->absPosition.distance( otherNeuron->absPosition );
		}

		// No neurons in receptive field.
		if( nfieldNeurons == 0 )
			continue;

		// ---
		// --- Determine number of synapses to create
		// ---
		float meanDistance = totalDistance / nfieldNeurons;
		float probabilitySynapse = _sheetsModel->getProbabilitySynapse( meanDistance );
		int nsynapses = round( probabilitySynapse * nfieldNeurons );

		if( nsynapses == 0 )
			continue;

		// ---
		// --- Create synapses
		// ---
		int stride = nfieldNeurons / nsynapses;

		for( int offset = 0; nsynapses > 0; offset += stride, nsynapses-- )
		{
			Neuron *otherNeuron = fieldNeurons[offset];
			Synapse *synapse = NULL;

			switch( role )
			{
			case Source:
				synapse = createSynapse( otherNeuron, currentNeuron );

				//Should never have synapses originating at a modulatory neuron
				if (otherNeuron->attrs.type > 2){
					assert(false);
				}
				
				break;
			case Target:
				synapse = createSynapse( currentNeuron, otherNeuron );
				
				//Should never have synapses originating at a modulatory neuron
				if (currentNeuron->attrs.type > 2) {
					assert(false);
				}
				
				break;
			default:
				assert( false );
			}

			if( synapse )
				synapseCreated( synapse );
		}
	}
}


vector<Neuron*> Sheet::getNeuronsInGasRange(Neuron *origin) {
	
		float radius  = origin->attrs.emissionRadius;
		float receptorStrengthLowerBound  = 0.1; //TODO -  GLOBAL CONSTANT for a lower bound on receptorStrength
		
		//used simply to access the distance between sheets - any neuron will do, since we just want it's planePos
		Neuron *target  = getNeuron(0,0);
		float neuronsDistanceFromSheet = 0;

		switch( _orientation )
            { 
            case PlaneXY: 
            	neuronsDistanceFromSheet = target->absPosition.z - origin->absPosition.z;
            	break;
            case PlaneXZ: 
            	neuronsDistanceFromSheet = target->absPosition.y - origin->absPosition.y;
            	break;
            case PlaneZY: 
            	neuronsDistanceFromSheet = target->absPosition.x - origin->absPosition.x;
            	break;
            default: 
            	assert( false );
            }		
		
		vector<Neuron*> inRangeNeurons;
		
		//only look at sheets that could have neurons within range of the originating gas-emitting neuron
		if ( fabs(neuronsDistanceFromSheet) < radius) {
		
			for (int i = 0; i < _nneurons; i++) {
				Neuron *otherNeuron = &_neurons[i];
				
				if (origin->absPosition.distance( otherNeuron->absPosition) < radius && origin != otherNeuron) {
					if (otherNeuron->attrs.receptorStrength >= receptorStrengthLowerBound ) {
						//meets all conditions to be included
						inRangeNeurons.push_back( otherNeuron );
						if (Brain::config.gasnetsDebugMode > 3) {
					        cout << "     GasnetsDebugMode[4]: " ;
							cout << "        radius: " << radius << "\n";
							cout << "        origin->absPosition    : " << origin->absPosition << "\n";
							cout << "        otherNeuron.absPosition: " << otherNeuron->absPosition << "\n";
							cout << "        distance between: " << origin->absPosition.distance( otherNeuron->absPosition) << "\n";
						}
					}
				}
			}
			
			if (Brain::config.gasnetsDebugMode > 3) {
				cout << "     GasnetsDebugMode[4]: " ;
				cout << "Number of gas neurons in sheet within range: " << inRangeNeurons.size() << "\n";
			}
		}
		return inRangeNeurons;
}



//---------------------------------------------------------------------------
// Sheet::addGasChannels
//---------------------------------------------------------------------------
int Sheet::addGasChannels( Neuron *origin,
							   SheetsModel *model)
{

	vector<Neuron*> inRangeNeuronsInAllSheets;
	for (int i=0; i < model->_allSheets.size(); i++) {
			
		//this was the easiest data structure to use to do this, maybe there is a better way? //todo fixme
		vector<Neuron*> inRangeNeurons = model->_allSheets[i]->getNeuronsInGasRange(origin);
		
		// Iterate over each neuron in range of this sheet's region.
        for (int j=0; j< inRangeNeurons.size(); j++)
		{
			Neuron* currentNeuron = inRangeNeurons[j];
			inRangeNeuronsInAllSheets.push_back(currentNeuron); //push pointer in vector of pointers
		}
		
	}
	
	if (Brain::config.gasnetsDebugMode > 2) {
        cout << "   GasnetsDebugMode[3]: " ;
		cout << "  Number of Neurons in Range:" <<inRangeNeuronsInAllSheets.size() << "\n";
	}
	
	for (int i=0; i< inRangeNeuronsInAllSheets.size(); i++) {
	    Neuron* otherNeuron = inRangeNeuronsInAllSheets[i];
	    createGasChannel( origin, otherNeuron );
	}
	
	return inRangeNeuronsInAllSheets.size();

}


//---------------------------------------------------------------------------
// Sheet::trim
//---------------------------------------------------------------------------
void Sheet::trim( float &center, float &size )
{
	float overflow;

	overflow = (center + size/2) - 1.0;
	if( overflow > 0 )
	{
		size -= overflow;
		center -= overflow / 2;

		assert( fabs(center + (size/2) - 1.0) < 1e-6 );
	}
	else
	{
		overflow = 0 - (center - size/2);
		if( overflow > 0 )
		{
			size -= overflow;
			center += overflow / 2;

			assert( fabs(center - (size/2)) < 1e-6 );
		}
	}
}

//---------------------------------------------------------------------------
// Sheet::createNeurons
//---------------------------------------------------------------------------
void Sheet::createNeurons( function<void (Neuron*)> &neuronCreated )
{
	_neurons = new Neuron[ _nneurons ];

	for( int i = 0; i < _neuronCount.a; i++ )
	{
		for( int j = 0; j < _neuronCount.b; j++ )
		{
			Neuron *neuron = getNeuron( i, j );

			neuron->sheet = this;
			neuron->nonCulledId = _sheetsModel->_numNonCulledNeurons++;
			neuron->id = -1;
			neuron->sheetIndex.set( i, j );
			neuron->cullState.touchedFromInput = false;
			neuron->cullState.touchedFromOutput = false;
			
			// ---
			// --- Compute position within sheet 
			// ---

			// 2D position within sheet
			neuron->sheetPosition.set( _neuronInsets.a + (i * _neuronSpacing.a),
									   _neuronInsets.b + (j * _neuronSpacing.b) );

			// ---
			// --- Compute absolute 3D position within model
			// ---


			// Scale 2D position by sheet size
			Vector2f position = neuron->sheetPosition;
			
			position.scale( _size );

			// Translate by sheet location
			position += _center - ( _size / 2 );

			// Create 3D position
			neuron->absPosition.set( _orientation, position, _slot );


			// Scale to model size.
			neuron->absPosition.scale( _sheetsModel->getSize() );
			
			trc( "Neuron[" << i << "," << j << "]: sheetPos=" << neuron->sheetPosition << ", absPos=" << neuron->absPosition );

			assert( neuron == getNeuron(neuron->sheetIndex) );

			neuronCreated( neuron );
			
		}
	}
}

//---------------------------------------------------------------------------
// Sheet::findNeurons
//---------------------------------------------------------------------------
NeuronSubset Sheet::findNeurons( const Vector2f &center,
								 const Vector2f &size )
{
	NeuronSubset result;

	Vector2f ul = center - (size / 2);
	Vector2f lr = center + (size / 2);

	result._begin.a = max( 0, (int)ceilf( (ul.a - _neuronInsets.a) / _neuronSpacing.a ) );
	result._end.a = min( _neuronCount.a - 1, (int)floorf( (lr.a - _neuronInsets.a) / _neuronSpacing.a ) );
	result._begin.b = max( 0, (int)ceilf( (ul.b - _neuronInsets.b) / _neuronSpacing.b ) );
	result._end.b =  min( _neuronCount.b - 1, (int)floorf( (lr.b - _neuronInsets.b) / _neuronSpacing.b ) );
	
	return result;
}

//---------------------------------------------------------------------------
// Sheet::findReceptiveFieldNeurons
//---------------------------------------------------------------------------
NeuronSubset Sheet::findReceptiveFieldNeurons( const Vector2f &neuronPosition,
											   const Vector2f &fieldOffset,
											   const Vector2f &fieldSize )
{	
	struct local
	{
		static void offsetCenter( float &center, float offset )
		{
			if( offset < 0 )
				center += center * offset;
			else
				center += (1 - center) * offset;
		}
	};

	vtrc( "Receptive neuron pos=" << neuronPosition << ", offset=" << fieldOffset << ", size=" << fieldSize );

	Vector2f center = neuronPosition;

	local::offsetCenter( center.a, fieldOffset.a );
	local::offsetCenter( center.b, fieldOffset.b );

	vtrc( "  Adjusted center=" << center );

	return findNeurons( center, fieldSize );
}


//---------------------------------------------------------------------------
// Sheet::createSynapse
//---------------------------------------------------------------------------
Synapse *Sheet::createSynapse( Neuron *from, Neuron *to )
{
	if( from == to )
		return NULL;

	if( from->synapsesOut[to] != NULL )
		return NULL;

	Synapse *synapse = new Synapse();
	synapse->from = from;
	synapse->to = to;

	from->synapsesOut[to] = synapse;
	to->synapsesIn[from] = synapse;

	trc( "Synapse [" << from->sheet->_id << "] " << from->absPosition << " --> [" << to->sheet->_id << "] " << to->absPosition );

	return synapse;
}

//---------------------------------------------------------------------------
// Sheet::createGasChannel
//---------------------------------------------------------------------------
GasChannel *Sheet::createGasChannel( Neuron *from, Neuron *to )
{
	if( from == to )
		return NULL;

	if( from->gasChannelsOut[to] != NULL )
		return NULL;

		
	GasChannel *gasChannel = new GasChannel();
	
	gasChannel->from = from;
	gasChannel->to = to;

	from->gasChannelsOut[to] = gasChannel;
	to->gasChannelsIn[from] = gasChannel;  

	trc( "GasChannel [" << from->sheet->_id << "] " << from->absPosition << " --> [" << to->sheet->_id << "] " << to->absPosition );

	return gasChannel;
}




//===========================================================================
// SheetsModel
//===========================================================================

//---------------------------------------------------------------------------
// SheetsModel::SheetsModel
//---------------------------------------------------------------------------
SheetsModel::SheetsModel( const Vector3f &size, float synapseProbabilityX )
: _numNonCulledNeurons( 0 )
, _size( size )
, _synapseProbabilityX( synapseProbabilityX )
{
}

//---------------------------------------------------------------------------
// SheetsModel::~SheetsModel
//---------------------------------------------------------------------------
SheetsModel::~SheetsModel()
{
	itfor( SheetVector, _allSheets, it )
		delete *it;
}

//---------------------------------------------------------------------------
// SheetsModel::createSheet
//---------------------------------------------------------------------------
Sheet *SheetsModel::createSheet( string name,
								 int id,
								 Sheet::Type type,
								 Orientation orientation,
								 float slot,
								 const Vector2f &center,
								 const Vector2f &size,
								 const Vector2i &neuronCount,
								 function<void (Neuron*)> &neuronCreated )
{
	if( id == -1 )
		id = (int)_allSheets.size();

	Sheet *sheet = new Sheet( this,
							  name,
							  id,
							  type,
							  orientation,
							  slot,
							  center,
							  size,
							  neuronCount,
							  neuronCreated );

	if( (int)_allSheets.size() <= id )
		_allSheets.resize( id + 1 );
	assert( _allSheets[id] == NULL );
	_allSheets[id] = sheet;

	switch( type )
	{
	case Sheet::Input:
		_inputSheets.push_back( sheet );
		break;
	case Sheet::Output:
		_outputSheets.push_back( sheet );
		break;
	case Sheet::Internal:
		_internalSheets.push_back( sheet );
		break;
	default:
		assert( false );
		break;
	}

	return sheet;
}

//---------------------------------------------------------------------------
// SheetsModel::getSheet
//---------------------------------------------------------------------------
Sheet *SheetsModel::getSheet( int id )
{
	if( id >= (int)_allSheets.size() )
		return NULL;

	return _allSheets[ id ];
}

//---------------------------------------------------------------------------
// SheetsModel::getSheets
//---------------------------------------------------------------------------
const SheetVector &SheetsModel::getSheets( Sheet::Type type )
{
	switch( type )
	{
	case Sheet::Input: return _inputSheets;
	case Sheet::Output: return _outputSheets;
	case Sheet::Internal: return _internalSheets;
	default: assert( false );
	}
}

//---------------------------------------------------------------------------
// SheetsModel::cull
//---------------------------------------------------------------------------
void SheetsModel::cull()
{
	assert( !_inputSheets.empty() );
	assert( !_outputSheets.empty() );

	for( Sheet *sheet : _inputSheets )
	{
		const Vector2i count = sheet->getNeuronCount();
		for( int a = 0; a < count.a; a++ )
			for( int b = 0; b < count.b; b++ )
				touchFromInput( sheet->getNeuron(a, b) );
	}

	for( Sheet *sheet : _outputSheets )
	{
		const Vector2i count = sheet->getNeuronCount();
		for( int a = 0; a < count.a; a++ )
			for( int b = 0; b < count.b; b++ )
				touchFromOutput( sheet->getNeuron(a, b) );
	}

	addNonCulledNeurons( _inputSheets );
	addNonCulledNeurons( _outputSheets );
	addNonCulledNeurons( _internalSheets );
}

//---------------------------------------------------------------------------
// SheetsModel::getNeurons
//---------------------------------------------------------------------------
NeuronVector &SheetsModel::getNeurons()
{
	return _neurons;
}

//---------------------------------------------------------------------------
// SheetsModel::getProbabilitySynapse
//---------------------------------------------------------------------------
float SheetsModel::getProbabilitySynapse( float distance )
{
	if( _synapseProbabilityX == 0.0 )
		return 1.0;
	else		
		return (1 / _synapseProbabilityX) * exp(-distance/_synapseProbabilityX);
}

//---------------------------------------------------------------------------
// SheetsModel::touchFromInput
//---------------------------------------------------------------------------
void SheetsModel::touchFromInput( Neuron *neuron )
{
	if( neuron->cullState.touchedFromInput )
		return;

	neuron->cullState.touchedFromInput = true;

	itfor( SynapseMap, neuron->synapsesOut, it )
	{
		Synapse *syn = it->second;
		touchFromInput( syn->to );
	}
}

//---------------------------------------------------------------------------
// SheetsModel::touchFromOutput
//---------------------------------------------------------------------------
void SheetsModel::touchFromOutput( Neuron *neuron )
{
	if( neuron->cullState.touchedFromOutput )
		return;

	neuron->cullState.touchedFromOutput = true;

	itfor( SynapseMap, neuron->synapsesIn, it )
	{
		Synapse *syn = it->second;
		touchFromOutput( syn->from );
	}
}


//---------------------------------------------------------------------------
// SheetsModel::removeNeuronSynapsesAndGasChannels
//---------------------------------------------------------------------------
void SheetsModel::removeNeuronSynapsesAndGasChannels( Neuron *neuron) 
{
	
	trc( "CULLED: " << "[" << sheet->getId() << "] " << neuron->absPosition );

	itfor( SynapseMap, neuron->synapsesIn, it_syn )
	{
		Synapse *syn = it_syn->second;
		syn->from->synapsesOut.erase( neuron );
		delete syn;
	}
	neuron->synapsesIn.clear();

	itfor( SynapseMap, neuron->synapsesOut, it_syn )
	{
		Synapse *syn = it_syn->second;
		syn->to->synapsesIn.erase( neuron );
		delete syn;
	}
	neuron->synapsesOut.clear();

	//GASCHANNEL DELETIONS
	itfor( GasChannelMap, neuron->gasChannelsIn, it_gc )
	{
		GasChannel *gc = it_gc->second;
		gc->from->gasChannelsOut.erase( neuron );
		delete gc;
	}
	neuron->gasChannelsIn.clear();

	itfor( GasChannelMap, neuron->gasChannelsOut, it_gc )
	{
		GasChannel *gc = it_gc->second;
		gc->to->gasChannelsIn.erase( neuron );
		delete gc;
	}
	neuron->gasChannelsOut.clear();

}





//---------------------------------------------------------------------------
// SheetsModel::addNonCulledNeurons
//---------------------------------------------------------------------------
void SheetsModel::addNonCulledNeurons( SheetVector &sheets )
{
	itfor( SheetVector, sheets, it )
	{
		Sheet *sheet = *it;
		const Vector2i count = sheet->getNeuronCount();
		for( int a = 0; a < count.a; a++ )
		{
			for( int b = 0; b < count.b; b++ )
			{
				Neuron *neuron = sheet->getNeuron( a, b );
				
				//Unfortunately we do not know the gasChannel sizes (hence the culling on synapses)
				//So we cannot necessarily cull every single possible neuron that might not be meaningful
				//TODO - find a way to cull all possible neurons by using the incoming/outgoing gaschannels
				if( neuron->cullState.touchedFromOutput && neuron->cullState.touchedFromInput )
				{
					neuron->id = (int)_neurons.size();
					_neurons.push_back( neuron );
					
				}
				else //at least one of the input/output was not culled
				{
					if (Brain::config.gasnetsEnabled)  //gasnets enabled
					{ 
						bool isGasNeuron = neuron->attrs.type > 2; //irrelevant
						bool isGasActivated = neuron->attrs.activatedByGas > 0.0;
						
						//if touched by input we can go ahead and add -no matter what
						if (neuron->cullState.touchedFromInput) {
						
							neuron->id = (int)_neurons.size();
							_neurons.push_back( neuron );
						
						//if not touched by input, but potentially activated by gas, we are still potentially OK as long as we have incoming GCs
						} else if ( isGasActivated  ) {
							
							neuron->id = (int)_neurons.size();
							_neurons.push_back( neuron );
							
						//if we are not touched by input and not potentially gas activated, we can never get activated!
						} else {
						
							removeNeuronSynapsesAndGasChannels( neuron);	
						}
						
					} else {
						//normal operation of code
						removeNeuronSynapsesAndGasChannels(neuron);
					}
				}
			}
		}
	}

}
