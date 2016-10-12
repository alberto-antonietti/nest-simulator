/*
 *  stdp_connection_sinexp.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
 
 /*
 
   Alberto Antonietti
   alberto.antonietti@polimi.it
  
   Cerebellar PF-PC Plasticity with an exp. sin. Kernel LTP and LTD
 
 */

#include "stdp_connection_sinexp.h"

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connector_model.h"
#include "event.h"
#include "kernel_manager.h"

// Includes from sli:
#include "dictdatum.h"


namespace nest
{
//
// Implementation of class STDPSinExpCommonProperties.
//

STDPSinExpCommonProperties::STDPSinExpCommonProperties()
  : CommonSynapseProperties()
  , vt_( 0 )
  , A_plus_( 1.0 )
  , A_minus_( 1.5 )
  , Wmin_( 0.0 )
  , Wmax_( 200.0 )
{
}

void STDPSinExpCommonProperties::get_status( DictionaryDatum& d ) const{
  CommonSynapseProperties::get_status( d );

  if ( vt_ != 0 )
    def< long_t >( d, "modulator", vt_->get_gid() );
  else
    def< long_t >( d, "modulator", -1 );

  def< double_t >( d, "A_plus", A_plus_ );
  def< double_t >( d, "A_minus", A_minus_ );
  def< double_t >( d, "Wmin", Wmin_ );
  def< double_t >( d, "Wmax", Wmax_ );
}

void STDPSinExpCommonProperties::set_status( const DictionaryDatum& d, ConnectorModel& cm ){
  CommonSynapseProperties::set_status( d, cm );

  long_t vtgid;
  if ( updateValue< long_t >( d, "vt", vtgid ) ){
    vt_ = dynamic_cast< volume_transmitter_alberto* >( kernel().node_manager.get_node( vtgid ) );

      

    if ( vt_ == 0 )
      throw BadProperty( "vt needs to be a Volume Transmitter" );
  }

  updateValue< double_t >( d, "A_plus", A_plus_ );
  updateValue< double_t >( d, "A_minus", A_minus_ );
  updateValue< double_t >( d, "Wmin", Wmin_ );
  updateValue< double_t >( d, "Wmax", Wmax_ );
}

Node* STDPSinExpCommonProperties::get_node(){
  if ( vt_ == 0 )
    throw BadProperty( "No neuron has been assigned as the modulator of the synapse." );
  else
    return vt_;
}

} // End of namespace nest
