/*
 *  stdp_connection_sinexp.h
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

#ifndef STDP_CONNECTION_SINEXP_H
#define STDP_CONNECTION_SINEXP_H

/* BeginDocumentation

   Name: 

   Description:
   
   Examples:
   
   Parameters:
           vt		 long   - ID of volume_transmitter collecting the spikes from the pool of
                              dopamine releasing neurons and transmitting the spikes
                              to the synapse. A value of -1 indicates that no volume
                              transmitter has been assigned.
     Common properties:
           A_plus    double - Amplitude of weight change for facilitation
           A_minus   double - Amplitude of weight change for depression
           Wmin      double - Minimal synaptic weight
           Wmax      double - Maximal synaptic weight

*/

#include "connection.h"
#include "spikecounter.h"
#include "volume_transmitter_alberto.h"
#include "numerics.h"
#include <math.h>

namespace nest{

/**
 * Class containing the common properties for all synapses of type dopamine connection.
 */
class STDPSinExpCommonProperties : public CommonSynapseProperties{
public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  STDPSinExpCommonProperties();

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );
  
  long get_vt_gid() const;

  double A_plus_;
  double A_minus_;
  double Wmin_;
  double Wmax_;
  volume_transmitter_alberto* vtC_;
};

inline long
STDPSinExpCommonProperties::get_vt_gid() const
{
    return -2;
}



/**
 * Class representing an STDPSinExpConnection with homogeneous parameters,
 * i.e. parameters are the same for all synapses.
 */
template < typename targetidentifierT > class STDPSinExpConnection : public Connection< targetidentifierT > {

public:

  Node* get_node();

  long get_vt_gid() const;

  volume_transmitter_alberto* vt_;
    
  std::vector<double> SpikeBuffer_;
  
  typedef STDPSinExpCommonProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  STDPSinExpConnection();

  STDPSinExpConnection( const STDPSinExpConnection& );

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay;
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void send( Event& e, thread t, double, const STDPSinExpCommonProperties& cp );

  void trigger_update_weight( thread t, const std::vector< spikecounter >& dopa_spikes, double t_trig, const STDPSinExpCommonProperties& cp);
    
  class ConnTestDummyNode : public ConnTestDummyNodeBase{
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport ){
      return invalid_port_;
    }
  };
  

  /*
   * This function calls check_connection on the sender and checks if the receiver
   * accepts the event type and receptor type requested by the sender.
   * Node::check_connection() will either confirm the receiver port by returning
   * true or false if the connection should be ignored.
   * We have to override the base class' implementation, since for STDP
   * connections we have to call register_stdp_pl_connection on the target neuron
   * to inform the Archiver to collect spikes for this connection.
   * Further, the STDP dopamine synapse requires a volume transmitter to be set before
   * any simulation is performed. Checking this satisfies ticket #926.
   *
   * \param s The source node
   * \param r The target node
   * \param receptor_type The ID of the requested receptor type
   * \param t_lastspike last spike produced by presynaptic neuron (in ms)
   */
  void check_connection( Node& s, Node& t, rport receptor_type, double t_lastspike, const CommonPropertiesType& cp ){
	  
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
  }

  void set_weight( double w ){
    weight_ = w;
  }

private:
  // update dopamine trace from last to current dopamine spike and increment index
  void update_dopamine_( const std::vector< spikecounter >& dopa_spikes,const STDPSinExpCommonProperties& cp );

  void update_weight_(double weight_change, const STDPSinExpCommonProperties& cp );
  
  void process_dopa_spikes_( const std::vector< spikecounter >& dopa_spikes, double t0, double t1, const STDPSinExpCommonProperties& cp );

  // data members of each connection
  double weight_;

  // dopa_spikes_idx_ refers to the dopamine spike that has just been processes
  // after trigger_update_weight a pseudo dopamine spike at t_trig is stored at index 0 and
  // dopa_spike_idx_ = 0
  index dopa_spikes_idx_;

  // time of last update, which is either time of last presyn. spike or time-driven update
  double t_last_update_;
};

//
// Implementation of class STDPSinExpConnection.
//

template < typename targetidentifierT > STDPSinExpConnection< targetidentifierT >::STDPSinExpConnection()
  : ConnectionBase()
  , vt_ ( 0 )
  , weight_( 1.0 )
  , dopa_spikes_idx_( 0 )
  , t_last_update_( 0.0 )
{
}

template < typename targetidentifierT > STDPSinExpConnection< targetidentifierT >::STDPSinExpConnection( const STDPSinExpConnection& rhs )
  : ConnectionBase( rhs )
  , vt_ ( rhs.vt_ )
  , weight_( rhs.weight_ )
  , dopa_spikes_idx_( rhs.dopa_spikes_idx_ )
  , t_last_update_( rhs.t_last_update_ )
{
}

template < typename targetidentifierT > void STDPSinExpConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const{

  // base class properties, different for individual synapse
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  if ( vt_ != 0 )
    def< long >( d, "modulator", vt_->get_gid() );
  else
    def< long >( d, "modulator", -1 );

}

template < typename targetidentifierT > long STDPSinExpConnection< targetidentifierT >::get_vt_gid( ) const{

    if ( vt_ != 0 ){
		return vt_->get_gid();
	}
	else
		return -1;

}

template < typename targetidentifierT > void STDPSinExpConnection< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm ){
  // base class properties
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  long vtgid;
  if ( updateValue< long >( d, "vt", vtgid ) ){
    vt_ = dynamic_cast< volume_transmitter_alberto* >( kernel().node_manager.get_node( vtgid ) );
    if ( vt_ == 0 )
      throw BadProperty( "vt needs to be a Volume Transmitter" );
  }
}


template < typename targetidentifierT > inline void STDPSinExpConnection< targetidentifierT >::update_dopamine_(
															  const std::vector< spikecounter >& dopa_spikes,
															  const STDPSinExpCommonProperties& cp ){
	// We enter here when there is a spike of the Volume Transmitter	
													  
	double minus_dt = dopa_spikes[ dopa_spikes_idx_+1].spike_time_-1;
   // std::cout << "We enter here when there is a spike of the Volume Transmitter"	<< std::endl;
	if (SpikeBuffer_.size()>0){
		//std::cout << SpikeBuffer_[0]  << " <-SpikeBuffer_[0] \t minus_dt-> " << minus_dt << std::endl;
		double LTD_amount = 0.0;
		for(int GR = 0; GR<SpikeBuffer_.size(); GR++){
			double sd= SpikeBuffer_[GR] - minus_dt;
			//std::cout << "SD: " << sd << std::endl;
			if (sd < -200){
				//SpikeBuffer_.erase(SpikeBuffer_.begin()+GR);
				//GR--;
			}
			if ( sd<0 && sd>=-200){
				LTD_amount += cp.A_minus_ * exp(-(sd-150)/1000.0)*pow((sin(2*3.1415*(sd-150)/1000.0)),20)/1.2848;
				//std::cout << "spike distance: " << sd << " LTD amout : " << exp(-(sd-150)/1000.0)*pow((sin(2*3.1415*(sd-150)/1000.0)),20)/1.2848 << std::endl;
			}
		}
		update_weight_(LTD_amount, cp);
	}
	
  ++dopa_spikes_idx_;
}




template < typename targetidentifierT > inline void STDPSinExpConnection< targetidentifierT >::update_weight_(double weight_change, const STDPSinExpCommonProperties& cp ){
  // LTP or LTD, depending on who calls this function
  weight_ = weight_+weight_change;
  
  if ( weight_ < cp.Wmin_ )
    weight_ = cp.Wmin_;
  if ( weight_ > cp.Wmax_ )
    weight_ = cp.Wmax_;
}

template < typename targetidentifierT > inline void STDPSinExpConnection< targetidentifierT >::process_dopa_spikes_(const std::vector< spikecounter >& dopa_spikes, double t0, double t1, const STDPSinExpCommonProperties& cp ){
  // process dopa spikes in (t0, t1]
  // propagate weight from t0 to t1
  if ( ( dopa_spikes.size() > dopa_spikes_idx_ + 1 ) && ( dopa_spikes[ dopa_spikes_idx_ + 1 ].spike_time_ <= t1 ) ){
    // A IO SPIKE IS DETECTED AT TIME T0, LTD happens with a different amplitude, it depends on the distance between IO SPIKE and PF spikes
    update_dopamine_( dopa_spikes, cp );
  }
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
template < typename targetidentifierT > inline void STDPSinExpConnection< targetidentifierT >::send( Event& e, thread t, double, const STDPSinExpCommonProperties& cp ){
  // t_lastspike_ = 0 initially

  Node* target = get_target( t );
  
  double t_spike = e.get_stamp().get_ms();
  //std::cout << "SEND" << std::endl;
  
  // LTP (of a factor A_plus) due to new pre-synaptic spike
  double t_spike_d = t_spike;
  SpikeBuffer_.push_back(t_spike_d);
  update_weight_(cp.A_plus_, cp);
  while(SpikeBuffer_[0]<t_spike-200.0){
	  SpikeBuffer_.erase(SpikeBuffer_.begin());
	  }  
  e.set_receiver( *target );
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  t_last_update_ = t_spike;
}

template < typename targetidentifierT > inline void STDPSinExpConnection< targetidentifierT >::trigger_update_weight(
																   thread t,
																   const std::vector< spikecounter >& dopa_spikes,
																   const double t_trig,
																   const STDPSinExpCommonProperties& cp
																   ){
  
  int Vid_Check = (dopa_spikes.back()).multiplicity_;
  std::vector< spikecounter > dopa_temp = dopa_spikes;
  dopa_temp.pop_back();
  const std::vector< spikecounter > dopa_temp2 = dopa_temp;
  //std::cout << Vid_Check << "<Check - vid>" << get_vt_gid() << std::endl;
  if (Vid_Check != get_vt_gid())
	return;

  // purely dendritic delay
  double dendritic_delay = get_delay();

  // get spike history in relevant range (t_last_update, t_trig] from postsyn. neuron
  std::deque< histentry >::iterator start;
  std::deque< histentry >::iterator finish;
  get_target( t )->get_history(t_last_update_ - dendritic_delay, t_trig - dendritic_delay, &start, &finish );
  
  //std::cout << get_target( t )->get_gid() << std::endl;
  
  // facilitation due to postsyn. spikes since last update
  double t0 = t_last_update_;

  // propagate weight, eligibility trace c, dopamine trace n and facilitation trace K_plus to time
  // t_trig
  // but do not increment/decrement as there are no spikes to be handled at t_trig
  //std::cout << "trigger_update_weight" << std::endl;
  process_dopa_spikes_( dopa_temp2, t0, t_trig, cp );
  
  t_last_update_ = t_trig;
  dopa_spikes_idx_ = 0;
}


template < typename targetidentifierT > inline Node* STDPSinExpConnection< targetidentifierT >::get_node(){
  if ( vt_ == 0 )
    throw BadProperty( "No neuron has been assigned as the modulator of the synapse." );
  else
    return vt_;
}



} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_SINEXP_H