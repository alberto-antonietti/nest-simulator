/*
 *  connector_model.h
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

#ifndef CONNECTOR_MODEL_H
#define CONNECTOR_MODEL_H

// C++ includes:
#include <cmath>
#include <string>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "event.h"
#include "nest_time.h"
#include "nest_types.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{
class ConnectorBase;
class CommonSynapseProperties;
class TimeConverter;
class Node;


inline ConnectorBase*
pack_pointer( ConnectorBase* p, bool has_primary, bool has_secondary )
{
  return reinterpret_cast< ConnectorBase* >(
    reinterpret_cast< unsigned long >( p ) | has_primary | ( has_secondary << 1 ) );
}

inline ConnectorBase*
validate_pointer( ConnectorBase* p )
{
  // erase 2 least significant bits to obtain the correct pointer
  return reinterpret_cast< ConnectorBase* >(
    ( reinterpret_cast< unsigned long >( p ) & ( -1l - 3l ) ) );
}

inline bool
has_primary( ConnectorBase* p )
{
  // the lowest bit is set, if there is at least one primary
  // connection
  return static_cast< bool >( reinterpret_cast< unsigned long >( p ) & 1 );
}
inline bool

has_secondary( ConnectorBase* p )
{
  // the second lowest bit is set, if there is at least one secondary
  // connection
  return static_cast< bool >( reinterpret_cast< unsigned long >( p ) & 2 );
}


class ConnectorModel
{

public:
  ConnectorModel( const std::string, bool is_primary, bool has_delay );
  ConnectorModel( const ConnectorModel&, const std::string );
  virtual ~ConnectorModel()
  {
  }

  /**
   * NAN is a special value in cmath, which describes double values that
   * are not a number. If delay or weight is omitted in an add_connection call,
   * NAN indicates this and weight/delay are set only, if they are valid.
   */
  virtual ConnectorBase* add_connection( Node& src,
    Node& tgt,
    ConnectorBase* conn,
    synindex syn_id,
    double_t delay = numerics::nan,
    double_t weight = numerics::nan ) = 0;
  virtual ConnectorBase* add_connection( Node& src,
    Node& tgt,
    ConnectorBase* conn,
    synindex syn_id,
    DictionaryDatum& d,
    double_t delay = numerics::nan,
    double_t weight = numerics::nan ) = 0;

  virtual ConnectorModel* clone( std::string ) const = 0;

  virtual void calibrate( const TimeConverter& tc ) = 0;

  virtual void get_status( DictionaryDatum& ) const = 0;
  virtual void set_status( const DictionaryDatum& ) = 0;

  virtual const CommonSynapseProperties& get_common_properties() const = 0;

  virtual SecondaryEvent* get_event() const = 0;

  virtual void set_syn_id( synindex syn_id ) = 0;

  virtual std::vector< SecondaryEvent* > create_event( size_t n ) const = 0;

  std::string
  get_name() const
  {
    return name_;
  }

  bool
  get_user_set_delay_extrema() const
  {
    return user_set_delay_extrema_;
  }

  Network&
  network() const
  {
    return net_;
  }

  bool
  is_primary() const
  {
    return is_primary_;
  }

  bool
  has_delay() const
  {
    return has_delay_;
  }

protected:
  std::string name_;
  bool default_delay_needs_check_; //!< Flag indicating, that the default delay must be checked
  bool is_primary_; //!< indicates, whether this ConnectorModel belongs to a primary connection
  bool has_delay_;  //!< indicates, that ConnectorModel has a delay

}; // ConnectorModel


template < typename ConnectionT >
class GenericConnectorModel : public ConnectorModel
{
private:
  typename ConnectionT::CommonPropertiesType cp_;
  typename ConnectionT::EventType*
    pev_; //!< used to create secondary events that belong to secondary connections

  ConnectionT default_connection_;
  rport receptor_type_;

public:
  GenericConnectorModel( const std::string name, bool is_primary, bool has_delay )
    : ConnectorModel( name, is_primary, has_delay )
    , receptor_type_( 0 )
  {
  }

  GenericConnectorModel( const GenericConnectorModel& cm, const std::string name )
    : ConnectorModel( cm, name )
    , cp_( cm.cp_ )
    , pev_( cm.pev_ )
    , default_connection_( cm.default_connection_ )
    , receptor_type_( cm.receptor_type_ )
  {
  }

  ConnectorBase* add_connection( Node& src,
    Node& tgt,
    ConnectorBase* conn,
    synindex syn_id,
    double_t weight,
    double_t delay );
  ConnectorBase* add_connection( Node& src,
    Node& tgt,
    ConnectorBase* conn,
    synindex syn_id,
    DictionaryDatum& d,
    double_t weight,
    double_t delay );

  ConnectorModel* clone( std::string ) const;

  void calibrate( const TimeConverter& tc );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

  typename ConnectionT::CommonPropertiesType const&
  get_common_properties() const
  {
    return cp_;
  }

  void set_syn_id( synindex syn_id );

  virtual typename ConnectionT::EventType*
  get_event() const
  {
    assert( false );
    return 0;
  }

  ConnectionT const&
  get_default_connection() const
  {
    return default_connection_;
  }

  virtual std::vector< SecondaryEvent* >
  create_event( size_t n ) const
  {
    // Should not be called for a ConnectorModel belonging to a primary
    // connection. Only required for secondary connection types.
    assert( false );
    std::vector< SecondaryEvent* > prototype_events;
    return prototype_events;
  }

private:
  void used_default_delay();

  ConnectorBase* add_connection( Node& src,
    Node& tgt,
    ConnectorBase* conn,
    synindex syn_id,
    ConnectionT& c,
    rport receptor_type );

}; // GenericConnectorModel

template < typename ConnectionT >
class GenericSecondaryConnectorModel : public GenericConnectorModel< ConnectionT >
{
private:
  typename ConnectionT::EventType*
    pev_; //!< used to create secondary events that belong to secondary connections

public:
  GenericSecondaryConnectorModel( const std::string name, bool has_delay )
    : GenericConnectorModel< ConnectionT >( name, /*is _primary=*/false, has_delay )
    , pev_( 0 )
  {
    pev_ = new typename ConnectionT::EventType();
  }

  GenericSecondaryConnectorModel( const GenericSecondaryConnectorModel& cm, const std::string name )
    : GenericConnectorModel< ConnectionT >( cm, name )
  {
    pev_ = new typename ConnectionT::EventType( *cm.pev_ );
  }


  ConnectorModel*
  clone( std::string name ) const
  {
    return new GenericSecondaryConnectorModel( *this, name ); // calls copy construtor
  }

  std::vector< SecondaryEvent* >
  create_event( size_t n ) const
  {
    std::vector< SecondaryEvent* > prototype_events( n, NULL );
    for ( size_t i = 0; i < n; i++ )
      prototype_events[ i ] = new typename ConnectionT::EventType();

    return prototype_events;
  }


  ~GenericSecondaryConnectorModel()
  {
    if ( pev_ != 0 )
      delete pev_;
  }

  typename ConnectionT::EventType*
  get_event() const
  {
    return pev_;
  }
};

} // namespace nest

#endif /* #ifndef CONNECTOR_MODEL_H */
