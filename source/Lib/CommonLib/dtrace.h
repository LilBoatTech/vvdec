/* -----------------------------------------------------------------------------
The copyright in this software is being made available under the Clear BSD
License, included below. No patent rights, trademark rights and/or 
other Intellectual Property Rights other than the copyrights concerning 
the Software are granted under this license.

The Clear BSD License

Copyright (c) 2018-2022, Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V. & The VVdeC Authors.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted (subject to the limitations in the disclaimer below) provided that
the following conditions are met:

     * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.

     * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

     * Neither the name of the copyright holder nor the names of its
     contributors may be used to endorse or promote products derived from this
     software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
THIS LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


------------------------------------------------------------------------------------------- */

/** \file     dtrace.h
 *  \brief    Implementation of trace messages support for debugging
 */

#pragma once

#include <stdio.h>

#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <cstdarg>

namespace vvdec
{

class CDTrace;

typedef std::string CType;

struct dtrace_channel
{
  int channel_number;
  std::string channel_name;
};

typedef std::vector<dtrace_channel> dtrace_channels_t;

class Condition
{
public:
    CType type;
    bool ( *eval )( int, int );
    int rval;

    Condition( CType t, bool ( *efunc )( int,int ), int refval )
    : type(t), eval(efunc), rval(refval)
    {}
};

class Channel
{
    typedef std::vector<Condition> Rule;
public:
    Channel() : rule_list(), _active(false), _counter(0) {}
    void update( std::map< CType, int > state );
    bool active() { return _active; }
    void add( Rule rule );
    void incrementCounter() { _counter++; }
    void decrementCounter() { _counter--  ; }
    int64_t getCounter() { return _counter; }
private:
    std::list< Rule > rule_list;
    bool _active;
    int64_t _counter;
};

class CDTrace
{
  typedef std::pair< CType, int > state_type;

  enum class Error
  {
    UnknownChannel = -3,
    BadRule        = -2,
    FileOpenFailed = -1,
    OK             = 0
  };


  //friend class Rules;
private:
  bool  copy         = false;
  FILE* m_trace_file = NULL;
  Error m_error_code = Error::OK;

  typedef std::string              Key;
  typedef std::vector<std::string> vstring;
  typedef std::map<Key, int>       channel_map_t;
  std::vector<Channel>             chanRules;
  std::set<CType>                  condition_types;
  std::map<CType, int>             state;
  std::map<Key, int>               deserializationTable;

public:
    CDTrace() : copy(false), m_trace_file(NULL) {}
    CDTrace( const char *filename, const vstring& channel_names );
    CDTrace( const char *filename, const dtrace_channels_t& channels );
    CDTrace( const std::string& sTracingFile, const std::string& sTracingRule, const dtrace_channels_t& channels );
    CDTrace( const CDTrace& other );
    CDTrace& operator=( const CDTrace& other );
    ~CDTrace();
    void swap         ( CDTrace& other );
    Error addRule     ( std::string rulestring );
    template<bool bCount>
    void dtrace       ( int, const char *format, /*va_list args*/... );
    void dtrace_repeat( int, int i_times, const char *format, /*va_list args*/... );
    bool update       ( state_type stateval );
    int  getLastError() { return (int)m_error_code;  }
    const char*  getChannelName( int channel_number );
    void getChannelsList( std::string& sChannels );
    std::string getErrMessage();
    int64_t getChannelCounter( int channel ) { return chanRules[channel].getCounter(); }
    void    decrementChannelCounter( int channel ) { chanRules[channel].decrementCounter(); }
};

}
