/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "util/errors.h"
#include <ostream>
#include <iostream>
#include <gsl_errno.h>
#include <cstring>

#ifdef __GNU_LIBRARY__
#include <execinfo.h>
#include <stdlib.h>
#include <cxxabi.h>
#endif

/** Standard exception classes for OpenMalaria. */
namespace OM { namespace util {

const char *Messages::switchDefault = "no valid option selected in a switch statement or if/else tree";

BaseException::BaseException(const string& msg, int code) :
    runtime_error(msg),
    m_errCode(code)
{}

TracedException::TracedException(const string& msg, const char *f, int l, int code, int s) :
    BaseException(msg,code),
    m_file(f), m_line(l), m_start(s)
{
#ifdef OM_GNU_STACK_TRACE
    // http://www.gnu.org/software/libc/manual/html_node/Backtraces.html
// probably 100 is big enough
#define MAX_STACK_SIZE 100
    void *array[MAX_STACK_SIZE];
    
    m_length = backtrace (array, MAX_STACK_SIZE);
    m_trace = backtrace_symbols (array, m_length);
#endif
}
TracedException::~TracedException() throw(){
#ifdef OM_GNU_STACK_TRACE
    free (m_trace);
#endif
}
ostream& operator<<(ostream& stream, const TracedException& e){
    stream<<"Call stack";
    if( e.m_file != 0 )
        stream << ", starting from "<<e.m_file<<':'<<e.m_line;
    stream<<":\n";
#ifdef OM_NO_STACK_TRACE
    stream << "backtrace capture disabled\n";
#elif defined OM_GNU_STACK_TRACE
    // demangle output: http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
    // spec: http://www.ib.cnea.gov.ar/~oop/biblio/libstdc++/namespaceabi.html
    
    //string binary;
    // Usually e.start=1: we skip the first frame (TracedException constructor)
    for( size_t i=e.m_start; i<e.m_length; ++i ){
        char *line = e.m_trace[i];
        char *lb = strchr(line, '(');
        char *plus = 0;
        char *rb = 0;
        char *demangled = 0;
        if( lb != 0 ){
            /* Print out the path to the binary? Seems pointless.
            if( binary != string( line, lb-line ) ){
                binary = string( line, lb-line );
                stream << "in " << binary << ":\n";
            }
            */
            *lb = '\0';
            plus = strchr(lb+1,'+');
            if( plus != 0 ){
                *plus = '\0';
                int status;
                demangled = abi::__cxa_demangle(lb+1, 0, 0, &status);
                /*if( status != 0 ){
                    cerr << "status: "<<status<<"\tstring:"<<lb+1<<endl;
                }*/
                rb = strchr(plus+1, ')');
            }else{
                rb = strchr(lb+1, ')');
            }
            if( rb != 0 ){
                *rb = '\0';
            }
        }
        if( lb == 0 ){
            stream << line;
        }else{
            if( plus != 0 ){
                stream << "+" << (plus+1);
            }
            if( demangled != 0 ){
                stream << '\t' << demangled;
                free(demangled);
            }else if( lb != 0 && rb > lb ){
                stream << '\t' << (lb+1) << "()";
            }
        }
        stream << '\n';
    }
#else
    stream << "sorry, no trace from this platform!\n";
#endif
    return stream;
}

UnimplementedException::UnimplementedException(const string& msg) :
        BaseException(msg, Error::NotImplemented)
{}

const string MsgPre = "Sorry, this feature hasn't been implemented by OpenMalaria yet: ";
vector<char> msgBuf;

const char *UnimplementedException::message() const{
    const char *baseMsg = BaseException::what();
    size_t baseLen = strlen( baseMsg );
    size_t required = baseLen + MsgPre.size();
    if( msgBuf.size() == 0 ){
        msgBuf.assign( required, '0' );
        memcpy( &msgBuf[0], MsgPre.c_str(), MsgPre.size() );
    }else if( msgBuf.size() < required ){
        msgBuf.resize( required );
    }
    memcpy( &msgBuf[MsgPre.size()], baseMsg, baseLen );
    return &msgBuf[0];
}

XmlScenarioError::XmlScenarioError(const string& msg) :
    BaseException(msg, Error::XmlScenario)
{}

FormatError::FormatError(const string& msg) :
    XmlScenarioError(msg)
{}

CheckpointError::CheckpointError(const string& msg) :
    TracedException(msg, 0, 0, Error::Checkpoint)
{}

CmdException::CmdException(const string& msg, int code) :
    BaseException(msg, code)
{}

void gslHandler( const char *reason,
                  const char *file,
                  int line,
                  int gslErrno){
    throw TracedException(reason, file, line, Error::GSL, 2);
}

void setGslHandler() {
    gsl_set_error_handler( &gslHandler );
}
   
} }
