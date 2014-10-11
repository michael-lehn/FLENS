/*
 *   Copyright (c) 2010, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cxxstd/algorithm.h>
#include <cxxstd/iostream.h>
#include <cxxstd/string.h>
#include <flens/debug/auxiliary/closurelog.h>
#include <flens/debug/auxiliary/closurelogstream.h>
#include <flens/auxiliary/macros.h>

namespace flens { namespace verbose {

bool
ClosureLog::_started = false;

int
ClosureLog::_indentLevel = -1;

std::string
ClosureLog::_tag;

std::ofstream
ClosureLog::_out;

VariablePool
ClosureLog::variablePool;

ClosureLogStream
ClosureLog::_closureLogStream = ClosureLogStream(ClosureLog::variablePool,
                                                 ClosureLog::_out);

void
ClosureLog::start(const char *filename, bool clearLog)
{
    std::ios_base::openmode mode = (clearLog) ? std::ios_base::trunc
                                              : std::ios_base::app;
    _out.open(filename, mode | std::ios_base::out);
    _started = true;
    _indentLevel = 0;
    _tag.assign("");

}

void
ClosureLog::stop()
{
    _out << std::endl;
    _out.close();
    _started = false;
}

bool
ClosureLog::started()
{
    return _started;
}

void
ClosureLog::separator()
{
    if (_started && _indentLevel==1) {
        std::string sep(80, '-');
        _out << std::endl << sep << std::endl;
    }
}

bool
ClosureLog::openEntry()
{
    if (_started) {
        return true;
    }
    return false;
}

bool
ClosureLog::createEntry()
{
    if (_started) {
        ++_indentLevel;
        return true;
    }
    return false;
}

void
ClosureLog::closeEntry()
{
    if (_started) {
        --_indentLevel;
    }
}


void
ClosureLog::setTag(const char *tag)
{
    _tag.assign(tag);
}

void
ClosureLog::unsetTag()
{
    _tag.assign("");
}

ClosureLogStream &
ClosureLog::append(bool startNewLine)
{
    ASSERT(_started);


    if (startNewLine) {
        std::string indent(std::max((_indentLevel-1)*4, 0), ' ');
        _out << std::endl;
        _out << indent << _tag;
    }
    return _closureLogStream;
}

} } // namespace verbose, namespace flens
