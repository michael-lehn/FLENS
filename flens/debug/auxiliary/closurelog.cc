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
ClosureLog::started_ = false;

int
ClosureLog::indentLevel_ = -1;

std::string
ClosureLog::tag_;

std::ofstream
ClosureLog::out_;

VariablePool
ClosureLog::variablePool;

ClosureLogStream
ClosureLog::closureLogStream_ = ClosureLogStream(ClosureLog::variablePool,
                                                 ClosureLog::out_);

void
ClosureLog::start(const char *filename, bool clearLog)
{
    std::ios_base::openmode mode = (clearLog) ? std::ios_base::trunc
                                              : std::ios_base::app;
    out_.open(filename, mode | std::ios_base::out);
    started_ = true;
    indentLevel_ = 0;
    tag_.assign("");

}

void
ClosureLog::stop()
{
    out_ << std::endl;
    out_.close();
    started_ = false;
}

bool
ClosureLog::started()
{
    return started_;
}

void
ClosureLog::separator()
{
    if (started_ && indentLevel_==1) {
        std::string sep(80, '-');
        out_ << std::endl << sep << std::endl;
    }
}

bool
ClosureLog::openEntry()
{
    if (started_) {
        return true;
    }
    return false;
}

bool
ClosureLog::createEntry()
{
    if (started_) {
        ++indentLevel_;
        return true;
    }
    return false;
}

void
ClosureLog::closeEntry()
{
    if (started_) {
        --indentLevel_;
    }
}


void
ClosureLog::setTag(const char *tag)
{
    tag_.assign(tag);
}

void
ClosureLog::unsetTag()
{
    tag_.assign("");
}

ClosureLogStream &
ClosureLog::append(bool startNewLine)
{
    ASSERT(started_);


    if (startNewLine) {
        std::string indent(std::max((indentLevel_-1)*4, 0), ' ');
        out_ << std::endl;
        out_ << indent << tag_;
    }
    return closureLogStream_;
}

} } // namespace verbose, namespace flens
