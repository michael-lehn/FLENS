#include <flens/examples/lu/scheduler.h>
#include <cassert>
#include <cstdlib>
#include <fstream>

namespace flens {

Scheduler::Scheduler(int numThreads_)
    : numThreads(numThreads_), numThreadsBusy(0), errorCode(0), stop(false)
{
    std::unique_lock<std::mutex>  lock(mtx);

    for (int i=0; i<numThreads; ++i) {
        auto t = std::thread( [=] { this->runThread(); } );
        t.detach();
    }
}

void
Scheduler::pause()
{
    std::unique_lock<std::mutex>  lock(mtx);

    stop = true;
}

void
Scheduler::resume()
{
    bool jobReady = false;
    {
        std::unique_lock<std::mutex>  lock(mtx);

        stop     = false;
        jobReady = (ready.size()>0);
    }
    if (jobReady) {
        cv_thread.notify_one();
    }
}

void
Scheduler::add(std::string keyStr, Task task,
               std::initializer_list<std::string> pre)
{
    bool jobReady = false;
    {
        std::unique_lock<std::mutex>  lock(mtx);

        Key  key      = getKey(keyStr);

        while (ready.size()>15) {
            cv_scheduler.wait(lock);
        }

        #ifndef NDEBUG
        if (scheduled.count(key)!=0) {
            std::cerr << "Already scheduled: " << keyStr << std::endl;
            assert(0);
        }
        #endif

        scheduled[key]         = task;

        for (auto p_ : pre) {
            auto p = getKey(p_);

            if (scheduled.count(p)>0) {
                successor[p].push_back(key);
                ++predecessorCount[key];
            } else {
                assert(scheduled.count(p)==0);
            }
        }

        if (predecessorCount[key]==0) {
            ready.push_back(key);
        }
        jobReady = (ready.size()>0);
    }
    if (jobReady) {
        cv_thread.notify_one();
    }
}

void
Scheduler::addArc(std::string preStr, std::string succStr)
{
    std::unique_lock<std::mutex>  lock(mtx);

    Key succ = getKey(succStr);
    Key pre  = getKey(preStr);
    // Dependencies for 'succ' must not be added if 'succ' is
    // already scheduled.

    if (scheduled.count(pre)==0) {
        return;
    }

    #ifndef NDEBUG
    if (scheduled.count(succ)!=0) {
        std::cerr << "Already scheduled and possibly running: "
                  << succStr << std::endl;
        std::cerr << preStr << " -> " << succStr << std::endl;
        assert(0);
    }
    #endif

    successor[pre].push_back(succ);
    ++predecessorCount[succ];
}

void
Scheduler::abort(int errorCode_)
{
    assert(errorCode_!=0);
    {
        std::unique_lock<std::mutex> lock(mtx);

        std::cerr << "abort: error code = " << errorCode_ << std::endl;

        errorCode = errorCode_;
        scheduled.clear();
    }
    cv_scheduler.notify_one();
}

void
Scheduler::dumpGraph(const char *filename)
{
    std::unique_lock<std::mutex> lock(mtx);

    std::ofstream out(filename);

    out << "digraph G {\nnode [shape=\"box\"];" << std::endl;
    for (auto t : tag) {
        out << t.first << "[label=\"" << t.second << "\"]"
            << std::endl;
    }
    for (auto s : successor) {
        for (auto succ : s.second) {
            out << "   \"" << s.first
                << "\" -> \"" << succ
                << "\";" << std::endl;
        }
    }
    out << "}" << std::endl;
    out.close();
}

int
Scheduler::join()
{
    std::unique_lock<std::mutex> lock(mtx);

    while (scheduled.size()+numThreadsBusy>0) {
        if (ready.size()+numThreadsBusy==0) {
            std::cerr << "Error: Task queue not empty but no job ready"
                         " or running." << std::endl;
            std::exit(666);
        }
        cv_scheduler.wait(lock);
    }

    // clean up
    ready.clear();
    keyString.clear();
    tag.clear();
    scheduled.clear();
    predecessorCount.clear();
    successor.clear();

    // restore error code and return old value
    int errorCode_ = errorCode;
    errorCode = 0;

    return errorCode_;
}

std::pair<Scheduler::Task, Scheduler::Key>
Scheduler::getTask()
{
    bool jobReady = false;
    Key  key;
    Task task;
    {
        std::unique_lock<std::mutex> lock(mtx);

        while (stop || ready.size()==0) {
            cv_thread.wait(lock);
        }

        key  = ready.front();
        task = scheduled[key];

        ready.pop_front();
        ++numThreadsBusy;

        if (ready.size()>0) {
            jobReady = true;
        }
    }
    if (jobReady) {
        cv_thread.notify_one();
    }
    return std::pair<Task, Key>(task,key);
}

void
Scheduler::runThread()
{
    while(1) {
        std::pair<Task,Key> job = getTask();
        job.first();
        taskDone(job.second);
    }
}

void
Scheduler::taskDone(Key key)
{
    bool allDone = false;

    {
        std::unique_lock<std::mutex> lock(mtx);

        if (errorCode!=0) {
            allDone= true;
        } else {
            for (auto succ : successor[key]) {
                --predecessorCount[succ];
                if (scheduled.count(succ)>0) {
                    if (predecessorCount[succ]==0) {
                        predecessorCount.erase(succ);
                        ready.push_back(succ);
                    }
                }
            }
            scheduled.erase(key);
            successor.erase(key);
        }
        --numThreadsBusy;

        if (scheduled.size()==0 && numThreadsBusy==0) {
            allDone = true;
        }
        if (ready.size()==10) {
            allDone = true;
        }
    }
    if (allDone) {
        cv_scheduler.notify_one();
    }

}

Scheduler::Key
Scheduler::getKey(std::string keyStr)
{
    if (keyString.count(keyStr)==0) {
        keyString[keyStr] = keyString.size();
    }
    tag[keyString[keyStr]] = keyStr;

    return keyString[keyStr];
}

} // namespace flens
